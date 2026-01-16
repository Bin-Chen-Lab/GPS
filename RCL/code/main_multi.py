# -*- coding:utf-8 -*-
import os
import torch 
import torch.nn as nn
import torch.nn.functional as F
from torch.autograd import Variable
import torchvision.transforms as transforms
from model import MLP
import argparse, sys
import numpy as np
import datetime
import shutil
# comment

from loss import loss_multi
from DrugGene import DrugGene
from sklearn.metrics import confusion_matrix, f1_score
import pickle

parser = argparse.ArgumentParser()
parser.add_argument('--lr', type = float, default = 0.001)
parser.add_argument('--result_dir', type = str, help = 'dir to save result txt files', default = 'results/')
parser.add_argument('--forget_rate', type = float, help = 'forget rate', default = None)
parser.add_argument('--num_gradual', type = int, default = 4, help='how many epochs for linear drop rate, can be 5, 10, 15. This parameter is equal to Tk for R(T) in Co-teaching paper.')
parser.add_argument('--exponent', type = float, default = 1, help='exponent of the forget rate, can be 0.5, 1, 2. This parameter is equal to c in Tc for R(T) in Co-teaching paper.')
parser.add_argument('--top_bn', action='store_true')
parser.add_argument('--cl', type = str, help = 'cell line', default='MCF7')
parser.add_argument('--n_epoch', type=int, default=25)
parser.add_argument('--seed', type=int, default=1)
parser.add_argument('--print_freq', type=int, default=200)
parser.add_argument('--num_workers', type=int, default=4, help='how many subprocesses to use for data loading')
parser.add_argument('--epoch_decay_start', type=int, default=40)
parser.add_argument('--fuzzy_exponent', type=int, default=4, help='exponent for the knowledge share fuzzy rate')
parser.add_argument('--stop_fuzzy', type=int, default=2, help='when to stop fuzzy compared to sample decay')
parser.add_argument('--num_networks', type=int, default=3, help='number of network branches')
parser.add_argument('--num_classes', type=int, default=3, help='number of classes')

args = parser.parse_args()

# Seed
def seed_everything(seed):
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)

    os.environ['PYTHONHASHSEED'] = str(seed)
    torch.backends.cudnn.deterministic = True

# Hyper Parameters
batch_size = 128
learning_rate = args.lr
input_size = 2131 #This parameter needs revision
num_classes = args.num_classes
args.top_bn = False


# load dataset
train_dataset = DrugGene(root='../data/output',
                         train=True,
                         num_classes=num_classes)
num_iter_per_epoch = int(len(train_dataset.train_labels)/batch_size)
print('num_iter_per_epoch {:d}'.format(num_iter_per_epoch))
test_dataset = DrugGene(root='../data/output',
                         train=False,
                         num_classes=num_classes)

if args.forget_rate is None:
    forget_rate=0.2
else:
    forget_rate=args.forget_rate

quality = train_dataset.train_quality

# Adjust learning rate and betas for Adam Optimizer
mom1 = 0.9
mom2 = 0.1
alpha_plan = [learning_rate] * args.n_epoch
beta1_plan = [mom1] * args.n_epoch
for i in range(args.epoch_decay_start, args.n_epoch):
    alpha_plan[i] = float(args.n_epoch - i) / (args.n_epoch - args.epoch_decay_start) * learning_rate
    beta1_plan[i] = mom2

def adjust_learning_rate(optimizer, epoch):
    for param_group in optimizer.param_groups:
        param_group['lr']=alpha_plan[epoch]
        param_group['betas']=(beta1_plan[epoch], 0.999) # Only change beta1
        
# define drop rate schedule
rate_schedule = np.ones(args.n_epoch)*forget_rate
print(rate_schedule)
rate_schedule[1:args.num_gradual+1] = np.linspace(0, forget_rate**args.exponent, args.num_gradual)
print(rate_schedule)

knowledge_fuzzy_rate = np.zeros(args.n_epoch)
knowledge_fuzzy_rate[1:args.num_gradual*args.stop_fuzzy+1] = 1 - (np.arange(0, args.num_gradual*args.stop_fuzzy, 1)/float(args.stop_fuzzy*args.num_gradual))**args.fuzzy_exponent
print(knowledge_fuzzy_rate)

save_dir = args.result_dir + args.cl +'/multi/'

if not os.path.exists(save_dir):
    os.system('mkdir -p %s' % save_dir)

model_str = args.cl + '_multi_net_'+str(args.num_networks) + '_forget_rate_' + str(args.forget_rate) + '_stop_fuzzy_' + str(args.stop_fuzzy) +  '_fuzzy_exp_' + str(args.fuzzy_exponent) +  '_seed_' + str(args.seed)

txtfile=save_dir+model_str+".txt"
nowTime=datetime.datetime.now().strftime('%Y-%m-%d-%H:%M:%S')
if os.path.exists(txtfile):
    os.system('mv %s %s' % (txtfile, txtfile+".bak-%s" % nowTime))  # rename exsist file for collison


def accuracy(logit, target, topk=(1,)):
    """Computes the precision@k for the specified values of k"""
    output = F.softmax(logit, dim=1)
    maxk = max(topk)
    batch_size = target.size(0)

    _, pred = output.topk(maxk, 1, True, True)
    pred = pred.t()
    correct = pred.eq(target.view(1, -1).expand_as(pred))

    res = []
    for k in topk:
        correct_k = correct[:k].view(-1).float().sum(0, keepdim=True)
        res.append(correct_k.mul_(100.0 / batch_size))
    return res

# Train the Model
def train(train_loader,epoch, all_models, all_optimizers):
    print 'Training %s...' % model_str
    pure_ratio_list=[]

    train_total=0
    train_correct=0
    for i, (images, labels, indexes) in enumerate(train_loader):
        # print('iteration===={:d}'.format(i))
        # print(indexes)

        ind=indexes.cpu().numpy().transpose()
        if i>num_iter_per_epoch:
            break
      
        images = Variable(images).cuda()
        labels = Variable(labels).cuda()
        
        # Forward + Backward + Optimize
        all_logits = dict()
        for k in range(args.num_networks):
            model_name = 'model' + str(k)
            logits_name = 'logits' + str(k)

            model = all_models[model_name]
            logits = model(images)
            all_logits[logits_name] = logits

            if k==0:
                prec, _ = accuracy(logits, labels, topk=(1, 2))
                train_total+=1
                train_correct+=prec

        all_losses, all_pure_ratios = loss_multi(all_logits, labels, rate_schedule[epoch], ind, quality, knowledge_fuzzy_rate[epoch], i)
        pure_ratio_list.append(100*all_pure_ratios[0])

        for k in range(args.num_networks):
            optimizer_name = 'optimizer' + str(k)
            optimizer = all_optimizers[optimizer_name]
            loss_name = 'loss' + str(k)
            loss = all_losses[loss_name]

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

        if (i+1) % args.print_freq == 0:
            print('Epoch [%d/%d], Iter [%d/%d] Training Accuracy1: %.4F, Loss1: %.4f Pure Ratio1 %.4f'
                  %(epoch, args.n_epoch, i+1, len(train_dataset)//batch_size, prec, all_losses['loss0'].data, np.sum(pure_ratio_list)/len(pure_ratio_list)))
            print(all_pure_ratios)
    train_acc=float(train_correct)/float(train_total)
    return train_acc, pure_ratio_list


def eval_metrics(labels_list, preds_list):
    """list of batch labels and batch preds"""
    # faltten first
    #labels_flatten = [item.data for sublist in labels_list for item in sublist]
    #preds_flatten = [item.data for sublist in preds_list for item in sublist]
    labels_flatten = [item for sublist in labels_list for item in sublist]
    preds_flatten = [item for sublist in preds_list for item in sublist]
    cm = confusion_matrix(labels_flatten, preds_flatten)
    f1 = f1_score(labels_flatten, preds_flatten, average='macro')
    return cm, f1


# Evaluate the Model
def evaluate(test_loader, all_models):
    print 'Evaluating %s...' % model_str

    all_accs = []
    all_f1s = []
    for k in range(args.num_networks):
        model_name = 'model' + str(k)
        model = all_models[model_name]

        model.eval()
        correct = 0
        total = 0
        labels_list = []
        pred_list = []
        output_list = []
        for images, labels, _ in test_loader:
            images = Variable(images).cuda()
            logits = model(images)
            outputs = F.softmax(logits, dim=1)
            _, pred = torch.max(outputs.data, 1)
            total += labels.size(0)
            correct += (pred.cpu() == labels).sum()

            pred_list.append(pred.cpu().numpy())
            labels_list.append(labels.cpu().numpy())
            output_list.append(outputs.detach().cpu().numpy())

        acc = 100*float(correct)/float(total)
        cm, f1 = eval_metrics(labels_list, pred_list)
        all_accs.append(acc)
        all_f1s.append(f1)
        if k==0:
            print(cm)
            print('f1_score {:.2f}'.format(f1))
            pred_final = pred_list
            labels_final = labels_list
            outputs_final = output_list
            res = [labels_final, pred_final, outputs_final]

    return all_accs, all_f1s, res


def main():
    seed_everything(args.seed)
    # Data Loader (Input Pipeline)
    print 'loading dataset...'
    train_loader = torch.utils.data.DataLoader(dataset=train_dataset,
                                               batch_size=batch_size, 
                                               num_workers=args.num_workers,
                                               drop_last=True,
                                               shuffle=True)
    
    test_loader = torch.utils.data.DataLoader(dataset=test_dataset,
                                              batch_size=batch_size, 
                                              num_workers=args.num_workers,
                                              drop_last=False,
                                              shuffle=False)
    # Define models
    print 'building model...'
    all_models = dict()
    all_optimizers = dict()
    for k in range(args.num_networks):
        model_name = 'model' + str(k)
        optimizer_name = 'optimizer' + str(k)

        model = MLP(input_size=input_size, n_outputs=num_classes)
        if torch.cuda.device_count() > 1:
            model = nn.DataParallel(model)
        model.cuda()
        optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
        # if k==0:
        print(model.parameters)

        all_models[model_name] = model
        all_optimizers[optimizer_name] = optimizer

    with open(txtfile, "a") as myfile:
        myfile.write('epoch: train_acc1 test_acc1 test_f1 pure_ratio1\n')

    epoch=0
    train_acc1=0
    mean_pure_ratio1=0
    # evaluate models with random weights
    test_accs, test_f1s, _=evaluate(test_loader, all_models)
    print('all test accs')
    print(test_accs)
    print('Epoch [%d/%d] Test Accuracy on the %s test images: Model1 %.4f %% F1 %.4f Pure ratio %.4f %% ' % (epoch+1, args.n_epoch, len(test_dataset), test_accs[0], test_f1s[0], mean_pure_ratio1))
    # save results
    with open(txtfile, "a") as myfile:
        myfile.write(str(int(epoch)) + ': '  + str(train_acc1)  +' ' + str(test_accs[0]) + " " + str(test_f1s[0]) + " " + str(mean_pure_ratio1) + "\n")

    # training
    for epoch in range(1, args.n_epoch):
        # train models
        for k in range(args.num_networks):
            model_name = 'model' + str(k)
            optimizer_name = 'optimizer' + str(k)

            model = all_models[model_name]
            optimizer = all_optimizers[optimizer_name]

            model.train()
            adjust_learning_rate(optimizer, epoch)

        train_acc1, pure_ratio_1_list = train(train_loader, epoch, all_models, all_optimizers)
        # evaluate models
        test_accs, test_f1s, res=evaluate(test_loader, all_models)
        print('all test accs')
        print(test_accs)
        mean_pure_ratio1 = sum(pure_ratio_1_list)/len(pure_ratio_1_list)
        # save results
        print('Epoch [%d/%d] Test Accuracy on the %s test images: Model1 %.4f %% F1 %.4f Pure ratio %.4f %% ' % (epoch+1, args.n_epoch, len(test_dataset), test_accs[0], test_f1s[0], mean_pure_ratio1))
        with open(txtfile, "a") as myfile:
            myfile.write(str(int(epoch)) + ': '  + str(train_acc1) +' '  + str(test_accs[0]) + " " + str(test_f1s[0]) + " " + str(mean_pure_ratio1) + "\n")

    torch.save(all_models, save_dir+model_str+'_model.pkl')
        # with open(save_dir+'res.pkl', 'wb') as f:
        #     pickle.dump(res, f)


if __name__=='__main__':
    main()

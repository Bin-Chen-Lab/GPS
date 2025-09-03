import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.autograd import Variable
import numpy as np



# Loss functions
def loss_coteaching(y_1, y_2, t, forget_rate, ind, quality):
    include_or_not_1 = np.zeros([len(y_1), 1], dtype=np.int8)
    include_or_not_2 = np.zeros([len(y_2), 1], dtype=np.int8)

    loss_1 = F.cross_entropy(y_1, t, reduce=False)
    ind_1_sorted = np.argsort(loss_1.cpu().data).cuda()  # np.argsort(loss_1.data).cuda()
    loss_1_sorted = loss_1[ind_1_sorted]

    loss_2 = F.cross_entropy(y_2, t, reduce=False)
    ind_2_sorted = np.argsort(loss_2.cpu().data).cuda()

    remember_rate = 1 - forget_rate
    num_remember = int(remember_rate * len(loss_1_sorted))
    # print('remember rate {:.3f} remember # {:d}'.format(remember_rate, num_remember))

    pure_ratio_1 = np.sum(quality[ind[ind_1_sorted.cpu()[:num_remember]]]) / float(num_remember)  # cpu
    pure_ratio_2 = np.sum(quality[ind[ind_2_sorted.cpu()[:num_remember]]]) / float(num_remember)

    ind_1_update = ind_1_sorted[:num_remember]
    ind_2_update = ind_2_sorted[:num_remember]

    # print(np.concatenate([ind_1_update.cpu().reshape(-1, 1), ind_2_update.cpu().reshape(-1, 1)], axis=1))
    # agree_n = len(np.intersect1d(ind_1_update.cpu(), ind_2_update.cpu()))
    # agree_p = len(np.intersect1d(ind_1_update.cpu(), ind_2_update.cpu())) / float(num_remember)
    # print('agree n={:d} p={:.3f}'.format(agree_n, agree_p))

    include_or_not_1[ind_1_update.cpu()] = 1
    include_or_not_2[ind_2_update.cpu()] = 1

    # exchange
    loss_1_update = F.cross_entropy(y_1[ind_2_update], t[ind_2_update])
    loss_2_update = F.cross_entropy(y_2[ind_1_update], t[ind_1_update])

    return torch.sum(loss_1_update) / num_remember, torch.sum(
        loss_2_update) / num_remember, pure_ratio_1, pure_ratio_2, include_or_not_1, include_or_not_2


def loss_decouple(y_1, y_2, t, forget_rate, ind, quality):
    
    outputs1 = F.softmax(y_1, dim=1)
    _, pred1 = torch.max(outputs1.data, 1)
    
    outputs2 = F.softmax(y_2, dim=1)
    _, pred2 = torch.max(outputs2.data, 1)
    
    update_index = np.where(pred1.cpu().numpy() != pred2.cpu().numpy())[0]
    remember_rate = 1 - forget_rate
    num_remember = int(remember_rate * len(t))
    # print('remember rate {:.3f} num_remember {:d} batch_size {:d}'.format(remember_rate, num_remember, len(loss_1_sorted)))
    
    pure_ratio_1 = np.sum(quality[ind[update_index]]) / float(len(update_index))  # cpu
    pure_ratio_2 = np.sum(quality[ind[update_index]]) / float(len(update_index))
    
    # exchange
    loss_1_update = F.cross_entropy(y_1[update_index], t[update_index])
    loss_2_update = F.cross_entropy(y_2[update_index], t[update_index])
    
    return torch.sum(loss_1_update) / len(update_index), torch.sum(loss_2_update) / len(update_index), pure_ratio_1, pure_ratio_2


def loss_spl(y, t, forget_rate, ind, quality):
    include_or_not = np.zeros([len(y), 1], dtype=np.int8)

    loss = F.cross_entropy(y, t, reduce=False)
    ind_sorted = np.argsort(loss.cpu().data).cuda()  # np.argsort(loss_1.data).cuda()
    loss_sorted = loss[ind_sorted]

    remember_rate = 1 - forget_rate
    num_remember = int(remember_rate * len(loss_sorted))
    # print('remember rate {:.3f} num_remember {:d} batch_size {:d}'.format(remember_rate, num_remember, len(loss_sorted)))

    pure_ratio = np.sum(quality[ind[ind_sorted.cpu()[:num_remember]]]) / float(num_remember)  # cpu
    ind_update = ind_sorted[:num_remember]
    loss_update = F.cross_entropy(y[ind_update], t[ind_update])

    include_or_not[ind_update.cpu()] = 1
    return torch.sum(loss_update) / num_remember, pure_ratio, include_or_not


def loss_triteaching(y_1, y_2, y_3, t, forget_rate, ind, noise_or_not):
    loss_1 = F.cross_entropy(y_1, t, reduce=False)
    ind_1_sorted = np.argsort(loss_1.cpu().data).cuda()  # np.argsort(loss_1.data).cuda()
    loss_1_sorted = loss_1[ind_1_sorted]

    loss_2 = F.cross_entropy(y_2, t, reduce=False)
    ind_2_sorted = np.argsort(loss_2.cpu().data).cuda()
    loss_2_sorted = loss_2[ind_2_sorted]

    loss_3 = F.cross_entropy(y_3, t, reduce=False)
    ind_3_sorted = np.argsort(loss_3.cpu().data).cuda()
    loss_3_sorted = loss_3[ind_3_sorted]

    remember_rate = 1 - forget_rate
    num_remember = int(remember_rate * len(loss_1_sorted))
    # print('remember rate {:.3f} remember # {:d}'.format(remember_rate, num_remember))

    pure_ratio_1 = np.sum(noise_or_not[ind[ind_1_sorted.cpu()[:num_remember]]]) / float(num_remember)  # cpu
    pure_ratio_2 = np.sum(noise_or_not[ind[ind_2_sorted.cpu()[:num_remember]]]) / float(num_remember)
    pure_ratio_3 = np.sum(noise_or_not[ind[ind_3_sorted.cpu()[:num_remember]]]) / float(num_remember)

    ind_1_update = ind_1_sorted[:num_remember]
    ind_2_update = ind_2_sorted[:num_remember]
    ind_3_update = ind_3_sorted[:num_remember]

    # print(np.concatenate([ind_1_update.cpu().reshape(-1, 1), ind_2_update.cpu().reshape(-1, 1)], axis=1))
    # agree_n = len(np.intersect1d(ind_1_update.cpu(), ind_2_update.cpu()))
    # agree_p = len(np.intersect1d(ind_1_update.cpu(), ind_2_update.cpu())) / float(num_remember)
    # print('agree n={:d} p={:.3f}'.format(agree_n, agree_p))

    # exchange
    loss_1_update = F.cross_entropy(y_1[ind_3_update], t[ind_3_update])
    loss_2_update = F.cross_entropy(y_2[ind_1_update], t[ind_1_update])
    loss_3_update = F.cross_entropy(y_3[ind_2_update], t[ind_2_update])

    return torch.sum(loss_1_update) / num_remember, torch.sum(loss_2_update) / num_remember, torch.sum(
        loss_3_update) / num_remember, pure_ratio_1, pure_ratio_2, pure_ratio_3


def update_index(ind_union, ind_intersect, fuzzy_rate, seed=1):
    np.random.seed(seed)
    n_intersect = len(ind_intersect)
    n_union = len(ind_union)

    if n_intersect == n_union:
        ind_update = ind_intersect
    else:
        n_in = int(np.floor((n_union-n_intersect)*fuzzy_rate))
        ind_diff = np.setdiff1d(ind_union, ind_intersect)
        ind_in = np.random.choice(ind_diff, n_in)
        ind_update = np.concatenate([ind_intersect, ind_in])
    return ind_update


def loss_tri_co(y_1, y_2, y_3, t, forget_rate, ind, fuzzy_rate):
    loss_1 = F.cross_entropy(y_1, t, reduce=False)
    ind_1_sorted = np.argsort(loss_1.cpu().data).cuda()  # np.argsort(loss_1.data).cuda()
    loss_1_sorted = loss_1[ind_1_sorted]

    loss_2 = F.cross_entropy(y_2, t, reduce=False)
    ind_2_sorted = np.argsort(loss_2.cpu().data).cuda()
    loss_2_sorted = loss_2[ind_2_sorted]

    loss_3 = F.cross_entropy(y_3, t, reduce=False)
    ind_3_sorted = np.argsort(loss_3.cpu().data).cuda()
    loss_3_sorted = loss_3[ind_3_sorted]

    remember_rate = 1 - forget_rate
    num_remember = int(remember_rate * len(loss_1_sorted))
    print('remember rate {:.3f} remember # {:d}'.format(remember_rate, num_remember))

    ind_1_update = ind_1_sorted[:num_remember]
    ind_2_update = ind_2_sorted[:num_remember]
    ind_3_update = ind_3_sorted[:num_remember]

    ind_1_union = np.union1d(ind_2_update.cpu(), ind_3_update.cpu())
    ind_2_union = np.union1d(ind_1_update.cpu(), ind_3_update.cpu())
    ind_3_union = np.union1d(ind_1_update.cpu(), ind_2_update.cpu())

    ind_1_intersect = np.intersect1d(ind_2_update.cpu(), ind_3_update.cpu())
    ind_2_intersect = np.intersect1d(ind_1_update.cpu(), ind_3_update.cpu())
    ind_3_intersect = np.intersect1d(ind_1_update.cpu(), ind_2_update.cpu())

    ind_1_ensemble = update_index(ind_1_union, ind_1_intersect, fuzzy_rate)
    ind_2_ensemble = update_index(ind_2_union, ind_2_intersect, fuzzy_rate)
    ind_3_ensemble = update_index(ind_3_union, ind_3_intersect, fuzzy_rate)

    n_agree_1 = len(ind_1_ensemble)
    n_agree_2 = len(ind_2_ensemble)
    n_agree_3 = len(ind_3_ensemble)

    print('n_included for 1 2 3 {:d} {:d} {:d}'.format(n_agree_1, n_agree_2, n_agree_3))

    # exchange
    loss_1_update = F.cross_entropy(y_1[ind_1_ensemble], t[ind_1_ensemble])
    loss_2_update = F.cross_entropy(y_2[ind_2_ensemble], t[ind_2_ensemble])
    loss_3_update = F.cross_entropy(y_3[ind_3_ensemble], t[ind_3_ensemble])

    return torch.sum(loss_1_update) / n_agree_1, torch.sum(loss_2_update) / n_agree_2, torch.sum(
        loss_3_update) / n_agree_3


def loss_multi(all_logits, t, forget_rate, ind, quality, fuzzy_rate, iter):

    n_nets = len(all_logits)
    remember_rate = 1 - forget_rate
    num_remember = int(remember_rate * len(t))

    all_inds = dict()
    all_pure_ratios = np.zeros(n_nets)
    for i in range(n_nets):
        logits_name = 'logits' + str(i)
        logits = all_logits[logits_name]

        loss = F.cross_entropy(logits, t, reduce=False)
        ind_sorted = np.argsort(loss.cpu().data).cuda()

        all_pure_ratios[i] = np.sum(quality[ind[ind_sorted.cpu()[:num_remember]]]) / float(num_remember)  # cpu

        ind_name = 'ind' + str(i)
        all_inds[ind_name] = ind_sorted[:num_remember]

    all_losses_update = dict()  # for storing loss as torch varibales
    all_sizes = []
    for i in range(n_nets):
        loss_name = 'loss' + str(i)
        # current ind_update from all the rest inds
        ind_rest = []
        for k in range(n_nets):
            if k != i:
                ind_rest.append(all_inds['ind'+str(k)].cpu().data)
        # print(ind_rest)

        ind_union = reduce(np.union1d, ind_rest)
        ind_intersect = reduce(np.intersect1d, ind_rest)
        ind_ensemble = update_index(ind_union, ind_intersect, fuzzy_rate)

        all_sizes.append(len(ind_ensemble))
        all_losses_update[loss_name] = F.cross_entropy(all_logits['logits'+str(i)][ind_ensemble], t[ind_ensemble])/len(ind_ensemble)

    if iter == 0:
        print('remember rate {:.3f} remember # {:d} include # {:d}'.format(remember_rate, num_remember, all_sizes[0]))
    return all_losses_update, all_pure_ratios



def loss_tri_co_v2(y_1, y_2, y_3, t, forget_rate, ind, noise_or_not):
    loss_1 = F.cross_entropy(y_1, t, reduce=False)
    ind_1_sorted = np.argsort(loss_1.cpu().data).cuda()  # np.argsort(loss_1.data).cuda()
    loss_1_sorted = loss_1[ind_1_sorted]

    loss_2 = F.cross_entropy(y_2, t, reduce=False)
    ind_2_sorted = np.argsort(loss_2.cpu().data).cuda()
    loss_2_sorted = loss_2[ind_2_sorted]

    loss_3 = F.cross_entropy(y_3, t, reduce=False)
    ind_3_sorted = np.argsort(loss_3.cpu().data).cuda()
    loss_3_sorted = loss_3[ind_3_sorted]

    remember_rate = 1 - forget_rate
    num_remember = int(remember_rate * len(loss_1_sorted))
    print('remember rate {:.3f} remember # {:d}'.format(remember_rate, num_remember))

    pure_ratio_1 = np.sum(noise_or_not[ind[ind_1_sorted.cpu()[:num_remember]]]) / float(num_remember)  # cpu
    pure_ratio_2 = np.sum(noise_or_not[ind[ind_2_sorted.cpu()[:num_remember]]]) / float(num_remember)
    pure_ratio_3 = np.sum(noise_or_not[ind[ind_3_sorted.cpu()[:num_remember]]]) / float(num_remember)

    ind_1_forward = ind_1_sorted[:num_remember]
    ind_2_forward = ind_2_sorted[:num_remember]
    ind_3_forward = ind_3_sorted[:num_remember]

    ind_1_ensemble_forward = np.intersect1d(ind_2_forward.cpu(), ind_3_forward.cpu())
    ind_2_ensemble_forward = np.intersect1d(ind_1_forward.cpu(), ind_3_forward.cpu())
    ind_3_ensemble_forward = np.intersect1d(ind_1_forward.cpu(), ind_2_forward.cpu())

    n_agree_1 = len(ind_1_ensemble)
    n_agree_2 = len(ind_2_ensemble)
    n_agree_3 = len(ind_3_ensemble)

    print('n_included for 1 2 3 {:d} {:d} {:d}'.format(n_agree_1, n_agree_2, n_agree_3))

    # exchange
    loss_1_update = F.cross_entropy(y_1[ind_1_ensemble], t[ind_1_ensemble])
    loss_2_update = F.cross_entropy(y_2[ind_2_ensemble], t[ind_2_ensemble])
    loss_3_update = F.cross_entropy(y_3[ind_3_ensemble], t[ind_3_ensemble])

    return torch.sum(loss_1_update) / num_remember, torch.sum(loss_2_update) / num_remember, torch.sum(
        loss_3_update) / num_remember, pure_ratio_1, pure_ratio_2, pure_ratio_3



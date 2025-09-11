from __future__ import print_function
import os
import os.path
import numpy as np
import sys
if sys.version_info[0] == 2:
    import cPickle as pickle
else:
    import pickle

import torch.utils.data as data
from utils import get_drug_fp_batch
import pandas as pd
import random


class DrugGene(data.Dataset):
    """`CIFAR10 <https://www.cs.toronto.edu/~kriz/cifar.html>`_ Dataset.

    Args:
        root (string): Root directory of dataset where directory
            ``cifar-10-batches-py`` exists or will be saved to if download is set to True.
        train (bool, optional): If True, creates dataset from training set, otherwise
            creates from test set.
        transform (callable, optional): A function/transform that  takes in an PIL image
            and returns a transformed version. E.g, ``transforms.RandomCrop``
        target_transform (callable, optional): A function/transform that takes in the
            target and transforms it.

    """

    # define train(test)_data, train(test)_label
    # define __getitem__ to return typical data point you want.
    def __init__(self, root, num_classes=3, train=True, random_seed=0):    #chage train=True to false for debug
        self.train = train  # training set or test set
        self.num_classes = num_classes
        self.random_seed=random_seed

        # now load the picked numpy arrays
        if self.train:
            fn = 'data_train.csv'
            table_train = pd.read_csv(fn)
            table_train["label"] = table_train["label"].astype("int64")
            print(table_train.shape)

            train_labels = table_train['label']
            train_quality = table_train['quality']

            if self.num_classes == 3:
                idx_in = self.down_sampling(train_labels)

                train_smiles = table_train['smiles'][idx_in]
                train_genes = table_train['gene'][idx_in]
                train_labels = np.asarray(train_labels[idx_in]) # be careful, label index need to be reset using np.array
                train_quality = np.asarray(train_quality[idx_in])

            print("get drug features")
            train_smiles_feature = get_drug_fp_batch(train_smiles).astype(np.float32)
            print("get gene features")
            train_genes_feature = self.get_gene_ft_batch(train_genes).astype(np.float32)
            train_data = np.concatenate([train_smiles_feature, train_genes_feature], axis=1)

            self.train_data, self.train_labels, self.train_quality = train_data, train_labels, train_quality

            unique, counts = np.unique(self.train_labels, return_counts=True)
            print(counts)

            print('train_data shape')
            print(self.train_data.shape)
            print(self.train_labels.shape)

        else:
            fn = 'data_test.csv'
            table_test = pd.read_csv(fn)
            table_test["label"] = table_test["label"].astype("int64")
            print(table_test.shape)
            test_smiles = table_test['smiles']
            test_genes = table_test['gene']
            test_labels = table_test['label']
            test_quality = table_test['quality']

            print("get drug features")
            test_smiles_feature = get_drug_fp_batch(test_smiles).astype(np.float32)
            print("get gene features")
            test_genes_feature = self.get_gene_ft_batch(test_genes).astype(np.float32)
            test_data = np.concatenate([test_smiles_feature, test_genes_feature], axis=1)

            self.test_data, self.test_labels, self.test_quality = test_data, test_labels, test_quality

            print('test_data shape')
            print(self.test_data.shape)
            print(self.test_labels.shape)


    def __getitem__(self, index):
        """
        Args:
            index (int): Index

        Returns:
            tuple: (image, target) where target is index of the target class.
        """
        if self.train:
            img, target = self.train_data[index], self.train_labels[index]
        else:
            img, target = self.test_data[index], self.test_labels[index]

        return img, target, index

    def __len__(self):
        if self.train:
            return len(self.train_data)
        else:
            return len(self.test_data)

    def down_sampling(self, y):
        unique, counts = np.unique(y, return_counts=True)
        max_idx = np.argmax(counts)
        max_value = unique[max_idx]
        max_counts = counts[max_idx]
        n_select = np.int((np.sum(counts)-max_counts)*0.5)
        print('max_value, max_counts, n_select')
        print(max_value, max_counts, n_select)

        random.seed(self.random_seed)
        idx_select = random.sample(np.where(y==max_value)[0], k=n_select)
        idx_final = np.concatenate([np.where(y==0)[0], idx_select, np.where(y==2)[0]])

        return idx_final

    def get_gene_ft_batch(self, gene):

        fn = '../data/output/go_fingerprints_l4.csv'
        gene_map = pd.read_csv(fn)
        gene_name = gene_map['gene']
        gene_map = gene_map.drop(columns='gene', axis=1)
        gene_map = gene_map.to_numpy()

        gene_features = []
        for g in gene:
            idx = np.where(gene_name==g)[0][0]
            gene_features.append(gene_map[idx])

        gene_features = np.array(gene_features)
        print(gene_features.shape)
        return gene_features







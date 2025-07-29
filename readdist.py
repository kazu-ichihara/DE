#!/usr/bin/env python
import os
import sys
import pandas as pd
import numpy as np
from tqdm import tqdm
import pysam
import pickle

filename = sys.argv[1]
outpath = sys.argv[2]
with open('./annotation/coding_exon_dict.pickle', 'rb') as f:
    exon_dict = pickle.load(f)

nascent = pd.read_csv('./gencode_DEmotif_nascent_window.csv', index_col=0)
nascent = nascent[nascent['start']>=150]
short = list(nascent[nascent['short']].index)
long = list(nascent[nascent['long']].index)

id_list = list(exon_dict.keys())
bamfile = pysam.AlignmentFile(filename, 'rb')

names = ['short', 'long']
lists = [short, long]

for n in [0,1]:
    key = names[n]
    meta = np.zeros(400)
    id_list = lists[n]
    for item in tqdm(range(len(id_list))):
        name = id_list[item]
        ex = exon_dict[name].values
        chrom = ex[0][1]
        strand = ex[0][5]
        exon = []
        m = nascent.loc[name]['start']
        left = m-150
        right = m+250
        if strand == '+':
            for e in ex:
                exon.append(np.arange(e[3], e[4] + 1))
            exon = np.concatenate(exon)
            exon_meta = exon[left:right]
            pos_to_index = {pos: i for i, pos in enumerate(exon_meta)}
            left_pos = exon_meta[0]
            right_pos = exon_meta[-1]
            if len(exon_meta) == 400:
                for read in bamfile.fetch(chrom, left_pos, right_pos):
                    ps = np.array(read.get_reference_positions())
                    ps = ps[np.isin(ps, exon_meta)]
                    indices = [pos_to_index[p] for p in ps if p in pos_to_index]
                    np.add.at(meta, indices, 1)
        else:
            for e in ex:
                exon.append(np.arange(e[3], e[4] + 1)[::-1])
            exon = np.concatenate(exon)
            exon_meta = exon[left:right]
            pos_to_index = {pos: i for i, pos in enumerate(exon_meta)}
            left_pos = exon_meta[0]
            right_pos = exon_meta[-1]
            if len(exon_meta) == 400:
                for read in bamfile.fetch(chrom, right_pos, left_pos):
                    ps = np.array(read.get_reference_positions())
                    ps = ps[np.isin(ps, exon_meta)]
                    indices = [pos_to_index[p] for p in ps if p in pos_to_index]
                    np.add.at(meta, indices, 1)
    np.save(outpath + '_read_' + key, meta)

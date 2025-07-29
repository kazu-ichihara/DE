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
with open('./annotation/coding_start_dict.pickle', 'rb') as f:
    start_dict = pickle.load(f)
with open('./annotation/coding_stop_dict.pickle', 'rb') as f:
    stop_dict = pickle.load(f)
with open('./annotation/length_dict.pickle', 'rb') as f:
    len_dict = pickle.load(f)

id_list = list(len_dict.keys())
accept = [x for x in id_list if (start_dict[x]>51)&(len_dict[x]-stop_dict[x]>50)]
bamfile = pysam.AlignmentFile(filename, 'rb')

res = {}
for item in tqdm(range(len(accept))):
    name = accept[item]
    ex = exon_dict[name].values
    chrom = ex[0][1]
    strand = ex[0][5]
    exon = []
    start = start_dict[name]
    stop = stop_dict[name]
    cs = np.zeros(5)
    if strand == '+':
        for e in ex:
            exon.append(np.arange(e[3], e[4] + 1))
        exon = np.concatenate(exon)
        start_exon = exon[start]
        stop_exon = exon[stop]
        start_left = exon[start-15]
        start_right = exon[start+15]
        stop_left = exon[stop-15]
        stop_right = exon[stop+15]
        for read in bamfile.fetch(chrom, exon[0], start_left):
            cs[0] += 1
        for read in bamfile.fetch(chrom, start_left, start_right):
            cs[1] += 1
        for read in bamfile.fetch(chrom, start_right, stop_left):
            cs[2] += 1
        for read in bamfile.fetch(chrom, stop_left, stop_right):
            cs[3] += 1
        for read in bamfile.fetch(chrom, stop_right, exon[-1]):
            cs[4] += 1
    else:
        for e in ex:
            exon.append(np.arange(e[3], e[4] + 1)[::-1])
        exon = np.concatenate(exon)
        start_exon = exon[start]
        stop_exon = exon[stop]
        start_left = exon[start-15]
        start_right = exon[start+15]
        stop_left = exon[stop-15]
        stop_right = exon[stop+15]
        for read in bamfile.fetch(chrom, start_left, exon[0]):
            cs[0] += 1
        for read in bamfile.fetch(chrom, start_right, start_left):
            cs[1] += 1
        for read in bamfile.fetch(chrom, stop_left, start_right):
            cs[2] += 1
        for read in bamfile.fetch(chrom, stop_right, stop_left):
            cs[3] += 1
        for read in bamfile.fetch(chrom, exon[-1], stop_right):
            cs[4] += 1

    res[name] = cs

df = pd.DataFrame(res).T
df.columns = ['5\'UTR','Start', 'CDS', 'Stop','3\'UTR']
df.to_csv(outpath+'_regioncount.csv')

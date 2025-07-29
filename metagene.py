#!/usr/bin/env python
import os
import sys
import pandas as pd
import numpy as np
from tqdm import tqdm
import pysam
import pickle
import cv2
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

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
accept = [x for x in id_list if (start_dict[x]>=100)&(len_dict[x]-stop_dict[x]>=100)]
bamfile = pysam.AlignmentFile(filename, 'rb')

meta_utr5 = np.zeros((1000,1))
meta_cds = np.zeros((1000,1))
meta_utr3 = np.zeros((1000,1))
for item in tqdm(range(len(list(accept)))):
    name = accept[item]
    start = start_dict[name]
    stop = stop_dict[name]
    ex = exon_dict[name].values
    chrom = ex[0][1]
    strand = ex[0][5]
    exon = []
    if strand == '+':
        for e in ex:
            exon.append(np.arange(e[3], e[4] + 1))
    else:
        for e in ex:
            exon.append(np.arange(e[3], e[4] + 1)[::-1])
    exon = np.concatenate(exon)
    pos_to_index = {pos: i for i, pos in enumerate(exon)}
    meta = np.zeros(len(exon) + 1)
    fetch_start = int(min(exon))
    fetch_end = int(max(exon))

    for read in bamfile.fetch(chrom, fetch_start, fetch_end):
        ps = np.array(read.get_reference_positions()) + 1
        ps = ps[np.isin(ps, exon)]
        indices = [pos_to_index[p] for p in ps if p in pos_to_index]
        np.add.at(meta, indices, 1)

    utr5 = cv2.resize(meta[:start], (1, 1000))
    cds = cv2.resize(meta[start:stop], (1, 1000))
    utr3 = cv2.resize(meta[stop:], (1, 1000))
    meta_utr5 += utr5
    meta_cds += cds
    meta_utr3 += utr3

conc = np.ravel(np.concatenate([meta_utr5, meta_cds, meta_utr3]))
np.save(outpath+'_metagene', conc)

pp = PdfPages(outpath+'_metagene.pdf')
plt.figure(figsize=(10, 5))
x = np.linspace(0, 3, 3000)
y = conc
plt.plot(x, y, c='black')
plt.xlim(0, 3)
plt.ylim(0,)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.xticks(np.arange(0, 4, 1))
plt.xticks(color="None")
pp.savefig()
pp.close()

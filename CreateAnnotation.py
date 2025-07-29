#!/usr/bin/env python
import os
import pandas as pd
import numpy as np
from Bio import SeqIO
import pickle

os.makedirs('./annotation', exist_ok=True)

gtf = pd.read_csv('./gencode.v31.annotation.gtf', sep='\t', header=None, comment='#')

gtf_transcript = gtf[gtf[2] != 'gene']
gtf_transcript = gtf_transcript.drop([1,5,7], axis=1)
gtf_transcript.columns = ['chr','type','start','end','strand','annotation']
gtf_transcript['length'] = gtf_transcript['end']-gtf_transcript['start']
gtf_transcript['id'] = [x.split(';')[1].split('"')[1] for x in gtf_transcript['annotation']]
gtf_exon = gtf_transcript[gtf_transcript['type'] == 'exon']
gtf_codon = gtf_transcript[(gtf_transcript['type'] == 'start_codon')|(gtf_transcript['type'] == 'stop_codon')]
codons = gtf_codon.groupby('id')['type'].sum()
transcripts = list(codons[(codons.str.contains('start'))&(codons.str.contains('stop'))].index)
gtf_coding = gtf_transcript[gtf_transcript['id'].isin(transcripts)]
gtf_coding['symbol'] = [x.split(';')[3].split('"')[1] for x in gtf_coding['annotation']]
symbol = gtf_coding[['id','symbol']].drop_duplicates('id').set_index('id')

gtf_exon_coding = gtf_exon[gtf_exon['id'].isin(transcripts)]
exon_dict = {name:group.reset_index() for name, group in gtf_exon_coding.groupby('id')}
with open('./annotation/coding_exon_dict.pickle', 'wb') as f:
    pickle.dump(exon_dict, f)
exon_rel_dict = {k:np.cumsum(np.array(v['length'])) for k, v in exon_dict.items()}
with open('./annotation/coding_exon_rel_dict.pickle', 'wb') as f:
    pickle.dump(exon_dict, f)

gtf_start = gtf_coding[gtf_coding['type'] == 'start_codon']
start_dict = {}
for name in range(len(transcripts)):
    starts = gtf_start[gtf_start['id'] == name]
    strand = starts['strand'].values[0]
    ex = exon_dict[name]
    if strand == '+':
        start = starts['start'].values[0]
        start_exon = ex[(ex['start'] <= start) & (ex['end'] >= start)]
        if start_exon.index == 0:
            start_dict[name] = start - start_exon['start'].values[0]
        else:
            start_dict[name] = start - start_exon['start'].values[0] + exon_rel_dict[name][start_exon.index[0]-1] + start_exon.index[0]
    else:
        start = starts['end'].values[0]
        start_exon = ex[(ex['start'] <= start) & (ex['end'] >= start)]
        if start_exon.index == 0:
            start_dict[name] = start_exon['end'].values[0] - start
        else:
            start_dict[name] = start_exon['end'].values[0] - start + exon_rel_dict[name][start_exon.index[0]-1] + start_exon.index[0]
with open('./annotation/coding_start_dict.pickle', 'wb') as f:
    pickle.dump(start_dict, f)

gtf_stop = gtf_coding[gtf_coding['type'] == 'stop_codon']
stop_dict = {}
for name in range(len(transcripts)):
    stops = gtf_stop[gtf_stop['id'] == name]
    strand = stops['strand'].values[0]
    ex = exon_dict[name]
    if strand == '+':
        stop = stops['end'].values[-1]
        stop_exon = ex[(ex['start'] <= stop) & (ex['end'] >= stop)]
        if stop_exon.index == 0:
            stop_dict[name] = stop - stop_exon['start'].values[0] + 1
        else:
            stop_dict[name] = stop - stop_exon['start'].values[0] + exon_rel_dict[name][stop_exon.index[0]-1] + stop_exon.index[0] + 1
    else:
        stop = stops['start'].values[0]
        stop_exon = exon_dict[name][(ex['start'] <= stop) & (ex['end'] >= stop)]
        if stop_exon.index == 0:
            stop_dict[name] = stop_exon['end'].values[0] - stop + 1
        else:
            stop_dict[name] = stop_exon['end'].values[0] - stop + exon_rel_dict[name][stop_exon.index[0]-1] + stop_exon.index[0] + 1
with open('./annotation/coding_stop_dict.pickle', 'wb') as f:
    pickle.dump(stop_dict, f)

fasta_in = './gencode.v31.transcripts.fa'
seq_dict = {}
for record in SeqIO.parse(fasta_in, 'fasta'):
    id_part = record.id
    name = id_part.split('|')[0]
    seq = record.seq
    if name in transcripts:
        seq_dict[name] = seq
len_dict = {k:len(v) for k, v in seq_dict.items()}
with open('./annotation/length_dict.pickle', 'wb') as f:
    pickle.dump(len_dict, f)

def getSlidingWindow(a, size, step=1):
    if not ((type(size) == type(0)) and (type(step) == type(0))):
        raise Exception("type(size) and type(step) must be int.")
    if step > size:
        raise Exception("step must not be larger than size.")
    if size > len(a):
        raise Exception("size must not be larger than array length.")
    a = np.asarray(a)
    h = ((len(a)-size)/step)+1
    indexer = np.arange(size)[None, :] + step*np.arange(int(h))[:, None]
    window = a[indexer]
    return window

def all_concat(array):
    iterate = int(len(array) / 2)
    conc = np.zeros(0)
    for i in range(iterate):
        start = array[i *2]
        stop = array[i * 2 + 1]
        conc = np.concatenate([conc, np.arange(start, stop+1)])
    return conc

pro_dict = {k:str(seq_dict[k][start_dict[k]:].translate(to_stop=True)) for k in transcripts}
wins = [10,9,8,7,6,5]
hits = {}
for w in wins:
    for name in transcripts:
        s = pro_dict[name]
        s_a = np.array(list(s))
        s_a = np.where((s_a=='D')|(s_a=='E'), 1, 0)
        s_a_w = getSlidingWindow(s_a, 10)
        starts = np.where((s_a_w[:,0]==1)&(s_a_w.sum(axis=1)>=w))[0]
        if len(starts>0):
            hits.setdefault(name, [starts[0], w, s[starts[0]:starts[0]+10]])
df = pd.DataFrame(hits, index=['motif', 'score', 'window']).T
df['gene_name'] = symbol.reindex(df.index)['symbol']
df['start'] = df['motif']*3+[start_dict[x] for x in df.index]
df['short'] = df['motif']<=40
df['long'] = df['motif']>=100
df.to_csv('./annotation/gencode_DEmotif_nascent_window.csv')

hits = []
w = 7
for name in transcripts:
    s = pro_dict[name]
    s_a = np.array(list(s))
    s_a = np.where((s_a=='D')|(s_a=='E'), 1, 0)
    s_a_w = getSlidingWindow(s_a, 10)
    starts = np.where((s_a_w[:,0]==1)&(s_a_w.sum(axis=1)>=w))[0]
    if len(starts>0):
        for st in starts:
            hits.append([name, st, str(s[st:st+10])])
df = pd.DataFrame(hits)
df.columns = ['Transcript ID', 'motif', 'window']
df['start'] = [start_dict[x] for x in df['Transcript ID']]
df['motif start'] = df['start']+df['motif']*3
df['motif end'] = df['motif start'] + 30
arrs = {}
for name, d in df.groupby('Transcript ID'):
    a = np.unique(all_concat(d[['motif start', 'motif end']].values.flatten()))
    split_indices = np.where(np.diff(a) != 1)[0] + 1
    split_arrays = np.split(a, split_indices)
    arr = np.array([[int(x[0]),int(x[-1])] for x in split_arrays])
    arrs[name] = arr
with open('./annotation/DE_position.pickle', 'wb') as f:
    pickle.dump(arrs, f)

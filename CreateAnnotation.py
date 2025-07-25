#!/usr/bin/env python
import os
import pandas as pd
import numpy as np
from Bio import SeqIO
import pickle

os.makedirs('./annotation', exist_ok=True)

gtf = pd.read_csv('./gencode.v31.annotation.gtf', sep='\t', header=None)

gtf_transcript = gencode[gencode[2] != 'gene']
gtf_transcript = gtf_transcript.drop([1,5,7], axis=1)
gtf_transcript.columns = ['chr','type','start','end','strand','annotation']
gtf_transcript['length'] = gtf_transcript['end']-gtf_transcript['start']
gtf_transcript['id'] = [x.split(';')[1].split('"')[1] for x in gtf_transcript['annotation']]
gtf_exon = gtf_transcriptf[gtf_transcript['type'] == 'exon']
gtf_codon = gtf_transcript[(gtf_transcript['type'] == 'start_codon')|(gtf_transcript['type'] == 'stop_codon')]
codons = gtf_codon.groupby('id')['type'].sum()
transcripts = list(codons[(codons.str.contains('start'))&(codons.str.contains('stop'))].index)
gtf_coding = gtf_transcripts[gtf_transcripts['id'].isin(transcripts)]
gtf_coding['symbol'] = [x.split(';')[3].split('"')[1] for x in gtf_coding['annotation']]
symbol = gtf_coding[['id','symbol']].drop_duplicates('id').set_index('id')

gtf_exon_coding = gtf_exon[gtf_exon['id'].isin(transcripts)]
exon_dict = {}
for name, group in gtf_exon_coding.groupby('id'):
    exon_dict[name] = group
    exon_dict[name] = exon_dict[name].reset_index()
with open('./annotation/coding_exon_dict.pickle', 'wb') as f:
    pickle.dump(exon_dict, f)
exon_rel_dict = {}
for item in exon_dict.keys():
    exon_rel_dict [item] = np.cumsum(np.array(exon_dict[item]['length']))

gtf_start = gtf_coding[gtf_coding['type'] == 'start_codon']
start_dict = {}
for name in range(len(transcripts)):
    starts = gtf_start[gtf_start['id'] == name]
    strand = starts['strand'].values[0]
    if strand == '+':
        start = starts['start'].values[0]
        start_exon = exon_dict[name][(exon_dict[name]['start'] <= start) & (exon_dict[name]['end'] >= start)]
        if start_exon.index == 0:
            start_dict[name] = start - start_exon['start'].values[0]
        else:
            start_dict[name] = start - start_exon['start'].values[0] + exon_rel_dict[name][start_exon.index[0]-1] + start_exon.index[0]
    else:
        start = starts['end'].values[0]
        start_exon = exon_dict[name][(exon_dict[name]['start'] <= start) & (exon_dict[name]['end'] >= start)]
        if start_exon.index == 0:
            start_dict[name] = start_exon.iloc[0,4] - start
        else:
            start_dict[name] = start_exon.iloc[0,4] - start + exon_rel_dict[start_exon.index[0]-1] + start_exon.index[0]
with open('./annotation/coding_start_dict.pickle', 'wb') as f:
    pickle.dump(start_exon_pos_dict, f)

fasta_in = './gencode.v31.transcripts.fa'
seq_dict = {}
for record in SeqIO.parse(fasta_in, 'fasta'):
    id_part = record.id
    name = id_part.split('|')[0]
    seq = record.seq
    if name in transcripts:
        seq_dict[id_part] = seq
with open('./annotation/sequence_dict.pickle', 'wb') as f:
    pickle.dump(seq_dict, f)

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

pro_dict = {k:str(seq_dict[k][start_dict[k]:].translate(to_stop=True)) for k in transcripts}
wins = [10,9,8,7]
hits = {}
for w in wins:
    for name in range(len(transcripts)):
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
df.to_csv('./gencode_DEmotif_nascent_window.csv')

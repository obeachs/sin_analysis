import sys
import Bio
import re
import itertools
from Bio import SeqIO
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.Seq import UnknownSeq
from Bio.Blast import NCBIXML
import numpy as np
import pandas as pd
import itertools
import argparse
import subprocess
import os
from os import listdir
from Bio import pairwise2

parser = argparse.ArgumentParser(description = 'fastax_len_filter.py Parameters\n')
parser.add_argument('--refcdna', required = False, default = '/Volumes/sesame/joerecovery/genomes/TAIR10_cdna_20101214_updated', type = str)
parser.add_argument('--queryfasta',required = False, default = '/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/sinapis_temp/sinapis_alba_var_s2_gc0560-79.maingenome.fasta', type = str)
parser.add_argument('--blastoutput', required = True, default = None, type= str)
parser.add_argument('--geneids', required = False, default = '/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/sinapis_temp/trichome_ids.txt', type= str)
parser.add_argument('--noextend', required=False, default=0,type=int)
args = parser.parse_args()


def fasta_condenser(fasta, tair=0):
    '''Preferably use this on the fasta file that has less info/annotation
    with it - usually the query fasta'''
    namelist = []
    seqlist = []
    if tair==0:
        with open (fasta,'r') as fa:
            for seq in SeqIO.parse(fa,'fasta'):
                namelist.append(seq.id)
                seqlist.append(str(seq.seq))
        df = pd.DataFrame(list(zip(namelist, seqlist)),columns =['ID', 'Seq'])
    else:
        with open (fasta,'r') as fa:
            for seq in SeqIO.parse(fa,'fasta'):
                namelist.append(seq.id[0:9])
                seqlist.append(str(seq.seq))
        df = pd.DataFrame(list(zip(namelist, seqlist)),columns =['TAIR_ID', 'TAIR_Seq'])
    return df





table = pd.read_table(args.blastoutput)
table['tair_id'] =table['tair_id'].replace(r'\..*', '', regex=True)
tair_fasta = args.refcdna
sinapis_genome = args.queryfasta
geneids = args.geneids
outfile = args.blastoutput
outfile = outfile.split('.')[0] + '_parsed.fa'
sinapis = fasta_condenser(sinapis_genome)
tair = fasta_condenser(tair_fasta,tair=1)



with open(geneids,'r') as ids:
    id_list = ids.read().splitlines() 

print(id_list)
tri_tair = tair[tair['TAIR_ID'].isin(id_list)]
print(table)
matches_table = table[table['tair_id'].isin(id_list)]
print(matches_table)
tri_sinapis = sinapis[sinapis['ID'].isin(matches_table['sinapis_id'])]
print(tri_sinapis)
merged = matches_table.merge(tri_tair, how='left',left_on='tair_id', right_on='TAIR_ID')
merged = merged.merge(tri_sinapis, how='left',left_on='sinapis_id', right_on='ID')
merged['set_seq'] = merged['Seq'].str.slice(merged['sinapis_start'],merged['sinapis_end'])
print(merged)
merged = merged[merged['Seq'].notna()]
sinapis_id = merged['sinapis_id'].to_list()
sinapis_seq = merged['Seq'].to_list()
sinapis_start = merged['sinapis_start'].to_list()
sinapis_end = merged['sinapis_end'].to_list()
tair_ids = merged['tair_id'].to_list()
tair_seq = merged['TAIR_Seq'].to_list()


print(sinapis_seq)


for i in range(len(sinapis_seq)):
    if args.noextend !=0:
        sinapis_seq[i] = sinapis_seq[i][sinapis_start[i]:sinapis_end[i]]
    else:
        try:
            sinapis_seq[i] = sinapis_seq[i][sinapis_start[i]-500:sinapis_end[i]-500]
        except:
            sinapis_seq[i] = sinapis_seq[i][sinapis_start[i]:sinapis_end[i]]


with open(outfile,'w+') as out: 
    for i in range(len(tair_ids)):
        print('>' + tair_ids[i] + '\n' + tair_seq[i] + '\n' + '>' + sinapis_id[i] + '\n' + sinapis_seq[i] + '\n')
        out.write('>' + tair_ids[i] + '\n' + tair_seq[i] + '\n' + '>' + sinapis_id[i] + '\n' + sinapis_seq[i] + '\n')








import sys
import Bio
import re
import itertools
from Bio import SeqIO
from Bio import SearchIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Seq import UnknownSeq
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW
import numpy as np
import pandas as pd
import itertools
import argparse
import subprocess
import sys
import os
from os import listdir
from Bio import pairwise2
import random


parser = argparse.ArgumentParser(description = 'fastax_len_filter.py Parameters\n')
parser.add_argument('--fasta', required = True, default = None, type = str)
parser.add_argument('--align', required = False, default=1,type= int)
args =parser.parse_args()


id_list = []
seq_list = []


fasta = args.fasta


with open(fasta,'r') as fa:
    for seq in SeqIO.parse(fa,'fasta'):
        id_list.append(str(seq.id))
        seq_list.append(str(seq.seq))

# function to get unique values
def unique(list1):
  
    # initialize a null list
    unique_list = []
      
    # traverse for all elements
    for x in list1:
        # check if exists in unique_list or not
        if x not in unique_list:
            unique_list.append(x)
    return unique_list


uni = unique(id_list)
'''This works well - finds the ID you're looking for then will scan forward and
take any other ids after that (the matches) until it reaches the next Arabidopsis
ID at which point it stops'''

for check in uni:
    if len(check) == 9 and check[0:1] =='A':
        for i in range(len(id_list)):
            if check == id_list[i]:
                print(check)
                print(id_list[i])
                newseq = []
                newid = []
                newseq.append(seq_list[i])
                newseq.append(seq_list[i+1])
                newid.append(id_list[i])
                newid.append(id_list[i+1])
                newseq_write = unique(newseq)
                newid_write = unique(newid)
                print(newid_write)
                result_fasta = fasta.split('.')[0] + '_' + check + '_matches.fa'
                with open(result_fasta,'w+') as result:
                    for i in range(len(newid_write)):
                        result.write('>' + newid_write[i] + '\n' + newseq_write[i] + '\n')



if args.align ==1:
    align = AlignIO.read(result_fasta, "fasta")
    with open(aligned_fasta,'w+') as align_final:
        SeqIO.write(align, align_final, 'nexus')

# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 22:43:59 2022

@author: Jing Li, Small steps make changes. dnt_seq@163.com
"""
import os
import pandas as pd
import numpy as np
pd.set_option('display.max_columns',1000)
pd.set_option('display.max_colwidth',1000)
pd.set_option('display.max_rows',50000)
pd.set_option('display.width', None)  # 设置字符显示宽度
np.set_printoptions(threshold=10000)
################

amino_table = ['I', 'D', 'M', 'H', 'E', 'W', 'R', 'L', 'Y', 'Q', 'G', 'A', 'S', 'P', 'C', 'T', 'V', 'F', 'N', 'K']

################
serial_num = 1

print (serial_num)
#####################

def cds_translator(seq):
    amino_list = ['F','F',\
                  'L','L','L','L','L','L',\
                  'S','S','S','S','S','S',\
                  'Y','Y',\
                  '*','*',\
                  'C','C',\
                  '*',\
                  'W',\
                  'P','P','P','P',\
                  'H','H',\
                  'Q','Q',\
                  'R','R','R','R','R','R',\
                  'I','I','I',\
                  'M',\
                  'T','T','T','T',\
                  'N','N',\
                  'K','K',\
                  'V','V','V','V',\
                  'A','A','A','A',\
                  'D','D',\
                  'E','E',\
                  'G','G','G','G']
    codon_table1 = ['ttt', 'ttc',\
                    'tta', 'ttg','ctt', 'ctc', 'cta', 'ctg',\
                    'tct', 'tcc', 'tca', 'tcg', 'agt', 'agc',\
                    'tat', 'tac',\
                    'taa', 'tag',\
                    'tgt', 'tgc',\
                    'tga',\
                    'tgg',\
                    'cct', 'ccc', 'cca', 'ccg',\
                    'cat', 'cac',\
                    'caa', 'cag',\
                    'cgt', 'cgc', 'cga', 'cgg','aga', 'agg',\
                    'att', 'atc', 'ata',\
                    'atg',\
                    'act', 'acc', 'aca', 'acg',\
                    'aat', 'aac',\
                    'aaa', 'aag',\
                    'gtt', 'gtc', 'gta', 'gtg',\
                    'gct', 'gcc', 'gca', 'gcg',\
                    'gat', 'gac',\
                    'gaa', 'gag',\
                    'ggt', 'ggc', 'gga', 'ggg']
    
    # print (len(amino_list))
    # print (len(codon_table1))
    codon_dict = dict(zip(codon_table1, amino_list))
    # print (codon_dict)
    seqlen = len(seq)
    # if seqlen % 3 != 0:
    #     return 'It is not a cds sequence!'
    # else:
        # seq = seq.lower()
    Protein = ''
    for codon_i in range(0,seqlen,3):
        codon = seq[codon_i:codon_i+3]
        if codon in codon_table1:
            Protein = Protein+codon_dict[codon]
        # else:
        #     Protein = Protein+'?'
    return Protein[:-1]

#####################

#####################

def count_codon (cds, amino_table = amino_table, Num4 = 20, freq = True):
    seq_len = len(cds)
    cds = cds.lower()

#################################################   amino acid counting
    count_amino = np.zeros(Num4)
    protein = cds_translator(cds.lower())
    seqlen = len(protein)
    for i in range(0, seqlen, 1):
        cut = protein[i:i+1]
        if cut in amino_table:
            amino = amino_table.index(cut)
            count_amino[amino] += 1

    if freq:
        return  (count_amino*20*3/seq_len)
    else:
        return count_amino

#####################


df_fasta1 = pd.read_csv (r'xxx.csv')
seq_num = df_fasta1.shape[0]

seq_list = df_fasta1.seq_Spike.tolist()
id_list = df_fasta1.seq_id.tolist()


lst_decom_all = []

for seq_i in range(10000):
    if seq_i % 1000 == 0:
        print (seq_i)
    seqID = id_list[seq_i]
    my_seq = seq_list[seq_i]
    seqlen = len(my_seq)
    # print (seqlen)
    decom_list = []
    for l in range(3, (seqlen), 3):
        if l <192:
            start_i1 = l-192
            start_i2 = 0
        else:
            start_i1 = seqlen
            start_i2 = l-192
        if l <=seqlen:
            end_i3 = 0
        else:
            end_i3 = l-seqlen

        seq_cut1 = my_seq[start_i1:]
        seq_cut2 = my_seq[start_i2:l]
        seq_cut3 = my_seq[:end_i3]
        seq_cut = seq_cut1 + seq_cut2 + seq_cut3
        decomposing = count_codon(seq_cut)
        lst_decom = decomposing.tolist()

        decom_list.append(lst_decom)
    lst_decom_all.append(decom_list)
array_decom_all = np.array(lst_decom_all)
print (array_decom_all.shape)
#print (array_decom_all[0,:,:])

df_lineage = pd.DataFrame()
df_lineage['seq_id'] = id_list
df_lineage['lineage'] = df_fasta1.lineage.tolist()

df_lineage.to_csv ('df_id_lineage.csv')

np.save ('array_sampled_lineage_oduageo_nt999_275121_0' + str(serial_num) + '.npy', array_decom_all)

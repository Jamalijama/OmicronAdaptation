# -*- coding: utf-8 -*-
"""
Created on Thu Aug 26 17:06:06 2021

@author: Jing Li, Small steps make changes. dnt_seq@163.com
"""
######################## import modules
import pandas as pd
import numpy as np
import torch
import torch.utils.data as Data
from torch.autograd import Variable
import torch.nn as nn
import torch.nn.functional as F
from sklearn.utils import shuffle
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
import os
from sklearn.preprocessing import label_binarize
from sklearn.utils.multiclass import unique_labels


path1 = '../decomposing/' 
           



###############################################################  CNN structure
###############################################################  CNN structure

class CNN (nn.Module):
    def __init__ (self):
        super (CNN, self).__init__()
        self.conv1 = nn.Sequential ( 
            nn.Conv2d ( in_channels = 1,
                        out_channels = 8,
                        kernel_size = (3, 3), # kernel only for 2d data
                        stride =(1,1),
                        padding = (1,1),
                        bias = True
                        ),            # 
            nn.ReLU (),
            nn.AvgPool2d (kernel_size = (2,2)) 
        )
        self.conv2 = nn.Sequential ( 
            nn.Conv2d ( in_channels = 8,
                        out_channels = 16,
                        kernel_size = (3, 3),# kernel only for 2d data
                        stride =(1,1),
                        padding = (1,1),
                        bias = True
                        ),            # 
            nn.ReLU (),
            nn.AvgPool2d (kernel_size = (2, 2)) # Max or Avg
        )
        self.conv3 = nn.Sequential (
            nn.Conv2d ( in_channels = 16,
                        out_channels = 32,
                        kernel_size = (3, 3),# kernel only for 2d data
                        stride =(1,1),
                        bias = True,
                        padding = (1,1)
                        ),
            nn.ReLU (),
            nn.AvgPool2d (kernel_size = (2, 2)) #MaxPool3d
        )
        
        self.conv4 = nn.Sequential (
            nn.Conv2d ( in_channels = 32,
                        out_channels = 64,
                        kernel_size = (3, 3),# kernel only for 2d data
                        stride =(1,1),
                        bias = True,
                        padding = (1,1)
                        ),
            nn.ReLU (),
            nn.AvgPool2d (kernel_size = (2, 2)) #MaxPool3d
        )

        self.fc1 = nn.Linear (5760, 960)  
        self.fc2 = nn.Linear (960, 160)  
        self.fc3 = nn.Linear (160, 3)

    def forward (self, x):
        x = self.conv1(x)
        x = self.conv2(x)
        x = self.conv3(x)
        x = self.conv4(x)

        x = x.view (x.size(0), -1) # flat x, similar to reshape of numpy
        fc_full = x
        x = F.sigmoid (self.fc1(x))  # to activate x   
        x = F.sigmoid (self.fc2(x))  # to activate x   
        prob_ = F.softmax (self.fc3(x))
        pred_ = (self.fc3(x))
        return pred_, prob_, fc_full
#        return x.shape

###############################################################  CNN structure
###############################################################  CNN structure

########################################### to load trained CNN model_
########################################### to load trained CNN model_

cnn = CNN()
if_use_gpu = 1
if if_use_gpu:
    cnn = cnn.cuda()
model_name = 'model_CNN2d_AAs_SARS2.txt'
b = 2
epoch_num_list = [30]
cnn = torch.load(model_name)

########################################################    prepare data
########################################################    prepare data

fasta_file = 'df_sampled3_lineage_oduageo_nt999_275121_concatenated_AAs_160000.csv'
df_fasta = pd.read_csv (path1 + fasta_file, index_col = 'seq_id')
print (df_fasta.shape)

df_fasta_test = df_fasta [df_fasta.lineage == 'Omicron']
print (df_fasta_test.shape)

index_list = df_fasta_test.seq_id.tolist()
 

ZZ = df_fasta_test.iloc[:,20:-1]
test_num = ZZ.shape[0]

y_test_list2 = ZZ.shape[0]*[2]
ZZ_array = np.array (ZZ)

ZZ_array = ZZ_array.reshape ([test_num, 1272, 20])
ZZ_array = ZZ_array.reshape ([test_num, 159, 160])
ZZ_test_tensor = torch.tensor(ZZ_array)

cnn = cnn.cpu()
test_pred_prob = cnn (Variable(ZZ_test_tensor)) #.cuda()
test_pred_matrix = test_pred_prob[0].cpu().detach().numpy()
test_pred_list = torch.max(valid_pred_prob[0],1)[1].data.cpu().detach().numpy()
test_pred_array = label_binarize (test_pred_list, classes = [0, 1, 2])

df_pred = pd.DataFrame (data = valid_pred_array, index = index_list)

df_pred.to_csv ('dataframe_predicted.csv')


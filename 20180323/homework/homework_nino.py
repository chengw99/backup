# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 18:33:29 2017

@author: DELL
"""

import pandas as pd
import os
import numpy as np
import sys
os.chdir('e:\data')
def re(var1,var2,sy=11):
    if len(var1)==len(var2):
        yy = len(var1)
    else:
        print('Error when calculating sliding correlation coefficients!')
        print('The input data has different len.')
        sys.exit()
    r = []
    for rnum in range(yy-sy+1):
        start_y = rnum
        end_y = rnum + sy -1
        p1,p21,p22 = 0.0,0.0,0.0
        var1_avg = np.average(var1[start_y:end_y+1])
        var2_avg = np.average(var2[start_y:end_y+1])
        
        for y in range(sy):
            p1 += (var1[start_y+y]-var1_avg)*(var2[start_y+y]-var2_avg)
            p21 += (var1[start_y+y]-var1_avg)**2
            p22 += (var2[start_y+y]-var2_avg)**2
        
        r.append(p1/np.sqrt(p21*p22))
    return r

prec_dir = './monthly_precip_1951_2007/'
prec_fn = ['Jan.txt','Jul.txt','Aug.txt']
sst_fn = 'nino_1951_2007.txt'
rm_fn = 'readme.txt'

output_fn = 'r_01.txt'

sta_lon = []
sta_lat = []
sta_num = []

f = open(prec_dir + rm_fn,'r') #f.name可以得到所读txt文件的文件名
tmp_data = f.readlines()
for i in tmp_data[21:]:
    tmp_list = i.strip().split()
    sta_lon.append(float(tmp_list[-1]))
    sta_lat.append(float(tmp_list[-2]))
    sta_num.append(int(tmp_list[-3]))
del(tmp_list)

prec_tmp = np.zeros((160,57,3),dtype=np.float)
prec_data = np.zeros((160,57),'float')

for mon,fname in enumerate(prec_fn):
    f = open(prec_dir + fname,'r')
    tmp_data = f.readlines()
    for year in range(57):
        for i in range(8):
            prec_tmp[i*20:(i+1)*20,year,mon] = tmp_data[year*8+i].split()
prec_data = np.average(prec_tmp,axis = -1)

del(prec_tmp)

#计算海温
sst_data = np.zeros(57,'float')

f = open(sst_fn,'r')
tmp_data = f.readlines()

for year in range(57):
    for month in range(3):
        sst_data[year] += float(tmp_data[12 +year*12+month].split()[4])
sst_data =sst_data/3.
            
r_prec_sst = np.zeros((160,47),dtype=np.float)

for sta in range(160):
    r_prec_sst = re(sst_data[:],prec_data[sta,:],11)

std_sta = np.zeros(160,'float')

for sta in range(160):
    std_sta = np.std(r_prec_sst[sta,:], ddof=1)
    
f=open(output_fn,'w')
for i in range(160):
    line = str(sta_num[i])+''+str(sta_lon[i])+''+str(sta_lat[i])+''+str(std_sta[i])+'\n'
    f.write(line)
    

  
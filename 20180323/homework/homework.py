# -*- coding: utf-8 -*-
"""
Created on Sun Oct  8 17:41:22 2017

@author: DELL
"""

import pandas as pd
import os
import numpy as np

os.chdir("e://data")
#a=pd.read_table('Jun.txt')

#中国夏季降水资料
a=np.loadtxt('Jun.txt')
b=np.loadtxt('Jul.txt')
c=np.loadtxt('Aug.txt')
summer_ave=((a+b+c)/3.).reshape(10720,1) #在此修改数组的维度
d=[]
for i in range(57):
    d.append(summer_ave[i*160:i*160+160]) #得到的summer_ave是二维数组，所以这样索引是不对的,需要修改上面的数组维度
summer = np.array(d).reshape(160,57)


#前期冬季海温
data=pd.read_table('nino_1951_2007.txt',delim_whitespace=True)
f=data['NINO3']
winter=[]
for j in range(57):
    winter.append(sum(f[12*j+11:12*j+14])/3.)

#滑动平均公式
'''def mov_ave(var1,var2,n=11):
    r=[]
    for k in range'''

x_mean=[]
winterf=[]
for k in range(6):
    y=np.array(d[k*11:(k+1)*11])-np.array(sum(d[11*k:11*(k+1)])/11.)#-sum(d[11*k:11*k+11])/11.
    h=np.array(winter[k*11:k*11+11])-np.array(sum(winter[k*11:k*11+11])/11.)
    x_mean.append(y)
    winterf.append(h)
    
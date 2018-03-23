# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 17:19:11 2018

@author: DELL
"""

import os
#import glob
import shutil

#------------可遍历指定目录，并生成对应文件夹--------------------
csv_filenames = glob.glob('e:\\py_output\\*.csv')
os.chdir('e:\\py_output')
for filename in csv_filenames:
    isExists = os.path.exists(filename[-6]) #os.path.exists()检验指定的对象是否存在，可用来判断文件或文件夹的存在与否
    if not isExists:                        #os.path.isfile()只能判断对象是否为文件
        os.mkdir(filename[-6]) #创建文件夹
    else:
        pass
#----------------------------------------------------------

#------此程序用于批量移动文件，适用文件类型：ibutton字母+数字-分离成不同日期的文件---------
#此处分析的文件格式是 A1.csv    
#如果是其他格式，需要更改下面的切片数字
#简单来说就是筛选特定格式的文件存入指定的文件夹
folderaddress = 'e:\\py_output'      #目标路径
file_all = os.listdir(folderaddress) #路径下所有文件名
folderlist = []                      #存放文件夹的命名
os.chdir('e:\\py_output')            #修改工作目录

#创建文件夹
for name in file_all:
    if name[-3:] == 'csv':            #此处可设定移动文件的类型 
        folderlist.append(name[:2])  #此处获取所需对应的文件夹名
        isExists = os.path.exists(name[0:2]) #判断文件夹是否存在
        if not isExists:                     #加一个判断语句，就不会因为文件已存在而报错
            os.mkdir(folderaddress+'\\'+name[0:2])
        else:
            pass

#移动文件    
for i in file_all:                    #遍历所有目标文件
    for j in folderlist:              #遍历文件夹 
        #判断文件是否存在以及文件名与文件夹名是否一致，进行定向筛选
        if os.path.isfile(folderaddress+'\\'+i) == True and i[0:2] == j:
            #移动文件到目标文件夹
            shutil.copyfile(folderaddress+'\\'+i,folderaddress+'\\'+j+'\\'+i)
            #删除原来位置的文件
            os.remove(folderaddress+'\\'+i)

#-------------------移动某个文件夹到另一个文件夹------------------------------------           
filenames = os.listdir('e:\\py_output')
for file in filenames:
    shutil.move('e:\\py_output\\'+file,r'e:\\研一') #  r''可以处理含有中文的目录 

#------------------复制文件到指定文件夹-------------------------------------------
file_data = os.listdir(r'C:\Users\DELL\Desktop\处理数据\源数据\ibutton数据\5月份')
for file in file_data:
    if file[-3:] == 'csv':
        shutil.copyfile(r'C:\Users\DELL\Desktop\处理数据\源数据\ibutton数据\5月份\\'+file,'e:\\py_data\\ibutton\\'+file) #

#-----------------在每个文件夹里生成一个命名为每小时平均的文件夹-----------------------
filenames = os.listdir(r'C:\Users\DELL\Desktop\处理数据\分析数据\ibutton')
for file in filenames:
    csv_filenames = os.listdir(r'C:\Users\DELL\Desktop\处理数据\分析数据\ibutton\\'+file)
    isExist = os.path.exists(r'C:\Users\DELL\Desktop\处理数据\分析数据\ibutton\\'+file+'\\'+r'每小时平均')
    if not isExist:
        #os.mkdir(r'C:\Users\DELL\Desktop\处理数据\分析数据\ibutton\\'+file+'\\'+r'每小时平均')#

#----------------批量重命名文件 -------------------------------------------------    
filename = os.listdir('e:\py_data\pic') #遍历目标文件夹
for file in filename:
    #重命名文件
    shutil.move('e:\py_data\pic\\'+file,'e:\py_data\pic\\'+file[-5]+'.csv') #（原文件名，新文件名）

#可以保留原文件  
filename = os.listdir('e:\py_data\pic') #遍历目标文件夹
for file in filename:
    #复制并重命名新文件
    shutil.copy('e:\py_data\pic\\'+file,'e:\py_data\\'+file[-5]+'.csv')  #（原目录\原文件名，新目录\新文件名）
    
    
    



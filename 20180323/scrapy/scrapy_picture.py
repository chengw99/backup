# -*- coding: utf-8 -*-
"""
Created on Wed Dec 20 17:19:17 2017

@author: chen
"""

import re 
import requests
import os
#import time

word = input('Input key word :')
url = 'http://image.baidu.com/search/flip?tn=baiduimage&ie=utf-8&word='+word+'&ct=201326592&v=flip'
html = requests.get(url).text
pic_url = re.findall('"objURL":"(.*?)",',html,re.S)
#查看是否存在存放目录
try:
    if not os.path.exists('E:\\pics\\'):
        os.mkdir('E:\\pics\\')
except:
    pass
#爬取
t = 0
i = 0
while t<20:
    for each in pic_url:
        print(each)
        try:
            pic = requests.get(each,timeout=10)
            #time.sleep(1)
        except requests.exceptions.ConnectionError:
            print('【错误】当前图片无法下载')
            continue
        string = 'E:\\pics\\'+str(i)+'.jpg'
        with open(string,'wb') as f:
            f.write(pic.content)
            f.close()

        i+=1
        t+=4

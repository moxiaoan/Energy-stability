#!/usr/bin/python
# coding=utf-8

import matplotlib.pyplot as plt
import numpy as np
from numpy import *
import pylab as pl

filename = 'phi_ac_3_energy_four_long_1.txt'
#filename = 'phi_ac_3_energy.txt'

X,Y = [],[]

with open(filename, 'r') as f:#1
    lines = f.readlines()#2
    for line in lines:#3
        value = [float(s) for s in line.split()]#4
        X.append(value[0])#5
        Y.append(value[1]) 
  
#print(X)
#print(Y)
#plt.title('Energy Stability curve of β to α',loc='center')
plt.title('Energy Stability curve of β to ω',loc='center')
plt.xlabel('step')
plt.ylabel('energy')
#plt.axis([200,6800])
#plt.xlim((200,6800))
plt.ylim((110000,670000))
ax = pl.gca()  # 获取当前图像的坐标轴信息
ax.xaxis.get_major_formatter().set_powerlimits((0,1)) # 将坐标轴的base number设置为一位。1是指科学计数法时的位数
ax.yaxis.get_major_formatter().set_powerlimits((0,1)) # 将坐标轴的base number设置为一位。1是指科学计数法时的位数

plt.vlines(200, 110000, 610362, colors='r', linestyle="dashed")
plt.text(200, 610362, "   start evolution")

plt.annotate("nucleation", (50,110000), xycoords='data',
             xytext=(200, 150000), 
             arrowprops=dict(arrowstyle='->'))

j = len(X)

for i in range (j):
    if (i <= 6700 and (i+200) % 1000 == 0 or i == 0):
    #if (i <= 6500 and i % 1000 == 0 or i == 0):
        k = X[i]
        l = Y[i]
        plt.plot(k, l,marker = "h", color = 'r', markersize = 5)
    else:
        #plt.plot(X, Y, linestyle = 'dashed', color = 'g')
        plt.plot(X, Y, linestyle = '-', color = 'g', linewidth=1)
#plt.plot(X,Y)
plt.savefig("test1.png", dpi=600,format="png")
#plt.show()


# -*- coding: utf-8 -*-
"""
Created on Sun Nov 26 09:46:56 2017

@author: DELL
"""
import math
import matplotlib.pyplot as plt
import numpy as np

#---------------------------初始系数--------------------------#
L = 1.
N = 50
h = L/N
c = 1.
sigma = 0.1
k_wave = math.pi/sigma

#-------------------------------------------------------------#
#---------------------FTCS------------------------------------#
def v(t):
    X=[]
    A=[]
    tau = 0.002
    coeff = -c*tau/(2.*h)
    for i in range(1,N+1):
        x = (i-1)*h - L/2.
        a = math.cos(k_wave*x)*math.exp(-x**2/(2.*sigma**2))
        X.append(x)
        A.append(a)
    
    for i in range(int(t/tau)):
        A_new = []
        for j in range(49):
            a_new = A[j] + coeff*(A[j+1] - A[j-1])
            A_new.append(a_new)
        A_new_bottom = A[49] + coeff*(A[0] - A[48])
        A_new.append(A_new_bottom)
        A = A_new
    return np.array(A)

X = np.linspace(-0.5,0.48,50)    #另一种表达形式

plt.figure(1)

plt.plot(X,v(0.0),'k-',label = r'$φ_{0}(x)$')
plt.plot(X,v(0.1),'g--',label = r'$φ(x,t=0.1)$')
plt.plot(X,v(0.6),'r--',label = r'$φ(x,t=0.6)$')
plt.plot(X,v(1.0),'y--',label = r'$φ(x,t=1.0)$')
plt.xlabel('x',fontsize = 18)
plt.ylabel('φ(x,t)',fontsize = 18)
plt.title('FTCS  scheme τ=0.002',fontsize = 20)
plt.xlim(-0.6,0.6)
plt.xticks([-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6])
plt.ylim(-2,2.4)
plt.yticks([-2,-1,0,1,2])
plt.legend()
plt.show()
#plt.savefig('FTCs scheme τ=0.002.png')

#----------------------------------------------------------#
#-------------------------Lax------------------------------#
def w(t,tau):
    X=[]
    B=[]
    coeff = -c*tau/(2.*h)  
    for i in range(1,N+1):
        x = (i-1)*h - L/2.
        b = math.cos(k_wave*x)*math.exp(-x**2/(2.*sigma**2))
        X.append(x)
        B.append(b)
    
    for i in range(int(t/tau)):
        B_new = []
        for j in range(49):
            b_new = 1./2*(B[j+1]+B[j-1]) + coeff*(B[j+1] - B[j-1])
            B_new.append(b_new)
        B_new_bottom = 1./2*(B[0]+B[48]) + coeff*(B[0] - B[48])
        B_new.append(B_new_bottom)
        B = B_new
    return B

X = np.linspace(-0.5,0.48,50)  #.....#

plt.figure(figsize=(16,16))
#-----------------------------------------------------------
plt.subplot(221)
plt.plot(X,w(0.0,0.001),'k-',label = r'$φ_{0}(x)$')
plt.plot(X,w(0.1,0.001),'g--',label = r'$φ(x,t=0.1)$')
plt.plot(X,w(0.6,0.001),'r--',label = r'$φ(x,t=0.6)$')
plt.plot(X,w(1.0,0.001),'y--',label = r'$φ(x,t=1.0)$')
plt.title('Lax  scheme τ=0.001',fontsize = 20)
plt.xlabel('x',fontsize = 18)
plt.ylabel('φ(x,t)',fontsize = 18)
plt.xlim(-0.6,0.6)
plt.xticks([-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6])
plt.ylim(-0.8,1.2)
plt.yticks([-0.8,0,0.8,1.2])
plt.legend()

#---------------------------------------------------------
plt.subplot(222)
plt.plot(X,w(0.0,0.005),'k-',label = r'$φ_{0}(x)$')
plt.plot(X,w(0.1,0.005),'g--',label = r'$φ(x,t=0.1)$')
plt.plot(X,w(0.6,0.005),'r--',label = r'$φ(x,t=0.6)$')
plt.plot(X,w(1.0,0.005),'y--',label = r'$φ(x,t=1.0)$')
plt.title('Lax  scheme τ=0.005',fontsize = 20)
plt.xlabel('x',fontsize = 18)
plt.ylabel('φ(x,t)',fontsize = 18)
plt.xlim(-0.6,0.6)
plt.xticks([-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6])
plt.ylim(-0.8,1.2)
plt.yticks([-0.8,0,0.8,1.2])
plt.legend()
#---------------------------------------------------------

plt.subplot(223)
plt.plot(X,w(0.0,0.01),'k-',label = r'$φ_{0}(x)$')
plt.plot(X,w(0.1,0.01),'g--',label = r'$φ(x,t=0.1)$')
plt.plot(X,w(0.6,0.01),'r--',label = r'$φ(x,t=0.6)$')
plt.plot(X,w(1.0,0.01),'y--',label = r'$φ(x,t=1.0)$')
plt.title('Lax  scheme τ=0.01',fontsize = 20)
plt.xlabel('x',fontsize = 18)
plt.ylabel('φ(x,t)',fontsize = 18)
plt.xlim(-0.6,0.6)
plt.xticks([-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6])
plt.ylim(-0.8,1.2)
plt.yticks([-0.8,0,0.8,1.2])
plt.legend()
#--------------------------------------------------------
plt.subplot(224)
plt.plot(X,w(0.0,0.02),'k-',label = r'$φ_{0}(x)$')
plt.plot(X,w(0.1,0.02),'g--',label = r'$φ(x,t=0.1)$')
plt.plot(X,w(0.6,0.02),'r--',label = r'$φ(x,t=0.6)$')
plt.plot(X,w(1.0,0.02),'y--',label = r'$φ(x,t=1.0)$')
plt.title('Lax  scheme τ=0.02',fontsize = 20)
plt.xlabel('x',fontsize = 18)
plt.ylabel('φ(x,t)',fontsize = 18)
plt.xlim(-0.6,0.6)
plt.xticks([-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6])
plt.ylim(-0.8,1.2)
plt.yticks([-0.8,0,0.8,1.2])
plt.legend()
#plt.savefig('Lax sheme.png')

#plt.show()
#plt.savefig('Lax scheme τ=0.002.png')'''
#--------------------Lax-Wendroff-----------------------------------------#
def u(t,tau):
    X=[]
    C=[]
    coeff = -1.*tau/(2.*h)  
    for i in range(1,N+1):
        x = (i-1)*h - L/2.
        c = math.cos(k_wave*x)*math.exp(-x**2/(2.*sigma**2))
        X.append(x)
        C.append(c)
    
    for i in range(int(t/tau)):
        C_new = []
        for j in range(49):
            c_new = C[j] + coeff*(C[j+1] - C[j-1]) + 2*(coeff**2)*(C[j+1]-2*C[j]+C[j-1])
            C_new.append(c_new)
        C_new_bottom = C[49] + coeff*(C[0] - C[48]) +2*(coeff**2)*(C[0]-2*C[49]+C[48])
        C_new.append(C_new_bottom)
        C = C_new
    return C

X = np.linspace(-0.5,0.48,50)  #.....#

plt.figure(figsize=(16,16))
#-------------------------------------------------
plt.subplot(221)
plt.plot(X,u(0.0,0.001),'k-',label = r'$φ_{0}(x)$')
plt.plot(X,u(0.1,0.001),'g--',label = r'$φ(x,t=0.1)$')
plt.plot(X,u(0.6,0.001),'r--',label = r'$φ(x,t=0.6)$')
plt.plot(X,u(1.0,0.001),'y--',label = r'$φ(x,t=1.0)$')
plt.title('Lax-Wendroff  scheme τ=0.001',fontsize = 20)
plt.xlabel('x',fontsize = 18)
plt.ylabel('φ(x,t)',fontsize = 18)
plt.xlim(-0.6,0.6)
plt.xticks([-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6])
plt.ylim(-0.8,1.2)
plt.yticks([-0.8,0,0.8,1.2])
plt.legend()
#-------------------------------------------------
plt.subplot(222)
plt.plot(X,u(0.0,0.005),'k-',label = r'$φ_{0}(x)$')
plt.plot(X,u(0.1,0.005),'g--',label = r'$φ(x,t=0.1)$')
plt.plot(X,u(0.6,0.005),'r--',label = r'$φ(x,t=0.6)$')
plt.plot(X,u(1.0,0.005),'y--',label = r'$φ(x,t=1.0)$')
plt.title('Lax-Wendroff  scheme τ=0.005',fontsize = 20)
plt.xlabel('x',fontsize = 18)
plt.ylabel('φ(x,t)',fontsize = 18)
plt.xlim(-0.6,0.6)
plt.xticks([-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6])
plt.ylim(-0.8,1.2)
plt.yticks([-0.8,0,0.8,1.2])
plt.legend()
#-------------------------------------------------------
plt.subplot(223)
plt.plot(X,u(0.0,0.01),'k-',label = r'$φ_{0}(x)$')
plt.plot(X,u(0.1,0.01),'g--',label = r'$φ(x,t=0.1)$')
plt.plot(X,u(0.6,0.01),'r--',label = r'$φ(x,t=0.6)$')
plt.plot(X,u(1.0,0.01),'y--',label = r'$φ(x,t=1.0)$')
plt.title('Lax-Wendroff  scheme τ=0.01',fontsize = 20)
plt.xlabel('x',fontsize = 18)
plt.ylabel('φ(x,t)',fontsize = 18)
plt.xlim(-0.6,0.6)
plt.xticks([-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6])
plt.ylim(-0.8,1.2)
plt.yticks([-0.8,0,0.8,1.2])
plt.legend()
#-------------------------------------------------------
plt.subplot(224)
plt.plot(X,u(0.0,0.02),'k-',label = r'$φ_{0}(x)$')
plt.plot(X,u(0.1,0.02),'g--',label = r'$φ(x,t=0.1)$')
plt.plot(X,u(0.6,0.02),'r--',label = r'$φ(x,t=0.6)$')
plt.plot(X,u(1.0,0.02),'y--',label = r'$φ(x,t=1.0)$')
plt.title('Lax-Wendroff  scheme τ=0.02',fontsize = 20)
plt.xlabel('x',fontsize = 18)
plt.ylabel('φ(x,t)',fontsize = 18)
plt.xlim(-0.6,0.6)
plt.xticks([-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6])
plt.ylim(-0.8,1.2)
plt.yticks([-0.8,0,0.8,1.2])
plt.legend()
#plt.savefig('Lax-Wendroff scheme.png')
# -*- coding: utf-8 -*-
"""
Created on Tue May 14 18:03:18 2019

@author: Cristina
"""
import math
import os
import numpy as np
from matplotlib import rcParams, rc
import matplotlib.pyplot as pl
from scipy.optimize import curve_fit
rcParams['font.family'] = 'serif'
rcParams['font.sans-serif'] = ['serif']


################################# Required ####################################
W =[226.1455545, 124.5538737, 84.99446781]
H_total = [2.0, 4.0, 6.0]
nlayers = list(range(1,4))

base_element = 1
L = 1.0#*math.cos(math.pi/6) # [mm] ,0.0,5.0
H = 1.0#math.cos(math.pi/6) #[mm]
shape=''# number of layers to test
leng = L*10*nlayers[-1]




#################___________CHOOSE-> BASE ELEMENT___________###################
be_dictionary ={0:'Line', 1:'Triangle', 2:'DoubleTriangle', 3:'SemiHoneycomb', 
                4:'Rectangle', 5:'SixTetCube', 6:'Mixedcube', 11:'SymTriangle'}


#################----> Creating results folder
# Set current direction
parent_folder = os.getcwd()
def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)
# Creates a folder in the current directory called data
# Jump inside the results folder
h = 0.1#[mm]        
folder = './' + be_dictionary[base_element] + shape +'_T_h' + str(h)+'mm_it'+str(nlayers[-1])+'/'
createFolder(folder)
os.chdir(folder)

E = 2000.0; # [MN/mm2]
nu = 0.25
G = E/(2*(1-nu)); # [MN/mm2]
rho = 1040.0 # [Kg/mm3] weight
b = 0.1 #[mm]
A = b*h # [mm2]
Iyy = (b*h*h*h)/12 # [mm4] moment of inertia about z axis
ky  = 5.0/6.0 # rigid in shear
Aboundary = b*h/2
Iyyboundary = (b*(h/2)**3)/12

force = (53.0)# [MN]

###########----> File writting###############
file = open('Results_COMSOL.txt','w') 
D = []
file.write('\n')
file.write('[height        [mm]; average_displ [mm]; D/E             ]\n')
for i in range(len(H_total)):
    file.write('%s ;' %H_total[i] )
    file.write('%s ;' %W[i])
    D.append(4.0*(float(leng/H_total[-1])**3)*force/(W[i]*b))
    file.write('%s ]\n' %D[i])
file.write('Curve fitting with least squares: \n')
# Fitting curve preparation
C = [G*ky*b*H_total[-1], (leng/H_total[-1])]
file.write('Timoshenko gradient constant = %6.3f' %C[0])
file.write('l/a = %6.3f' %C[1])
file.write('Resulting parameters\n')
def funcT(h, g, EE): #Timoshenko
    h = np.asanyarray(h)
    return (4*C[1]**2*C[0]*(12*g**2 + h**2))*EE/(h*(4*C[1]**2*C[0]*h+12*EE*b*g**2+EE*b*h**2))#(12*C[0]*(1+12*(g*g)/(h*h))/(4*C[0]*C[1]*h**2+(EE*b*h**2)*(1+12*(g*g)/(h*h))*EE
def func(h, g, EE): #Euler-Bernoulli
    h = np.asanyarray(h)
    return (1+12*(g*g)/(h*h))*EE
# Fitting resulting curve
poptT, pcovT = curve_fit(funcT, np.asanyarray(H_total[:]), np.asanyarray(D[:]))
popt, pcov = curve_fit(func, np.asanyarray(H_total[:]), np.asanyarray(D[:]))
DT = [d/poptT[1] for d in D]
DB = [d/popt[1] for d in D]
file.close() 
# Plotting
labelEcB = 'Results Bernoulli Ec = ' + '%.3f' %popt[1] + 'MPa'
labelEcT ='Results Timoshenko Ec = ' + '%.3f' %poptT[1] + 'MPa'
pl.plot(H_total, DB, 'ro', label=labelEcB)
pl.plot(H_total, DT, 'bx', label=labelEcT)
hint = np.arange(H_total[0],H_total[-1],1.0) ##hint = np.arange(H_total[0],1000,1)
#pl.plot(hint, funcT(hint, *poptT)/poptT[1], 'g:',label='fit: g=%.3f' % poptT[0] + ' & Ec=%.3f MPa (TG)' % popt[1], lw=1)
#pl.plot(hint, func(hint, *popt)/popt[1], 'g:',label='fit: g=%.3f' % popt[0] + ' & Ec=%.3f MPa (BG)' % popt[1], lw=1)
pl.xlabel('total height [mm]')
pl.ylabel('D/D_0')
pl.legend()
pl.plot(H_total, [1 for d in D],color='grey',ls='--', label='', lw=1)
pl.ylabel('$D/D_{0}$')
pl.autoscale()
pl.title(be_dictionary[base_element])
pl.savefig('COMSOL_NBD_' +be_dictionary[base_element] +'.pdf')
pl.show()

#x = np.arange(1, 10)
#y = x.reshape(-1, 1)
#h = x * y
#delta = 0.025
#X = np.arange(0.001,0.6,0.001)
#Y = np.arange(0,1, 0.0001)
#x,y = np.meshgrid(X,Y)
#Z = (1+12*(popt[0]/(np.multiply(x,x))))/(1+2e9*(np.divide(x,(np.multiply(y,y)*C[0])))*(1+12*(popt[0])/(np.multiply(x,x))))
#fig, ax = pl.subplots()
#CS = ax.contour(X, Y, Z,20)
#pl.autoscale()
#ax.clabel(CS, inline=1, fontsize=10)
#ax.set_title('Simplest default with labels')
#cs.cmap.set_over('red')
#cs.cmap.set_under('blue')
#cs.changed()

# Get back to the working directory
os.chdir(parent_folder)
print('------------------------------END---------------------------------')

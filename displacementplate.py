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
folder = './PLATE_T_h' +'mm_it'+str(nlayers[-1])+'/'
createFolder(folder)
os.chdir(folder)
D = np.array([[226.1455545,124.5538737,84.99446781,64.32184361,50.94557941,43.18000672]])
D = D/np.min(D)
H_total = np.reshape(np.arange(2,14,2.),(1,6))
# Plotting
pl.plot(H_total, D, 'ro', label='')
hint = np.arange(H_total[0,0],H_total[0,-1],1.0) ##hint = np.arange(H_total[0],1000,1)
#pl.plot(hint, funcT(hint, *poptT)/poptT[1], 'g:',label='fit: g=%.3f' % poptT[0] + ' & Ec=%.3f MPa (TG)' % popt[1], lw=1)
#pl.plot(hint, func(hint, *popt)/popt[1], 'g:',label='fit: g=%.3f' % popt[0] + ' & Ec=%.3f MPa (BG)' % popt[1], lw=1)
pl.xlabel('Total height [mm]')
pl.ylabel('Normalized displacement')
pl.legend()
#pl.plot(H_total, [1 for d in D],colo1r='grey',ls='--', label='', lw=1)
pl.autoscale()
pl.title('Plate. Octet-truss')
pl.savefig('COMSOL_NBD_' +'Plate. Octet-truss'+'.pdf')
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

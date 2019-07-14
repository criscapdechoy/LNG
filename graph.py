# -*- coding: utf-8 -*-
"""
Created on Sun Jul 14 20:36:49 2019

@author: jccap
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams, rc

import pandas as pd
rcParams['font.family'] = 'serif'
rcParams['font.sans-serif'] = ['serif']


data = pd.read_csv("data.csv",sep=';')
n = np.array(data['n'])
e = np.array(data['e'])	
f = np.array(data['f'])
t = np.array(data['t'])


#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
### for Palatino and other serif fonts use:
##rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)
plt.loglog(n,t,'o-r',linewidth=2, markersize=4, label = 'nodes')
plt.loglog(e,t,'o-b',linewidth=2, markersize=4, label = 'edges')
plt.loglog(f,t,'o-g',linewidth=2, markersize=4, label = 'faces')
plt.loglog(np.array([n[0],n[-1]]),np.array([t[0],t[-1]]),linewidth=1,color='grey', label = 'O(n)')
plt.legend()
plt.axis('scaled')
plt.xlabel('number of elements', fontsize = 12)
plt.ylabel('running time [s]', fontsize = 12)
plt.show()
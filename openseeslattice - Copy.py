#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 18:47:00 2019

@author: cristinacapdevilachoy
"""
import openseespy.opensees as op
import sys
import base
import timeit
import math
import os
import numpy as np
from matplotlib import rcParams
import matplotlib.pyplot as pl
from scipy.optimize import curve_fit




#################___________CHOOSE-> BASE ELEMENT___________###################


#################----> Creating results folder
# Set current direction



####################___________MODEL_OPENSEES___________#######################
## ------------------------------
## Start of model generation
## -----------------------------
# remove existing model
op.wipe()

# set modelbuilder

## generate a material
Ecc = 12.5 # [MN/mm2] Core Young's modulus
Ec = 12.5; # [MN/mm2] Core Young's modulus
E = 2000.0; # [MN/mm2]
nu = 0.25
G = E/(2*(1-nu)); # [MN/mm2]
rho = 1040 # [Kg/mm3] weight

## generate a cross section inside
b = 1.0 #[mm]
A = b*h # [mm2]
Iyy = (b*h*h*h)/12 # [mm4] moment of inertia about z axis
#Izz = 0.002 # [mm4] moment of inertia about y axis
#Ipp = 0.003 # [mm4] torsional moment of inertia
ky  = 5.0/6.0 # rigid in shear
#kz  = 0.0 # rigid in shear

# generata the cross section outside (upper&bottom edges)
Aboundary = b*h/2
Iyyboundary = (b*(h/2)**3)/12

# longitudinal dimensions of each beam bottom
L = 5.0 # [mm]
H = 5.0*math.cos(math.pi/6) #[mm]
# load case:
force = 4.0e-6 # [MN]
op.wipe()
op.model('basic', '-ndm', 2,'-ndf', 3)
op.geomTransf('Linear', 1)
op.node(0, *[0.0,0.0])
op.node(1, *[5.0,0.0])
op.element('ElasticTimoshenkoBeam', 0, 0,1, E, G, Aboundary, Iyyboundary, Aboundary*ky, 1) #2D
op.fix(0, 1,1,1) #0free 1fix
force_p = float(force)
op.timeSeries('Constant',1)
op.pattern('Plain',1,1)
op.load(1, 0.0,20*force_p,0.0)
op.constraints('Transformation')
op.numberer('RCM')
op.system('FullGeneral')
op.integrator('LoadControl',1.0)
op.algorithm('Linear')
op.analysis('Static')
op.analyze(1)
u0, w0 ,rot0 = 0.0, 0.0, 0.0
u, w, rot = op.nodeDisp(1)
print("results")
print("-T-")
print('Average displacements at x = L: \n')
print("  ux= %.16f[mm]" % u)
print("  uz= %.16f[mm]" % w) 
print("  rot= %.16f[mm]" % rot) 
##############################################################
op.wipeAnalysis()
op.wipe()
op.model('basic', '-ndm', 2,'-ndf', 3)
op.geomTransf('Linear', 1)
op.node(0, *[0.0,0.0])
op.node(1, *[5.0,0.0])
op.element('elasticBeamColumn', 0,0,1, Aboundary, E, Iyyboundary, 1)
op.fix(0, 1,1,1) #0free 1fix
force_p = float(force)
op.timeSeries('Constant',1)
op.pattern('Plain',1,1)
op.load(1, 0.0,force_p,0.0)
op.constraints('Transformation')
op.numberer('RCM')
 #   op.test('NormUnbalance', 1e-20, 100)
op.system('FullGeneral')
op.integrator('LoadControl',1.0)
op.algorithm('Linear')
op.analysis('Static')
op.analyze(1)
 
# calculate member forces and displacements
u0, w0 ,rot0 = 0.0, 0.0, 0.0
#    filedisp.write('Iteration %s\n' % i) 

#        filedisp.write('%s; ' % j) 
#        filedisp.write('%s; %s;\n'%(nodelist[j].getDisplacement(0).x(), nodelist[j].getDisplacement(0).z()))
u, w, rot = op.nodeDisp(1)

# Print in console and file
print ("---")
print ("Iteration:"+ str(i))
print ("---")
# Dimensions of the sandwhich beam
print("---")
print("results")
print("-EB-")


print('Average displacements at x = L: \n')
print("  ux= %.16f[mm]" % u)
print("  uz= %.16f[mm]" % w) 
print("  rot= %.16f[mm]" % rot) 

# Opensees file model creation


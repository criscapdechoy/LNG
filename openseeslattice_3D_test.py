# -*- coding: utf-8 -*-

"""

Created on Sun Apr 28 18:41:04 2019



@author: Cristina

"""



import openseespy.opensees as op

import LNG_builder, LNG_engine, LNG_main

import timeit
import math
import os

import numpy as np

from matplotlib import rcParams, rc

import matplotlib.pyplot as pl

from scipy.optimize import curve_fit

#################___________CHOOSE-> BASE ELEMENT___________###################

be_dictionary ={0:'Bar', 1:'Triangle', 2:'DoubleTriangle', 3:'SemiHoneycomb',

                4:'Rectangle', 5:'Cube', 6:'X', 7:'hexa_tetra', 8:'octet', 9:'FCC',10:'CFCC',11:'BFCC'}

# RTV, TypeElem:

    # 2D:

        # bar: 0

        # half-triangle: 1 or 11 filled 2.5,0 , 5.0*math.cos(math.pi/6)

        # double triangle: 2

        # semihoney: 3 5*math.cos(math.pi/6),0.0,5.0

        # rectangle: 4

    # 3D:

        # cube: 5

        # X shape: 6

        # hexa_tetra: 7

        # octet-truss: 8

    # Shape of the global volume, shape:

        # cuboid: ''

        # cylinder: 'cylinder'

        # sphere: 'sphere'

base_element = 8

shape=''

perc = 0.0

ilay = 2

n_G = [20,20,6]

#################----> Creating results folder

# Set current direction

parent_folder = os.getcwd()

# Creates a folder in the current directory called data

# Jump inside the results folder

h = 0.1#[mm]

folder = './SIMPLE_STRUCTURE_ANALYSIS'

LNG_engine.createFolder(folder)

os.chdir(folder)

folder = './TEST' + shape +'t'+str(base_element)+'s'+ ''.join(map(str, n_G[:])) + 'r' +str(perc)

LNG_engine.createFolder(folder)

os.chdir(folder)

####################___________MODEL_OPENSEES___________#######################

## ------------------------------

## Start of model generation

## -----------------------------

# remove existing model

op.wipe()

# set modelbuilder

## generate a material

Ecc = 262.5 # [MN/mm2] Core Young's modulus

Ec = 121.0 # [MN/mm2] Core Young's modulus

E = 2000.0# [MN/mm2]

nu = 0.25

G = E/(2*(1-nu)); # [MN/mm2]

rho = 1040.0 # [Kg/mm3] weight



## generate a cross section inside

b = h #[mm]

A = b*h # [mm2]

Iyy = (b*h**3)/12 # [mm4] moment of inertia about z axis

Izz = (h*b**3)/12 # [mm4] moment of inertia about y axis

J = b*h*(h*h+b*b)/12 # [mm4] torsional moment of inertia

ky  = 5.0/6.0 # rigid in shear

kz  = 5.0/6.0 # rigid in shear

# generata the cross section outside (upper&bottom edges)

h = h/2

Abound = b*h # [mm2]

Iyybound = (b*h**3)/12 # [mm4] moment of inertia about z axis

Izzbound = (h*b**3)/12 # [mm4] moment of inertia about y axis

Jbound = b*h*(h*h+b*b)/12 # [mm4] torsional moment of inertia

kybound  = 5.0/6.0 # rigid in shear

kzbound  = 5.0/6.0 # rigid in shear

h=2*h

# longitudinal dimensions of each beam bottom

L = 1.0 # [mm]

H = 1.0#0*math.cos(math.pi/6) #[mm]



# load case:

force = 1.0 # [MN]

## ------------------------------

## Generate summary file .txt

## ------------------------------

file = open('Results.txt','w')

filedisp = open('Displacements.txt','w')

filedisp.write('u; v; w; Fx; Fy; Fz\n')

file.write('MODEL COMMON PARAMETERS\n')

file.write('-  Material data:\n')

file.write('      E = %E [MN/mm2]\n' % E)

file.write('      nu = %E \n' % nu)

file.write('      G = %E [MN/mm2]\n' % G)

file.write('      rho %E [MN/mm3]\n' % rho)

file.write('-  Longitudinal dimension data:\n')

file.write('      Length of the base element, L = %E [mm]\n' % L)

file.write('      Height of the base element, H = %E [mm]\n' % H)

file.write('-  Cross section data:\n')

file.write('      b = %E [mm]\n' % b)

file.write('      h = %E [mm]\n' % h)

file.write('      A = %E [mm2]\n' % A)

file.write('      Iyy = %E [mm4]\n' % Iyy)

file.write('      Izz = %E [mm4]\n' % Izz)

file.write('      ky = %E \n' % ky)

file.write('      kz = %E \n' % kz)

file.write('-  Load:\n')

file.write('      Point load, f = %E [MN]\n\r' % force)

file.write('\n\r')

rcParams['font.family'] = 'serif'

rcParams['font.sans-serif'] = ['serif']



# Starting the iteration



W = []

D=[]

H_total = []



# Start timer

start_time = timeit.default_timer()

op.wipe()

op.model('basic', '-ndm', 3,'-ndf', 6)



#    op.geomTransf('Linear', 1, 1.0,1.0,0.0) # X beam

op.geomTransf('Linear', 1, 0.0,0.0,1.0) # Z beam

op.geomTransf('Linear', 2, 0.0,1.0,0.0) # Y beam



# Call domain generation function base_element, shape, [18*i,1,i],

LS = LNG_builder.main(base_element, shape, n_G, meshSize=[], sizeXYZ=[L,L,L], Export=False, View=False)

# nodes of the domain

N  = len(LS.nodes)

# generate beam elements from edges inside

ii = -1

for n0 in LS.edges:

    try:

        op.node(n0[0], *LS.nodes[n0[0],:])

    except:

        pass

    for n1 in n0[1:]:

        try:

            op.node(n1, *LS.nodes[n1,:])

        except:

            pass

        transfTag = 2 - (sum([LS.nodes[n0[0],o]-LS.nodes[n1,o] for o in [0,1]]) != 0) # = 1 if not column or = 2 for column

        ii= ii+1

        op.element('elasticBeamColumn', ii, n0[0],n1, A, E, G, J, Iyy,Izz, int(transfTag))#int(transfTag)) #3D

#        if not(((n0[0] in LS.boundary[4])and(n1 in LS.boundary[4]))or((n0[0] in LS.boundary[5])and(n1 in LS.boundary[5]))or(

#                (n0[0] in LS.boundary[0])and(n1 in LS.boundary[0]))or((n0[0] in LS.boundary[1])and(n1 in LS.boundary[1]))or(

#                (n0[0] in LS.boundary[2])and(n1 in LS.boundary[2]))or((n0[0] in LS.boundary[3])and(n1 in LS.boundary[3]))):

#            op.element('elasticBeamColumn', ii, n0[0],n1, A, E, G, J, Iyy,Izz, int(transfTag))#int(transfTag)) #3D

##        else:

#            op.element('elasticBeamColumn', ii, n0[0],n1, Abound, E, G, Jbound, Iyybound,Izzbound, int(transfTag))#int(transfTag)) #3D



# generate beam elements from edges conforming the boundary



# generate boundary conditions



for node in LS.boundary[2]:

    op.fix(node, 1,1,1,1,1,1) #0free 1fix

# Load

force_p = force/len(LS.boundary[3]) # divide load to be equally distributed at each right node

op.timeSeries('Constant',1)

op.pattern('Plain',1,1)

i=0

while i < len(LS.boundary[3])-1:

    op.load(LS.boundary[3][i], 0.0,0.0,force_p,0.0,0.0,0.0)

#        op.rigidDiaphragm('beam', LS.boundary[1][-1],  LS.boundary[1][i])

    i = i+1

op.load(LS.boundary[3][-1], 0.0,0.0,force_p,0.0,0.0,0.0)

#    op.rigidDiaphragm(1, LS.boundary[1][-1],  *LS.boundary[1][:-1])



# set up analysis procedure

op.constraints('Transformation')

op.numberer('RCM')

 #   op.test('NormUnbalance', 1e-20, 100)

op.system('BandGen')

op.integrator('LoadControl',1.0)

op.algorithm('Linear')

op.analysis('Static')

op.analyze(1)

op.reactions()

# calculate member forces and displacements

u0, v0, w0 = 0.0, 0.0, 0.0

Fill = []

Disps = []

#filedisp.write('Iteration %s\n' % i)

for node in LS.boundary[3]:

    U = op.nodeDisp(node)

    u0 = U[0] + u0

    v0 = U[1] + v0

    w0 = U[2] + w0

for node in range(N):

    U = op.nodeDisp(node)

    F = op.nodeReaction(node)

    Fill.append(F)

    Disps.append(U)

    filedisp.write(str(U[0]) + ';'+str(U[1]) + ';'+str(U[2]) + ';'+str(F[0]) + ';'+str(F[1]) + ';'+ str(F[2]) + ';\n')

del(U)

u = u0/len(LS.boundary[3])

v = u0/len(LS.boundary[3])

w = w0/len(LS.boundary[3])

# Plotting





# Print in console and file

print ("---")

print ("Iteration:"+ str(ilay))

print ("---")

# Dimensions of the sandwhich beam

print("---")

print("results")

print("---")

#file.write('Iteration: %s\n' % i)

#file.write('   Total lenght of the hole beam: %.2f [mm]\n' % float(n_G[0]*L))

#file.write('   Total height of the hole beam: %.2f [mm]\n' % float(i*H))

#file.write('RESULTS: \n')



print('Average displacements at x = L: \n')

print("  ux= %.16f[mm]" % u)

print("  uz= %.16f[mm]" % v)

print("  uz= %.16f[mm]" % w)

file.write('   average displacement at x = L: \n')

file.write('     ux= %.16f[mm]\n' % u)

file.write('     uy= %.16f[mm]\n' % v)

file.write('     uz= %.16f[mm]\n' % w)

# Opensees file model creation





#    H_total.append(i*H)

# Finish timer

elapsed = timeit.default_timer() - start_time

print("elapsed time: %E seconds" % elapsed )

file.write('elapsed time: %E seconds \n' % elapsed)

file.write('############################################################## \n\r')

filedisp.write('############################################################## \n\r')

op.wipeAnalysis()

del(elapsed, w0, w, v0, v, u0, u, ii)



filedisp.close()

# Write Curve Data to Results.txt file

file.close()

#

## Get back to the working directory

os.chdir(parent_folder)

print('------------------------------END---------------------------------')





#[op.nodeCoord(i) for i in range(len(coordinate))]

#nodeDOFs(nodeTag)

    #reactions('-dynamic', '-rayleight')

    #nodeResponse(nodeTag, dof, responseID) Disp = 1

#• Vel = 2#• Accel = 3#• IncrDisp = 4#• IncrDeltaDisp = 5#• Reaction = 6#• Unbalance = 7#• RayleighForces = 8



#1. fix command

#2. fixX command

#3. fixY command

#4. fixZ command



#Create constraints for multiple dofs of multiple nodes.

#1. equalDOF command

#2. equalDOF_Mixed command

#3. rigidDiaphragm command

#4. rigidLink command





#rigidLink(type, rNodeTag, cNodeTag)

#

#LS.nodes('line', tag, numnodes, *ndtags, id, ndf, LS.nodessize, eleType=”, *eleArgs=[])

#LS.nodes('tri', tag, numlines, *ltags, id, ndf, LS.nodessize, eleType=”, *eleArgs=[])

#element('ElasticTimoshenkoBeam', eleTag, *eleNodes, E, G, A, Iz, Avy, transfTag[, '-mass', massDens ][,'-cMass ]) #2D

#

#element('Tri31', eleTag, *eleNodes, thick, type, matTag[, pressure, rho, b1, b2 ])

#nDMaterial('ElasticIsotropic', matTag, E, v, rho=0.0)


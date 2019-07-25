#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 18:47:00 2019

@author: cristinacapdevilachoy
"""
import openseespy.opensees as op
import base
import base_func_flowchart as base_func
import timeit
import math
import os
import numpy as np
from matplotlib import rcParams, rc
import matplotlib.pyplot as pl
from scipy.optimize import curve_fit
#################___________CHOOSE-> BASE ELEMENT___________###################
be_dictionary ={0:'Bar', 1:'Triangle', 2:'DoubleTriangle', 3:'SemiHoneycomb', 
                4:'Rectangle', 5:'Cube', 6:'X', 7:'hexa_tetra', 8:'octet', 11:'SymTriangle'}
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
base_element = 4
shape=''

nlayers = 7
n_G = [18*nlayers,1,nlayers]
#################----> Creating results folder
# Set current direction
parent_folder = os.getcwd()

# Creates a folder in the current directory called data

# Jump inside the results folder
h = 0.5#[mm]        
folder = './SIMPLE_STRUCTURE_ANALYSIS'
base_func.createFolder(folder)
folder = './TEST' + shape +'t'+str(base_element)+'s'+ ''.join(map(str, n_G[:]))
base_func.createFolder(folder)
os.chdir(folder)
####################___________MODEL_OPENSEES___________#######################
## ------------------------------
## Start of model generation
## -----------------------------
# remove existing model
op.wipe()

# set modelbuilder
## generate a material

E = 2000.0; # [MN/mm2]
nu = 0.25
G = E/(2*(1-nu)); # [MN/mm2]
rho = 1040.0 # [Kg/mm3] weight
## generate a cross section inside
b = 1.0 #[mm]
A = b*h # [mm2]
Iyy = (b*h*h*h)/12 # [mm4] moment of inertia about z axis
ky  = 5.0/6.0 # rigid in shear
# generata the cross section outside (upper&bottom edges)
Aboundary = b*h/2
Iyyboundary = (b*(h/2)**3)/12
# longitudinal dimensions of each beam bottom
L = 2.5#*math.cos(math.pi/6) # [mm] ,0.0,5.0
H = 5.0*math.cos(math.pi/6) #[mm]

# load case:
force = 4.0e-3 # [MN]

###########----> File writting###############
file = open('Results.txt','w') 
filedisp = open('Displacements.txt','w') 
filedisp.write('DISPLACEMENTS (u,w) at each node [mm]\n\r') 
file.write( str(base_element) + '_shape_' + shape+ '_height_' + str(h) + '[mm] \n') 
file.write('Timoshenko beam with (Openseespy)\n') 
file.write('\nAmount of simulations: '+ str(nlayers)) 
file.write('\n\rMODEL COMMON PARAMETERS \n') 
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
file.write('      ky = %E \n' % ky)
file.write('-  Cross section data bread beams:\n')
file.write('      A_bread = %E [mm2]\n' % Aboundary)
file.write('      Iyy_bread = %E [mm4]\n' % Iyyboundary) 
file.write('-  Load:\n') 
file.write('      Point load, f = %E [MN]\n\r' % force)
file.write('\n\r') 
file.write('###################--STARTING--ITERATION--###################\n\r')  
rcParams['font.family'] = 'serif'
rcParams['font.sans-serif'] = ['serif']
#rc('font',**{'family':'serif'})
rc('text', usetex=False)
# Starting the iteration
W =[]
D = []
H_total = []

#for i in1 nlayers:
i = nlayers
# Start timer
start_time = timeit.default_timer()
op.wipe() #cler model
op.model('basic', '-ndm', 2,'-ndf', 3) #initialize model
op.geomTransf('Linear', 1) #define geometric transformation

# Call domain generation function
mesh, e_i, faces_i, B = base.main(base_element, shape, n_G, meshSize=[h/2, h*(2**0.5)/2], sizeXYZ=(L,0,H), t=0, ExtractGeom=False, View=False)
  #  op.nDMaterial('ElasticIsotropic', 1, E, nu, rho)

# nodes of the domain
ii = -1
ilayer = 0
ii0 = 0
for node in mesh:
    ii = ii +1
    op.node(ii, *mesh[ii,[0,2]])
totN = ii
#        try:
#            if mesh[ii+1,2] > mesh[ii,2]:
#                ii0 = ii0 + 1 #int((ii0 in  B[0]))
#                print(ii0)
#                for iii in range(ii0, ii):
#                    
#                    op.equalDOF(ii, iii,*[1,2])
#                ii0 = ii+1
#                print(ii0)
#        except:
#            pass
# generate beam elements from edges inside
ii= -1
Bele=set()
for n0 in e_i:
    for n1 in n0[1:]:  
        ii = ii+1
        if ((n0[0] in B[4])and(n1 in B[4])) or ((n0[0] in B[5])and(n1 in B[5])):
             op.element('elasticBeamColumn', ii, n0[0],n1, Aboundary, E, Iyyboundary, 1)
             Bele.add(ii)
        else:
            op.element('elasticBeamColumn', ii, n0[0],n1, A, E, Iyy, 1)
#            op.element('ElasticTimoshenkoBeam', ii, n0[0],n1, E, G, A, Iyy, A*ky, 1) #2D

# generate beam elements from edges belonging to the boundary
#for BB in B[4:]:
#    for Bb in BB:
#        for Bcon in e_i[Bb][1:]:
#            if Bcon in BB:
#                ii = ii +1
##                op.element('ElasticTimoshenkoBeam', ii, e_i[Bb][0], Bcon, E, G, Aboundary, Iyyboundary, Aboundary*ky, 1) #2D
#                op.element('elasticBeamColumn', ii, e_i[Bb][0], Bcon, Aboundary, E, Iyyboundary, 1)    
del(e_i)

# generate boundary conditions
for node in B[0]: 
    op.fix(node, 1,1,1) #0free 1fix
    modB = divmod(len(B[1]),2)
#    if modB[1]!=0:
#        op.fix(B[1][modB[0]], 1,0,0) #0free 1fix
# Load
force_p = float(float(force)/float(len(B[1])))
# divide load to be equally distributed at each right node
op.timeSeries('Constant',1)
op.pattern('Plain',1,1)
k = 0
for k in range(len(B[1])-1):
    op.load(B[1][k], 0.0 ,force_p, 0.0)
#        op.equalDOF(B[1][-1],  B[1][k], *[2,3])
    op.rigidLink('beam', B[1][-1],  B[1][k])
op.load(B[1][-1], 0.0,force_p,0.0)
##        if k != modB[0]: 
#            op.equalDOF(B[1][modB[0]],  B[1][k], *[2,3])
##        k = k+1
##    op.rigidLink('bar', B[1][-1],  B[1][0])
#    B[0].sort()
#    nlayer = B[0][1]-B[0][0]
#    for k in range(1,len(B[4])):
#        for kk in range(1,i+1):
#            op.equalDOF(B[4][k],  B[4][k]+kk*nlayer, *[2,3])
##    

#    op.load(B[1][-1], 0.0,force_p,0.0)
#     set up analysis procedure 
op.constraints('Transformation')
op.numberer('RCM')
 #   op.test('NormUnbalance', 1e-20, 100)
op.system('FullGeneral')
op.integrator('LoadControl',1.0)
op.algorithm('Linear')
op.analysis('Static')
op.analyze(1)
op.reactions()
# calculate member forces and displacements
Fill= []
Disps = []
u0, w0 ,rot0 = 0.0, 0.0, 0.0
#    filedisp.write('Iteration %s\n' % i) 
stress=[]
for node in range(len(mesh)):
    U = op.nodeDisp(node)
    F = op.nodeReaction(node)
    Fill.append(F)
    Disps.append(U)
    if node in Bele:
        stress.append(Disps[-1][2]*h/(4*Iyyboundary))
    else:
        stress.append(Disps[-1][2]*h/(2*Iyy))
    filedisp.write(str(U[0]) + ';'+str(U[1]) + ';'+str(F[0]) + ';'+str(F[1]) + ';\n')
    if node in B[1]:
        u0 = float(U[0] + u0)
        w0 = float(U[1] + w0)

u = u0/float(len(B[1]))
w = w0/float(len(B[1]))
Fill = np.asarray(Fill)
Disps = np.asarray(Disps)

filedisp.write('##########################################') 
               
               
# Plotting
pl.scatter(mesh[:,0], Disps[:,1]*10e2+mesh[:,2], c=Disps[:,1])
pl.show()
#pl.scatter(mesh[:,0], Disps[:,1]*10e2+mesh[:,2], c=Disps[:,1])
#pl.show()
# Print in console and file
print ("---")
print ("Number of layers:"+ str(i))
print ("---")
# Dimensions of the sandwhich beam
print("---")
print("results")
print("---")
file.write('Iteration: %s\n' % i) 
file.write('   Total lenght of the hole beam: %.8f [mm]\n' % float(n_G[0]*L)) 
file.write('   Total height of the hole beam: %.8f [mm]\n' % float(i*H))
file.write('RESULTS: \n') 
file.write(str(n_G))
file.write('\n RESULTS: \n') 

print('Average displacements at x = L: \n')
print("  ux= %.16f [mm]" % u)
print("  uz= %.16f [mm]" % w) 
file.write('   average displacement at x = L: \n') 
file.write('     ux= %.16f [mm]\n' % u) 
file.write('     uz= %.16f [mm]\n' % w) 
# Opensees file model creation
# D coef:
leng = float(L*n_G[0])
height = float(H*i)
Iyy_total = b*(height**3)/12
W.append(w)
#    D.append([force/w,3*Ec*Iyy_total*(leng**-3)])
#    DD.append([4.0*18*(float(n_G[0])/n_G[1])*force/(w*Ec*b)])
H_total.append(i*H)
# Finish timer
elapsed = timeit.default_timer() - start_time
print("elapsed time: %E seconds" % elapsed )
file.write('elapsed time: %E seconds \n' % elapsed) 
file.write('############################################################## \n\r') 
filedisp.write('############################################################## \n\r') 
op.wipeAnalysis()
del(  elapsed, i, ii)
    
filedisp.close() 
file.close()
# Write Curve Data to Results.txt file
#file.write('\n')
#file.write('[height        [mm]; average_displ [mm]; D/E             ]\n')

#for i in range(len(H_total)):
#    file.write('%s ;' %H_total[i] )
#    file.write('%s ;' %W[i])
#    D.append(4.0*(float(leng/H_total[-1])**3)*force/(W[i]*b))
#    file.write('%s ]\n' %D[i])
#file.write('Curve fitting with least squares: \n')
## Fitting curve preparation
#C = [G*ky*b*H_total[-1], (leng/H_total[-1])]
#file.write('Timoshenko gradient constant = %6.3f' %C[0])
#file.write('l/a = %6.3f' %C[1])
#file.write('Resulting parameters\n')
#
#def funcT(h, g, EE): #Timoshenko
#    h = np.asanyarray(h)
#    return (4*C[1]**2*C[0]*(12*g**2 + h**2))*EE/(h*(4*C[1]**2*C[0]*h+12*EE*b*g**2+EE*b*h**2))#(12*C[0]*(1+12*(g*g)/(h*h))/(4*C[0]*C[1]*h**2+(EE*b*h**2)*(1+12*(g*g)/(h*h))*EE
#def func(h, g, EE): #Euler-Bernoulli
#    h = np.asanyarray(h)
#    return (1+12*(g*g)/(h*h))*EE
## Fitting resulting curve
#poptT, pcovT = curve_fit(funcT, np.asanyarray(H_total[:]), np.asanyarray(D[:]))
#    
#popt, pcov = curve_fit(func, np.asanyarray(H_total[:]), np.asanyarray(D[:]))
#DT = [d/poptT[1] for d in D]
#DB = [d/popt[1] for d in D]
## Writting
##file.write('gB = %6.6f' %poptT[0])
##file.write('gT = %6.6f' %popt[0])
##file.write('EcB = %6.6f' %poptT[1])
##file.write('EcT = %6.6f' %popt[1])
#file.close() 
#
## Plotting
#labelEcB = 'Results Bernoulli Ec = ' + '%.3f' %popt[1] + 'MPa'
#labelEcT ='Results Timoshenko Ec = ' + '%.3f' %poptT[1] + 'MPa'
#pl.plot(H_total, DB, 'ro', label=labelEcB)
#pl.plot(H_total, DT, 'bx', label=labelEcT)
##l.plot(H_total, divvv, 'c-', label='drawing: g=%5.3f' % gg)
##pl.plot(H_total, funcT(H_total, *popt)/popt[1], 'g--',label='fit: g=%5.3f' % poptT[0] + ' & Ec=%5.4f MPa' % poptT[1])
#hint = np.arange(H_total[0],H_total[-1],1.0) ##hint = np.arange(H_total[0],1000,1)
#pl.plot(hint, funcT(hint, *poptT)/poptT[1], 'g:',label='fit: g=%.3f' % poptT[0] + ' & Ec=%.3f MPa (TG)' % popt[1], lw=1)
#pl.plot(hint, func(hint, *popt)/popt[1], 'g:',label='fit: g=%.3f' % popt[0] + ' & Ec=%.3f MPa (BG)' % popt[1], lw=1)
#
#pl.xlabel('total height [mm]')
#pl.ylabel('D/D_0')
#pl.legend()
#pl.plot(H_total, [1 for d in D],color='grey',ls='--', label='', lw=1)
#pl.ylabel('$D/D_{0}$')
#pl.autoscale()
#pl.title(be_dictionary[base_element])
#pl.savefig('NBD_' +be_dictionary[base_element] +'.pdf')
#pl.show()

# Get back to the working directory
os.chdir(parent_folder)
print('------------------------------END---------------------------------')


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 18:47:00 2019

@author: cristinacapdevilachoy
"""
import openseespy.opensees as op
import base
import timeit
import math
import os
import numpy as np
from matplotlib import rcParams
import matplotlib.pyplot as pl
from scipy.optimize import curve_fit

#################___________CHOOSE-> BASE ELEMENT___________###################
be_dictionary ={0:'Line', 1:'Triangle', 2:'DoubleTriangle', 3:'SemiHoneycomb', 
                4:'Rectangle', 5:'SixTetCube', 6:'Mixedcube', 11:'SymTriangle'}
base_element = 1
shape=''
nlayers = list(range(1,11)) # number of layers to test
h = 0.25 #[mm]
n_G = []
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
folder = './' + be_dictionary[base_element] + shape +'_bar_h' + str(h)+'mm_it'+str(nlayers[-1])+'/'
createFolder(folder)
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
nu = 0.01
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
#L = 2.5 # [mm]
#H = 5.0*math.cos(math.pi/6) #[mm]
L = 5.0*math.cos(math.pi/6) # [mm]
H = 5.0 #[mm]
# load case:
force = 4.0*10.0e-6 # [MN]

###########----> File writting
file = open('Results.txt','w') 
#filedisp = open('Displacements.txt','w') 
#filedisp.write('DISPLACEMENTS (u,w) at each node [mm]\n\r') 
file.write( str(base_element) + '_shape_' + shape+ '_height_' + str(h) + '[mm] \n') 
file.write('Euler-Bernoulli beam with (Openseespy)\n') 
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
file.write('-  Cross section data bread beams:\n')
file.write('      A_bread = %E [mm2]\n' % Aboundary)
file.write('      Iyy_bread = %E [mm4]\n' % Iyyboundary) 
file.write('-  Load:\n') 
file.write('      Point load, f = %E [MN]\n\r' % force)
file.write('\n\r') 
file.write('###################--STARTING--ITERATION--###################\n\r')  
rcParams['font.family'] = 'serif'
rcParams['font.sans-serif'] = ['serif']

# Starting the iteration

W = []
D=[]
H_total = []

for i in nlayers:
    # Start timer
    start_time = timeit.default_timer()
    op.wipe()
    op.model('basic', '-ndm', 2,'-ndf', 3)
    op.geomTransf('Linear', 1)
    
    n_G = [i*2*18,1,i]
    # Call domain generation function
    mesh, e_i, faces_i, B = base.main(base_element, shape, n_G, meshSize=[h/2, h*(2**0.5)/2], sizeXYZ=(L,0,H), t=0, ExtractGeom=False, View=False)
  #  op.nDMaterial('ElasticIsotropic', 1, E, nu, rho)

    # nodes of the domain
    ii = -1
    for node in mesh:
        ii = ii +1
        op.node(ii, *mesh[ii,[0,2]])
        if (ii not in B[0]):
            op.fix(ii, 0,1,0) #0free 1fix
    del(mesh)
    # generate beam elements from edges inside
    ii= -1
    for n0 in e_i:
        for n1 in n0[1:]:
            if not(((n0[0] in B[4])and(n1 in B[4]))or((n0[0] in B[5])and(n1 in B[5]))or(
                    (n0[0] in B[0])and(n1 in B[0]))or((n0[0] in B[1])and(n1 in B[1]))):
                ii = ii+1
                op.element('elasticBeamColumn', ii, n0[0],n1, A, E, Iyy, 1) #2D
    # generate beam elements from edges conforming the bread

    for BB in B:
        for Bb in BB:
            for Bcon in e_i[Bb][1:]:
                if Bcon in BB:
                    ii = ii +1
                    op.element('elasticBeamColumn', ii, e_i[Bb][0], Bcon, Aboundary, E, Iyyboundary, 1) #2D
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
        op.load(B[1][k], -force_p, 0.0,0.0)
        op.rigidLink('beam', B[1][-1],  B[1][k])
    op.load(B[1][-1], -force_p,0.0,0.0)
        
#        if k != modB[0]: 
#            op.equalDOF(B[1][modB[0]],  B[1][k], *[2,3])
#        k = k+1
#    op.rigidLink('bar', B[1][-1],  B[1][0])
#    B[0].sort()
#    nlayer = B[0][1]-B[0][0]
#    for k in range(1,len(B[4])):
#        for kk in range(1,i+1):
#            op.equalDOF(B[4][k],  B[4][k]+kk*nlayer, *[2,3])
#    

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
    
    
    # calculate member forces and displacements
    u0, w0 ,rot0 = 0.0, 0.0, 0.0
#    filedisp.write('Iteration %s\n' % i) 
    for node in B[1]:

#        filedisp.write('%s; ' % j) 
#        filedisp.write('%s; %s;\n'%(nodelist[j].getDisplacement(0).x(), nodelist[j].getDisplacement(0).z()))
        u, w, rot = op.nodeDisp(node)
        u0 = float(u + u0)
        w0 = float(w + w0)
        rot0 = float(rot + rot0)
    del(u,w,rot)
    u = u0/float(len(B[1]))
    w = w0/float(len(B[1]))
    rot = rot0/len(B[1])
    # Print in console and file
    print ("---")
    print ("Iteration:"+ str(i))
    print ("---")
    # Dimensions of the sandwhich beam
    print("---")
    print("results")
    print("---")
    file.write('Iteration: %s\n' % i) 
    file.write('   Total lenght of the hole beam: %.2f [mm]\n' % float(n_G[0]*L)) 
    file.write('   Total height of the hole beam: %.2f [mm]\n' % float(i*H))
    file.write(str(n_G))
    file.write('RESULTS: \n') 

    print('Average displacements at x = L: \n')
    print("  ux= %.16f[mm]" % u)
    print("  uz= %.16f[mm]" % w) 
    file.write('   average displacement at x = L: \n') 
    file.write('     ux= %.16f[mm]\n' % u) 
    file.write('     uz= %.16f[mm]\n' % w) 
    # Opensees file model creation

    
    # D coef:
    leng = float(L*n_G[0])
    height = float(H*i)
    Iyy_total = b*(height**3)/12
    W.append(u)
    H_total.append(i*H)
    # Finish timer
    elapsed = timeit.default_timer() - start_time
    print("elapsed time: %E seconds" % elapsed )
    file.write('elapsed time: %E seconds \n' % elapsed) 
    file.write('############################################################## \n\r') 
#    filedisp.write('############################################################## \n\r') 
    op.wipeAnalysis()
    del( B, elapsed)
    
#filedisp.close() 
# Write Curve Data to Results.txt file
file.write('\n')
file.write('[height        [mm]; average_displ [mm]; D/E             ]\n')


## BEEEEEEEEEEEEEEEEEEAAM
#for i in range(len(H_total)):
#    file.write('%s ;' %H_total[i] )
#    file.write('%s ;' %W[i])
#    D.append(4.0*((leng/H_total[-1])**3)*force/(W[i]*b))
#    file.write('%s ]\n' %D[i])
#file.write('Curve fitting with least squares: \n')
## Fitting curve preparation
#C = [G*ky*b*H_total[-1], (leng/H_total[-1])]
#file.write('Timoshenko gradient constant = %6.3f' %C[0])
#file.write('l/a = %6.3f' %C[1])
#file.write('Resulting parameters\n')
#def funcT(h, g, EE): #Timoshenko
#    h = np.asanyarray(h)
#    return (4*C[1]**2*C[0]*(12*g**2 + h**2))*EE/(h*(4*C[1]**2*C[0]*h+12*EE*b*g**2+EE*b*h**2))#(12*C[0]*(1+12*(g*g)/(h*h))/(4*C[0]*C[1]*h**2+(EE*b*h**2)*(1+12*(g*g)/(h*h))*EE
def func(h, g, EE): #Euler-Bernoulli
    h = np.asanyarray(h)
    return (1+12*(g*g)/(h*h))*EE
## Fitting resulting curve
#poptT, pcovT = curve_fit(funcT, np.asanyarray(H_total[5:]), np.asanyarray(D[5:]))
#popt, pcov = curve_fit(func, np.asanyarray(H_total[5:]), np.asanyarray(D[5:]))
#poptT =1,1
#popt =1,1
#DT = [d/poptT[1] for d in D]
#DB = [d/popt[1] for d in D]
#labelEcB = 'Results Bernoulli Ec = ' + '%.3f' %popt[1] + 'MPa'
#labelEcT ='Results Timoshenko Ec = ' + '%.3f' %poptT[1] + 'MPa'
#pl.plot(H_total, DB, 'ro', label=labelEcB)
#pl.plot(H_total, DT, 'bx', label=labelEcT)
#hint = np.arange(H_total[0],H_total[-1],1.0) ##hint = np.arange(H_total[0],1000,1)
#pl.plot(hint, funcT(hint, *poptT)/poptT[1], 'g:',label='fit: g=%.3f' % poptT[0] + ' & Ec=%.3f MPa (TG)' % popt[1], lw=1)
#pl.plot(hint, func(hint, *popt)/popt[1], 'g:',label='fit: g=%.3f' % popt[0] + ' & Ec=%.3f MPa (BG)' % popt[1], lw=1)

# BAAAAAAAAAAAAAAAAAAAAAAAAR!
C = [G*ky*b*H_total[-1], (leng/H_total[-1])]
D_b=[]
for i in range(len(H_total)):
    file.write('%s ;' %H_total[i] )
    file.write('%s ;' %W[i])
    D_b.append(force*C[1]/(W[i]*H_total[i]))
    file.write('%s ]\n' %D_b[i])
file.write('Curve fitting with least squares: \n')
# Fitting curve preparation
C = [G*ky*b*H_total[-1], (leng/H_total[-1])]
file.write('Timoshenko gradient constant = %6.3f' %C[0])
file.write('l/a = %6.3f' %C[1])
file.write('Resulting parameters\n')
def func_b(h, g, EE):#bar((h*C[1])/g*force/(2*b*h*math.cosh(g/(h*C[1]))))*(math.exp(g/(h*C[1]))-math.exp(g/(h*C[1])))+force/(b*h)
    h = np.asanyarray(h)
    return EE/((((h*C[1])/g)/(2*np.cosh(g/(h*C[1]))))*(np.exp(-g/(h*C[1]))-np.exp(g/(h*C[1]))+1))
popt_b, pcov_b = curve_fit(func_b, np.asanyarray(H_total), np.asanyarray(D_b))
Dbar = [d/popt_b[1] for d in D_b]

labelEc_b ='Results Bar Ec = ' + '%.3f' %popt_b[1] + 'MPa'
pl.plot(H_total, Dbar, 'bx', label=labelEc_b)
hint = np.arange(H_total[0],H_total[-1],1.0) ##hint = np.arange(H_total[0],1000,1)
pl.plot(hint, func_b(hint, *popt_b)/popt_b[1], 'g:',label='fit: g=%.3f' % popt_b[0] + ' & Ec=%.3f MPa (BG)' % popt_b[1], lw=1)
pl.xlabel('total height [mm]')
pl.ylabel('D/D_0')
pl.legend()
pl.ylabel('$D/D_{0}$')
pl.autoscale()
pl.title(be_dictionary[base_element])
pl.savefig('NBD_bar' +be_dictionary[base_element] +'.pdf')
pl.show()

# Get back to the working directory
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
#mesh('line', tag, numnodes, *ndtags, id, ndf, meshsize, eleType=”, *eleArgs=[])
#mesh('tri', tag, numlines, *ltags, id, ndf, meshsize, eleType=”, *eleArgs=[])
#element('ElasticTimoshenkoBeam', eleTag, *eleNodes, E, G, A, Iz, Avy, transfTag[, '-mass', massDens ][,'-cMass ]) #2D
#
#element('Tri31', eleTag, *eleNodes, thick, type, matTag[, pressure, rho, b1, b2 ])
#nDMaterial('ElasticIsotropic', matTag, E, v, rho=0.0)
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
be_dictionary ={0:'Line', 1:'Triangle', 2:'DoubleTriangle', 3:'SemiHoneycomb', 
                4:'Rectangle', 5:'SixTetCube', 6:'Mixedcube', 11:'SymTriangle'}
base_element = 1
t=0.5
mini=t/2
maxi=t
shape=''

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
folder = './' + be_dictionary[base_element]+ '_homo' + '/'
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
E = 2000*(10**-6); # [MN/mm2]
nu = 0.25
G = E/(2*(1-nu)); # [MN/mm2]
rho = 1040*10**-6 # [Kg/mm3] weight

## generate a cross section inside
b = 0.5 #[mm]
h = 0.5 #[mm2]
A = b*h # [mm2]
Iyy = (b*h**3)/12 # [mm4] moment of inertia about z axis


# generata the cross section outside (upper&bottom edges)
Aboundary = b*h/2
Iyyboundary = (b*(h/2)**3)/12

# longitudinal dimensions of each beam bottom
L = 5#2.5 # [mm]
H = 5*math.cos(math.pi/6) #[mm]

# load case:
force = 1*10**(-9) # [MN]

###########----> File writting
file = open('Results.txt','w') 
filedisp = open('Displacements.txt','w') 
filedisp.write('DISPLACEMENTS (u,w) at each node [mm]\n\r') 
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


file.write('-  Load:\n') 
file.write('      Point load, f = %E [MN]\n\r' % force)
file.write('\n\r') 
file.write('###################--STARTING--ITERATION--###################\n\r')  
rcParams['font.family'] = 'serif'
rcParams['font.sans-serif'] = ['serif']

# Starting the iteration

D = []
H_total = []
i = -1
eleArgs=['Tri31',1.0, 'PlaneStress',1] #meshargs
for i in range(2):
    # Start timer
    start_time = timeit.default_timer()

    op.wipe()
    op.model('basic', '-ndm', 2,'-ndf', 3)
    op.geomTransf('Linear', 1)
    # Call domain generation function
    mesh, e_i, faces_i, B, n_G = base.main(base_element, shape, i, meshSize=[mini,maxi], sizeXYZ=(L,0,H), t=0.5, ExtractGeom=True, View=False)
    op.nDMaterial('ElasticIsotropic', 1, E, nu, rho)

    # nodes of the domain
    ii = -1
    for node in mesh:
        ii = ii +1
        op.node(ii, *mesh[ii,[0,2]])

    # generate triangle elements from face
    iii = -1


    for n0 in faces_i:
#        ntag.add(n0[0])
#        for n1 in n0[1:]:
#            ntag.add(n1)
#            op.mesh('line', ii, 2, n0[0], n1, 1,3,1.0)
#        ii = ii +1
#        op.mesh('line', ii, 2, n0[0], n1, 1,3,1.0)
        n0.sort(reverse=True)
        iii= iii +1
        op.element('Tri31', iii, *n0, 1.0, 'PlaneStress',1)

    # generate boundary conditions
    for node in B[4]+B[5]: 
        op.fix(node,0,1,0) #-> NO vertical movement
    for node in B[0]:
        op.setNodeDisp(node, 0, -L/2)

    # set up analysis procedure 
    op.constraints('Plain')
    op.numberer('Plain')
    op.test('NormUnbalance', 1e-9, 10)
    op.system('FullGeneral')
    op.integrator('LoadControl',1.0)
    op.algorithm('Linear')
    op.analysis('Static')
    op.analyze(1)
     
    # calculate member forces and displacements
    u0, w0 ,rot0 = 0, 0, 0
    filedisp.write('Iteration %s\n' % i) 
    for node in B[1]:
#        filedisp.write('%s; ' % j) 
#        filedisp.write('%s; %s;\n'%(nodelist[j].getDisplacement(0).x(), nodelist[j].getDisplacement(0).z()))
        u, w, rot = op.nodeDisp(node)
        u0 = u + u0
        w0 = w + w0
        rot0 = rot + rot0
    del(u,w,rot)
    u = u0/len(B[1])
    w = w0/len(B[1])
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
    file.write('RESULTS: \n') 

    print('Average displacements at x = L: \n')
    print("  ux= %.16f[mm]" % u)
    print("  uz= %.16f[mm]" % w) 
    file.write('   average displacement at x = L: \n') 
    file.write('     ux= %.16f[mm]\n' % u) 
    file.write('     uz= %.16f[mm]\n' % w) 
    # Opensees file model creation
    op.printModel('-file', 'Results_opensees', '-JSON', '-node')
    
    # D coef:
    leng = L*n_G[0]
    height = H*i
    Iyy_total = 0.5*(height**3)/12
    D.append([force/w,3*Ec*Iyy_total*(leng**-3)])
    H_total.append(i*H)
    # Finish timer
    elapsed = timeit.default_timer() - start_time
    print("elapsed time: %E seconds" % elapsed )
    file.write('elapsed time: %E seconds \n' % elapsed) 
    file.write('############################################################## \n\r') 
    filedisp.write('############################################################## \n\r') 
    del( B, elapsed)
    

# Write Curve Data to Results.txt file
file.write('D/D_0:\n')
div = []
g=[]
i=-1
for d,do in D:
    i = i+1
    division = d/do
    div.append(division)
    g.append(((division-1)*((1/12.0)*float(H_total[i])**2))**(0.5))
    file.write('%s' %division)
    file.write('%s;' %g[i])
file.write('\n')
file.write('Total heights [mm]:\n')
for t in H_total:
    file.write('%s ;' %t)
file.write('\n')
# Close writting files
file.close() 
filedisp.close() 

# Fit curve
def func(h, g):
    h = np.asanyarray(h)
    return 1+12*(g*g)/(h*h)

popt, pcov = curve_fit(func, np.asanyarray(H_total), np.asanyarray(div))

# Plotting resulting curve
pl.plot(H_total, div, 'ro-', label='simulation')
pl.plot(H_total, func(H_total, *popt), 'g--',label='fit: g=%5.3f' % popt[0])
pl.xlabel('total height of the sandwich')
pl.ylabel('D/D_0')
pl.legend()
#pl.plot(H_total, [1 for d in D],color='grey',ls='--', label='', lw=1)
#pl.plot(H_total, [d for d in div], color='r', ls='-',marker='o', label='Lattice', lw=1, )
#pl.xlabel('total height of the sandwich beam [mm]')
#pl.ylabel('D/D_0')
pl.autoscale()
pl.title(be_dictionary[base_element])
pl.savefig('NBR_' + be_dictionary[base_element] +'.pdf')
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
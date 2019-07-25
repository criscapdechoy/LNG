#!/usr/bin/env python

# -----------------------------------------------------------------------------
# simplysupportedbeam.py
# -----------------------------------------------------------------------------
import sys
import karamba
import base
import timeit
import math
import os
import numpy as np
from matplotlib import rcParams
import matplotlib.pyplot as pl
from scipy.optimize import curve_fit

# Set current direction
parent_folder = os.getcwd()

#################___________CHOOSE-> BASE ELEMENT___________###################
be_dictionary ={0:'Line', 1:'Triangle', 2:'DoubleTriangle', 3:'SemiHoneycomb', 
                4:'Rectangle', 5:'SixTetCube', 6:'Mixedcube', 11:'SymTriangle'}
base_element = 3
shape=''
#################----> Creating results folder
def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)
# Creates a folder in the current directory called data
# load the license
lic = karamba.License.Instance()
lic.loadLicense("LicensePRO.lic", "public.key")
if (lic.error_flg()):
    sys.exit(lic.error_msg()) 
# Jump inside the results folder
folder = './' + be_dictionary[base_element] + 'slim/'
createFolder(folder)
os.chdir(folder)
#############___________MODEL_KARAMBA3D___________#############################
# set number of load cases
lc_num= 1
model = karamba.Model(lc_num)

# generate a material
Ecc = 12.35*(10**-6); # [MN/mm2] Core Young's modulus
Ec = 12.3*(10**-6); # [MN/mm2] Core Young's modulus
E = 2000.0*(10**-6); # [MN/mm2]
nu = 0.25
G = E/(2*(1-nu)); # [MN/mm2]
gamma = 1040.0*9.81*10**-15 # [MN/mm3] weight
material = karamba.Material(E, G, gamma)
material.thisown = False
model.add(material)

# generate a cross section inside
b = 1
h = 0.5
A = b*h # [mm2]
Iyy = (b*h**3)/12 # [mm4] moment of inertia about z axis
Izz = 0.002 # [mm4] moment of inertia about y axis
Ipp = 0.003 # [mm4] torsional moment of inertia
ky  = 0.0 # rigid in shear
kz  = 0.0 # rigid in shear
crosec = karamba.Beam3DCroSec(A, Iyy, Izz, Ipp, ky, kz)
crosec.thisown = False
model.add(crosec)

# generata the cross section outside (upper&bottom edges)
Aboundary = b*h/2
Iyyboundary = (b*(h/2)**3)/12
crosec_bound = karamba.Beam3DCroSec(Aboundary, Iyyboundary, Izz, Ipp, ky, kz)
crosec_bound.thisown = False
model.add(crosec_bound)

# longitudinal dimensions of each beam
L = 25.0 # [mm]
H = 50.0*math.cos(math.pi/6.) #[mm]

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
file.write('      gamma %E [MN/mm3]\n' % gamma)
file.write('-  Longitudinal dimension data:\n')
file.write('      Length of the base element, L = %E [mm]\n' % L)
file.write('      Height of the base element, H = %E [mm]\n' % H)
file.write('-  Cross section data:\n') 
file.write('      b = %E [mm]\n' % b)
file.write('      h = %E [mm]\n' % h)
file.write('      A = %E [mm2]\n' % A)
file.write('      Iyy = %E [mm4]\n' % Iyy)
file.write('      Izz = %E [mm4]\n' % Izz)
file.write('      Ipp = %E [mm4]\n' % Ipp)
file.write('      ky = %E \n' % ky)
file.write('      kz = %E \n' % kz)
file.write('-  Cross section data bread beams:\n')
file.write('      A_bread = %E [mm2]\n' % Aboundary)
file.write('      Iyy_bread = %E [mm4]\n' % Iyyboundary) 
file.write('-  Load:\n') 
file.write('      Point load, f = %E [MN]\n\r' % force) 
rcParams['font.family'] = 'serif'
rcParams['font.sans-serif'] = ['serif']

# Starting the iteration
nlayers = list(range(1,2))
D = []
H_total = []
file.write('\n\r') 
file.write('###################--STARTING--ITERATION--###################\n\r') 
for i in nlayers:
    # Start timer
    start_time = timeit.default_timer()
    print ("---")
    print ("Iteration:"+ str(i))
    print ("---")
    file.write('Iteration: %s\n' % i) 
    
    # Call domain generation function
    m, e_i, faces_i, B, n_G = base.main(base_element, shape, i, meshSize=[h/2, h*(2**0.5)/2], sizeXYZ=(L,0,H), t=0, ExtractGeom=False, View=False)
    nX_G = n_G[0]
    # Initializing list of nodes and list of edges
    nodelist=list()
    elemlist=list()
    
    # Dimensions of the sandwhich beam
    file.write('   Total lenght of the sandwhich beam: %.2f [mm]\n' % float(nX_G*L)) 
    file.write('   Total height of the sandwhich beam: %.2f [mm]\n' % float(i*H)) 
    lc_num= 1
    model = karamba.Model(lc_num)
    crosec_bound = karamba.Beam3DCroSec(Aboundary, Iyyboundary, Izz, Ipp, ky, kz)
    crosec_bound.thisown = False
    model.add(crosec_bound)
    crosec = karamba.Beam3DCroSec(A, Iyy, Izz, Ipp, ky, kz)
    crosec.thisown = False
    model.add(crosec)
    material = karamba.Material(E, G, gamma)
    material.thisown = False
    model.add(material)
    # nodes of the domain
    for c1,c2,c3 in m:
        no=karamba.Node(c1,c2,c3)
        no.thisown = False
        nodelist.append(no)
        model.add(no)
    # generate beam elements from edges inside

    for n0 in e_i:
        for n1 in n0[1:]:
            if not(((n0[0] in B[4])and(n1 in B[4]))or((n0[0] in B[5])and(n1 in B[5]))) :
                element = karamba.Beam3D(nodelist[n0[0]], nodelist[n1], material, crosec)
                element.thisown = False
                model.add(element)
                elemlist.append(element)
    # generate beam elements from edges conforming the bread
    for BB in B[4:]:
        for Bb in BB:
            for Bcon in e_i[Bb][1:]:
                if Bcon in BB:

                    element = karamba.Beam3D(nodelist[e_i[Bb][0]], nodelist[Bcon], material, crosec_bound)
                    element.thisown = False
                    model.add(element)
                    elemlist.append(element)
                
    del(c1,c2,c3,n0,n1)
    # generate boundary conditions
    movement = 0
    for node in B[0]: 
        bc = karamba.BoundaryCondition(node, karamba.Node.x_t, karamba.BoundaryCondition.disp, movement)
        bc.thisown = False
        model.add(bc)
        bc = karamba.BoundaryCondition(node, karamba.Node.y_t, karamba.BoundaryCondition.disp, movement)
        bc.thisown = False
        model.add(bc)
        bc = karamba.BoundaryCondition(node, karamba.Node.z_t, karamba.BoundaryCondition.disp, movement)
        bc.thisown = False
        model.add(bc)
        bc = karamba.BoundaryCondition(node, karamba.Node.x_r, karamba.BoundaryCondition.disp, movement)
        bc.thisown = False
        model.add(bc)
        bc = karamba.BoundaryCondition(node, karamba.Node.y_r, karamba.BoundaryCondition.disp, movement)
        bc.thisown = False
        model.add(bc)
        bc = karamba.BoundaryCondition(node, karamba.Node.z_r, karamba.BoundaryCondition.disp, movement)
        bc.thisown = False
        model.add(bc)
    # Load
    force_p = force/len(B[1]) # divide load to be equally distributed at each right node
    for node in B[1]:
        bc = karamba.BoundaryCondition(node, karamba.Node.z_t, karamba.BoundaryCondition.force, force_p)
        bc.thisown = False
        model.add(bc)
    
        # set up analysis procedure 
    analysis = karamba.Deform(model)
    response = karamba.Response(analysis)
    
    # calculate member forces and displacements
    response.updateMemberForces()
    lc_ind = 0
    
    # maximum deflection in load case lc_ind
    maxDisp = response.maxDisplacement(lc_ind)
    print("---")
    print("results")
    print("---")
    
    print("maximum displacement %.16f [mm]" % maxDisp )
    file.write('RESULTS: \n') 
    file.write('   maximum displacement %.16f [mm]\n'% maxDisp) 
    u, w = 0, 0
    filedisp.write('Iteration %s\n' % i) 
    for node in B[1]:
#        filedisp.write('%s; ' % j) 
#        filedisp.write('%s; %s;\n'%(nodelist[j].getDisplacement(0).x(), nodelist[j].getDisplacement(0).z()))
        u = u+float(nodelist[node].getDisplacement(0).x())
        w = w+float(nodelist[node].getDisplacement(0).z())
    u = u/len(B[1])
    w = w/len(B[1])
    print('Average displacements at x = L: \n')
    print("  ux= %.16f[mm]" % u)
    print("  uz= %.16f[mm]" % w) 
    file.write('   average displacement at x = L: \n') 
    file.write('     ux= %.16f[mm]\n' % u) 
    file.write('     uz= %.16f[mm]\n' % w) 
    
    # D coef:
    leng = L*nX_G
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
    del( elapsed, elemlist, m, w, u, maxDisp, nodelist, model)

# Write Curve Data to Results.txt file
file.write('D/D_0:\n')
div = []
g=[]
divv=[]
i=-1
for d,do in D:
    i = i+1
    division = d/do
    div.append(division)
    divv.append((division*(Ec/Ecc)))
    #do=2*3*Ec*Iyy_total*(leng**-3)
    g.append((abs((division-1))*((1/12.0)*float(H_total[i])**2))**(0.5))
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
gg = 1.4 #g to draw the curve by hand
divvv=[1+12*gg**2/(h**2) for h in H_total]
pl.plot(H_total, div, 'ro', label='simulation Ec=246.7 MPa')
pl.plot(H_total, divv, 'bo', label='simulation Ec=262.5 MPa')
pl.plot(H_total, divvv, 'c-', label='drawing: g=%5.3f' % gg)
pl.plot(H_total, func(H_total, *popt), 'g--',label='fit: g=%5.3f' % popt[0])
pl.xlabel('total height of the sandwich (mm)')
pl.ylabel('D/D_0')
pl.legend()
pl.plot(H_total, [1 for d in D],color='grey',ls='--', label='', lw=1)
#pl.plot(H_total, [d for d in div], color='r', ls='-',marker='o', label='Lattice', lw=1, )
#pl.xlabel('total height of the sandwich beam [mm]')
#pl.ylabel('D/D_0')
pl.xlim(0,H_total[-1])
pl.title(be_dictionary[base_element])
pl.savefig('NBR_' + be_dictionary[base_element] +'.pdf')
pl.show()

# Get back to the working directory
os.chdir(parent_folder)
print('------------------------------END---------------------------------')

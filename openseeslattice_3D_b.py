# -*- coding: utf-8 -*-

"""

Created on Sun Apr 28 18:41:04 2019



@author: Cristina

"""
import openseespy.opensees as op
import LNG_builder
import LNG_engine
import LNG_main
import timeit
import math
import os
import numpy as np
from matplotlib import rcParams, rc
import matplotlib.pyplot as pl
from scipy.optimize import curve_fit
import copy
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib as mpl
from dxfwrite import DXFEngine as dxf
#################___________CHOOSE-> BASE ELEMENT___________###################

be_dictionary ={0:'Line', 1:'Triangle', 2:'DoubleTriangle', 3:'SemiHoneycomb',

                4:'Rectangle', 5:'Cube', 6:'X', 7:'hexa_tetra', 8:'octet', 9:'FCC',10:'CFCC',11:'BFCC'}
base_element = 8
shape=''
nlayers = list(range(1,2)) # number of layers to test
h = 0.10663773 #[mm]
siderep = 10
#################----> Creating results folder
# Set current direction
parent_folder = os.getcwd()
# Jump inside the results folder
folder = './'+'plate_8'+str(siderep) + be_dictionary[base_element] + 'h' +str(h)+ 'it_'+ str(nlayers[-1])+'_3D/'
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
Ec = 9.43 # [MN/mm2] Core Young's modulus
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
###########----> File writting
file = open('Results.txt','w')

#filedisp = open('Displacements.txt','w')
#filedisp.write('DISPLACEMENTS (u,w) at each node [mm]\n\r')
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
file.write('-  Cross section data bread beams:\n')

#file.write('      A_bread = %E [mm2]\n' % Aboundary)

#file.write('      Iyy_bread = %E [mm4]\n' % Iyyboundary)

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
B_total = []
for ilay in nlayers:
    ilay=3
    # Start timer
    start_time = timeit.default_timer()
    op.wipe()
    op.model('basic', '-ndm', 3,'-ndf', 6)
#    op.geomTransf('Linear', 1, 1.0,1.0,0.0) # X beam

    op.geomTransf('Linear', 1, 0.0,0.0,1.0) # Z beam
    op.geomTransf('Linear', 2, 0.0,1.0,0.0) # Y beam
    n_G = [ilay*2*10, 2*siderep, 2*ilay]
    # Call domain generation function base_element, shape, [18*i,1,i],
    LS = LNG_builder.main(base_element, shape, n_G, meshSize=[], sizeXYZ=[L,L,L], Export=False, View=False)
    # nodes of the domain
    # generate beam elements from edges inside
    ii = -1
    NODE=set()
    for n0 in LS.edges:

        if n0[0] not in NODE:
            NODE.add(n0[0])
            op.node(n0[0], *LS.nodes.astype('float64')[n0[0],:])
        for n1 in n0[1:]:
            if n1 not in NODE:
                NODE.add(n1)
                op.node(n1, *LS.nodes.astype('float64')[n1,:])
            transfTag = 2 - (sum([LS.nodes[n0[0],o]-LS.nodes[n1,o] for o in [0,1]]) != 0) # = 1 if not column or = 2 for column
            ii= ii+1
            if not(((n0[0] in LS.boundary[4])and(n1 in LS.boundary[4]))or((n0[0] in LS.boundary[5])and(n1 in LS.boundary[5]))or(

                    (n0[0] in LS.boundary[0])and(n1 in LS.boundary[0]))or((n0[0] in LS.boundary[1])and(n1 in LS.boundary[1]))or(

                    (n0[0] in LS.boundary[2])and(n1 in LS.boundary[2]))or((n0[0] in LS.boundary[3])and(n1 in LS.boundary[3]))):

                op.element('elasticBeamColumn', ii, n0[0],n1, A, E, G, J, Iyy,Izz, int(transfTag))#int(transfTag)) #3D
            else:
                op.element('elasticBeamColumn', ii, n0[0],n1, Abound, E, G, Jbound, Iyybound,Izzbound, int(transfTag))#int(transfTag)) #3D
    del(LS.nodes, NODE, LS.edges)
    # generate beam elements from edges conforming the boundary
    # generate boundary conditions
    for node in LS.boundary[0]:

        op.fix(node, 1,1,1,1,1,1) #0free 1fix

    # Load

    force_p = force/len(LS.boundary[1]) # divide load to be equally distributed at each right node

    op.timeSeries('Constant',1)

    op.pattern('Plain',1,1)

    i=0

    while i < len(LS.boundary[1])-1:

        op.load(LS.boundary[1][i], 0.0,0.0,force_p,0.0,0.0,0.0)

#        op.rigidDiaphragm('beam', LS.boundary[1][-1],  LS.boundary[1][i])

        i = i+1

    op.load(LS.boundary[1][-1], 0.0,0.0,force_p,0.0,0.0,0.0)

#    op.rigidDiaphragm(1, LS.boundary[1][-1],  *LS.boundary[1][:-1])

    bound = LS.boundary
    del(LS)
    # set up analysis procedure

    op.constraints('Transformation')

    op.numberer('RCM')

 #   op.test('NormUnbalance', 1e-20, 100)

    op.system('BandGen')

    op.integrator('LoadControl',1.0)

    op.algorithm('Linear')

    op.analysis('Static')

    op.analyze(1)



    # calculate member forces and displacements

    u0, v0, w0 = 0.0, 0.0, 0.0

    #filedisp.write('Iteration %s\n' % i)

    for node in bound[1]:

#        #filedisp.write('%s; ' % j)

#        #filedisp.write('%s; %s;\n'%(nodelist[j].getDisplacement(0).x(), nodelist[j].getDisplacement(0).z()))

        U = op.nodeDisp(node)

        u0 = U[0] + u0

        v0 = U[1] + v0

        w0 = U[2] + w0

    del(U)

    u = u0/len(bound[1])

    v = u0/len(bound[1])

    w = w0/len(bound[1])

    # Print in console and file

    print ("---")

    print ("Iteration:"+ str(ilay))

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

    print("  uz= %.16f[mm]" % v)

    print("  uz= %.16f[mm]" % w)

    file.write('   average displacement at x = L: \n')

    file.write('     ux= %.16f[mm]\n' % u)

    file.write('     uz= %.16f[mm]\n' % w)

    # Opensees file model creation





    # D coef:

    leng = float(L*n_G[0])

    W.append(w)

   # DD.append(4.0*100.0*force/(w*Ec*b))

    H_total.append(n_G[2]*H)

    B_total.append(n_G[1]*H)



#    H_total.append(i*H)

    # Finish timer

    elapsed = timeit.default_timer() - start_time

    print("elapsed time: %E seconds" % elapsed )

    file.write('elapsed time: %E seconds \n' % elapsed)

    file.write('############################################################## \n\r')

    #filedisp.write('############################################################## \n\r')

    op.wipeAnalysis()

    del( bound, elapsed)



#filedisp.close()

# Write Curve Data to Results.txt file

file.write('\n')

file.write('[height        [mm]; side       [mm] average_displ [mm]; D/E             ]\n')



for i in range(len(H_total)):

    file.write('%s ;' %H_total[i] )

    file.write('%s ;' %H_total[i] )

    file.write('%s ;' %W[i])

    D.append(4.0*(float(leng/H_total[-1])**3)*force/(W[i]*B_total[-1]))

    file.write('%s ]\n' %D[i])

file.write('Curve fitting with least squares: \n')

# Fitting curve preparation

C = [G*ky*b*H_total[-1], (leng/H_total[-1])]

file.write('Timoshenko gradient constant = %6.3f' %C[0])

file.write('l/a = %6.3f' %C[1])

file.write('Resulting parameters\n')

file.close()
#
##H_total[:]=H_total[1:]
##D[:]=D[1:]
## Fitting resulting curve
#def funcT(h, g): #Timoshenko
#    h = np.asanyarray(h)
#    return  (1+12*(g*g)/(h*h))
#
#def func(h, g, EE): #Euler-Bernoulli
#    h = np.asanyarray(h)
#    return (1+12*(g*g)/(h*h))*EE
## Fitting resulting curve
##poptT, pcovT = curve_fit(funcT, np.asanyarray(H_total[:]), np.asanyarray(D[:]))
#
#popt, pcov = curve_fit(func, np.asanyarray(H_total[:]), np.asanyarray(D[:]))
#poptT, pcovT = curve_fit(funcT, np.asanyarray(H_total[:]), np.asanyarray(D[:])/Ec)
#
#Da = [d/Ec for d in D]
#Db = [d/popt[1] for d in D]
## Plotting
#labelEcB = 'Result Lattice'
#fig = pl.figure()
#hint = np.arange(H_total[0],H_total[-1],0.2) ##hint = np.arange(H_total[0],1000,1)
#pl.semilogx([],[],alpha=0.0,label='E$^{homo}_{c}$=%.1f MPa' % popt[1])
#pl.semilogx(H_total, [1 for d in D],color='grey',ls='--', label='Classical Linear Theory', lw=1)
#pl.semilogx(hint, (func(hint, *popt)/popt[1]), 'g:',label='Gradient EB: g$^{fit}$=%.3f' % popt[0], lw=1)
#pl.semilogx(H_total, Db, 'ro', label=labelEcB)
#pl.xlabel('total height [mm]')
#pl.ylabel('$D/D_{0}$')
#pl.legend()
#pl.autoscale()
#pl.show()
#pl.savefig('NBD_' +be_dictionary[base_element] +'a.pdf')
#pl.close()
#
#gg=np.linspace(poptT[0]*0.85,poptT[0]*1.15,12)
#fig, (ax1, ax2) = pl.subplots(1,2, sharex=False,sharey=False,gridspec_kw={ 'width_ratios':[0.95, 0.05]})
##fig.subplots_adjust(left=0.6)
#cmap = mpl.cm.cool
#norm = mpl.colors.Normalize(vmin=min(gg), vmax=max(gg))
#cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=cmap,
#                                norm=norm,
#                                orientation='vertical')
#ggN=(np.max(gg)-gg)/(np.max(gg)-np.min(gg))
#ax2.set_title('g')
#
#for i in range(len(gg)):
#    ax1.plot(hint, (funcT(hint, gg[i])), c=cmap(ggN[i]),lw=0.85)
#ax1.semilogx([],[],alpha=0.0,label='E$^{homo}_{c}$=%.1f MPa' % Ec)
#ax1.semilogx(H_total, Da, 'kx', label=labelEcB)
#ax1.semilogx(H_total, [1 for d in D],color='grey',ls='--', label='Classical Linear Theory', lw=1)
#ax1.semilogx(hint, (funcT(hint, *poptT)), ':',color='black',label='Gradient EB: g$^{fit}$=%.3f' % poptT[0] , lw=1.5)#ax1.xlabel('total height [mm]')
#ax1.set_xlabel('total height [mm]')
#ax1.set_ylabel('$D/D_{0}$')
#ax1.legend()
#ax1.autoscale()
#fig.show()
#pl.savefig('NBD_' +be_dictionary[base_element] +'b.pdf')
## Get back to the working directory
os.chdir(parent_folder)


print('------------------------------END---------------------------------')
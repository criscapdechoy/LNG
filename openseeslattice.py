#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 18:47:00 2019

@author: cristinacapdevilachoy
"""

#SISIIIII
import openseespy.opensees as op
import LNG_builder
import LNG_main
import timeit
import math
import LNG_engine
import os
import numpy as np
from matplotlib import rcParams, rc
import matplotlib.pyplot as pl
from scipy.optimize import curve_fit
import matplotlib as mpl
#################___________CHOOSE-> BASE ELEMENT___________###################
be_dictionary ={0:'Line', 1:'Triangle', 2:'DoubleTriangle', 3:'SemiHoneycomb',
                4:'Rectangle', 5:'SixTetCube', 6:'Mixedcube', 11:'SymTriangle'}
base_element = 1
nlayers = list(range(1, 25)) # number of layers to test


#################----> Creating results folder
# Set current direction
parent_folder = os.getcwd()

# Creates a folder in the current directory called data
shape=''

# Jump inside the results folder
h = 0.5#[mm]
folder = './' + be_dictionary[base_element] + shape +'_T_h' + str(h)+'mm_it'+str(nlayers[-1])+'/'
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
Ec=246.7
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
L = 2.5#math.cos(math.pi/6) # [mm] ,0.0,
H = 5.0*math.cos(math.pi/6) #[mm]

# load case:
force = 1.0 # [MN]

###########----> File writting###############
file = open('Results.txt','w')
filedisp = open('Displacements.txt','w')
filedisp.write('DISPLACEMENTS (u,w) at each node [mm]\n\r')
file.write( str(base_element) + '_shape_' + shape + '_height_' + str(h) + '[mm] \n')
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

for i in nlayers:
    n_G = [18*2*i,1,i]
    # Start timer
    start_time = timeit.default_timer()
    op.wipe() #cler model
    op.model('basic', '-ndm', 2,'-ndf', 3) #initialize model
    op.geomTransf('Linear', 1) #define geometric transformation
    # Call domain generation function
    LS = LNG_builder.main(base_element, shape, n_G, sizeXYZ=[L,0,H],  Export = False, View = False)
  #  op.nDMaterial('ElasticIsotropic', 1, E, nu, rho)

    # nodes of the domain
    ii = -1
    ilayer = 0
    ii0 = 0

    for node in LS.nodes:
        ii = ii +1
        op.node(ii, *node.astype('float64')[[0,2]])
    totN = ii
#        try:
#            if mesh[ii+1,2] > mesh[ii,2]:
#                ii0 = ii0 + 1 #int((ii0 in  LS.boundary[0]))
#                print(ii0)
#                for iii in range(ii0, ii):
#
#                    op.equalDOF(ii, iii,*[1,2])
#                ii0 = ii+1
#                print(ii0)
#        except:
#            pass
    del(LS.nodes)
    # generate beam elements from edges inside
    ii= -1
    tot = -1
    for n0 in LS.edges:
        for n1 in n0[1:]:
            tot = tot+1
            if not(((n0[0] in LS.boundary[4])and(n1 in LS.boundary[4]))or((n0[0] in LS.boundary[5])and(n1 in LS.boundary[5]))) :
                ii = ii+1
                op.element('elasticBeamColumn', ii, n0[0],n1, A, E, Iyy, 1) #2D
#                op.element('ElasticTimoshenkoBeam', ii, n0[0],n1, E, G, A, Iyy, A*ky, 1) #2D
#    print(tot, ii)

   # generate beam elements from edges belonging to the boundary
    for BB in LS.boundary[4:]:
        for Bb in BB:
            for Bcon in LS.edges[Bb][1:]:
                if Bcon in BB:
                    ii = ii +1
                    op.element('elasticBeamColumn', ii, LS.edges[Bb][0], Bcon, Aboundary, E, Iyyboundary, 1) #2D
#                    op.element('ElasticTimoshenkoBeam', ii, LS.edges[Bb][0], Bcon, E, G, Aboundary, Iyyboundary, Aboundary*ky, 1) #2D
#    print(ii)
    del(LS.edges)

    # generate boundary conditions
    for node in LS.boundary[0]:
        op.fix(node, 1,1,1) #0free 1fix
        modB = divmod(len(LS.boundary[1]),2)
#    if modB[1]!=0:
#        op.fix(LS.boundary[1][modB[0]], 1,0,0) #0free 1fix
    # Load
    force_p = float(float(force)/float(len(LS.boundary[1])))
    # divide load to be equally distributed at each right node
    op.timeSeries('Constant',1)
    op.pattern('Plain',1,1)
    k = 0
    for k in range(len(LS.boundary[1])-1):
        op.load(LS.boundary[1][k], 0.0 ,force_p, 0.0)
#        op.equalDOF(LS.boundary[1][-1],  LS.boundary[1][k], *[2,3])
        op.rigidLink('beam', LS.boundary[1][-1],  LS.boundary[1][k])
    op.rigidLink('beam', LS.boundary[1][-1],  LS.boundary[1][0])

    op.load(LS.boundary[1][-1], 0.0,force_p,0.0)
##        if k != modB[0]:
#            op.equalDOF(LS.boundary[1][modB[0]],  LS.boundary[1][k], *[2,3])
##        k = k+1
##    op.rigidLink('bar', LS.boundary[1][-1],  LS.boundary[1][0])
#    LS.boundary[0].sort()
#    nlayer = LS.boundary[0][1]-LS.boundary[0][0]
#    for k in range(1,len(LS.boundary[4])):
#        for kk in range(1,i+1):
#            op.equalDOF(LS.boundary[4][k],  LS.boundary[4][k]+kk*nlayer, *[2,3])
##

#    op.load(LS.boundary[1][-1], 0.0,force_p,0.0)
#     set up analysis procedure
    op.constraints('Transformation')
    op.numberer('RCM')
 #   op.test('NormUnbalance', 1e-20, 100)
    op.system('BandGen')
    op.integrator('LoadControl',1.0)
    op.algorithm('Linear')
    op.analysis('Static')
    op.analyze(1)

    # calculate member forces and displacements
    u0, w0 ,rot0 = 0.0, 0.0, 0.0
#    filedisp.write('Iteration %s\n' % i)
    for node in range(totN+1):
        u, w, rot = op.nodeDisp(node)
        u00, w00 = op.nodeCoord(node)
        u00= u+ u00
        w00 = w+ w00
        filedisp.write('%s ;' %u )
        filedisp.write('%s ;' %w)
        filedisp.write('%s ;' %rot)
        filedisp.write('%s ; ' %u00)
        filedisp.write('%s ;\n' %w00)

        if node in LS.boundary[1]:
            u0 = float(u + u0)
            w0 = float(w + w0)
            rot0 = float(rot + rot0)
    u = u0/float(len(LS.boundary[1]))
    w = w0/float(len(LS.boundary[1]))
    rot = rot0/len(LS.boundary[1])
    filedisp.write('##########################################')

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
    del( LS.boundary, elapsed, i, ii)

filedisp.close()
# Write Curve Data to Results.txt file
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
file.write('l/a = %6.3f' %C[1])
file.write('Resulting parameters\n')
file.close()

def funcT(h, g): #Timoshenko
    h = np.asanyarray(h)
    return  (1+12*(g*g)/(h*h))

def func(h, g, EE): #Euler-Bernoulli
    h = np.asanyarray(h)
    return (1+12*(g*g)/(h*h))*EE
# Fitting resulting curve
#poptT, pcovT = curve_fit(funcT, np.asanyarray(H_total[:]), np.asanyarray(D[:]))

popt, pcov = curve_fit(func, np.asanyarray(H_total[:]), np.asanyarray(D[:]))
poptT, pcovT = curve_fit(funcT, np.asanyarray(H_total[:]), np.asanyarray(D[:])/Ec)

Da = [d/Ec for d in D]
Db = [d/popt[1] for d in D]
# Writting
#file.write('gB = %6.6f' %poptT[0])
#file.write('gT = %6.6f' %popt[0])
#file.write('EcB = %6.6f' %poptT[1])
#file.write('EcT = %6.6f' %popt[1])

# Plotting
labelEcB = 'Result Lattice'

#labelEcB = 'Gradient Bernoulli Ec = ' + '%.1f' %popt[1] + 'MPa'
#labelEcT ='Results Timoshenko Ec = ' + '%.1f' %poptT[1] + 'MPa'
#pl.plot(H_total, DT, 'bx', label=labelEcT)
#l.plot(H_total, divvv, 'c-', label='drawing: g=%5.3f' % gg)
#pl.plot(H_total, funcT(H_total, *popt)/popt[1], 'g--',label='fit: g=%5.3f' % poptT[0] + ' & Ec=%5.4f MPa' % poptT[1])
hint = np.arange(H_total[0],H_total[-1],0.2) ##hint = np.arange(H_total[0],1000,1)
#pl.plot(hint, funcT(hint, *poptT)/poptT[1], 'g:',label='Gradient T: g=%.3f' % popt[0] + ' & Ec=%.1f MPa' % popt[1], lw=1)
pl.semilogx([],[],alpha=0.0,label='E$^{homo}_{c}$=%.1f MPa' % popt[1])
pl.semilogx(H_total, [1 for d in D],color='grey',ls='--', label='Classical Linear Theory', lw=1)
pl.semilogx(hint, (func(hint, *popt)/popt[1]), 'g:',label='Gradient EB: g$^{fit}$=%.3f' % popt[0], lw=1)
pl.semilogx(H_total, Db, 'ro', label=labelEcB)
pl.xlabel('total height [mm]')
#pl.ylabel('D/D_0')
pl.ylabel('$D/D_{0}$')
pl.legend()
pl.autoscale()
#pl.title(be_dictionary[base_element])
pl.show()
pl.savefig('NBD_' +be_dictionary[base_element] +'a.pdf')
pl.close()

#fig, axs = pl.subplots(1, 1)
#cmap = pl.get_cmap('cool', 5)
#cmap.set_under('gray')
#colors = pl.cm.cool((max(gg)-gg)/(max(gg)-min(gg)))
#for i in range(len(gg)):
#    pl.plot(hint, (funcT(hint, gg[i])), c=gg[i],cmap=cmap,lw=0.85, label='g=%.3f' % gg[i])#,label='Gradient EB (a): g=%.3f' % g + ' & Ec=%.1f MPa' % Ec, lw=1)
##pl.semilogx(hint, (func(hint, *popt)/popt[1]), 'g:',label='Gradient EB (b): g=%.3f' % popt[0] + ' & Ec=%.1f MPa' % popt[1], lw=1)
#pl.semilogx(H_total, Da, 'kx', label=labelEcB)
#pl.xlabel('total height [mm]')
##pl.ylabel('D/D_0')
#pl.ylabel('$D/D_{0}$')
#pl.legend()
#pl.autoscale()
##cax = axs.scatter(x, y, c=z, s=100, cmap=colors, vmin=0.1, vmax=z.max())
#fig.colorbar(cax, extend='min')
#
#
#pl.title('E$^{homo}_{c}$=%.1f MPa' % Ec)
#fig.col
#pl.show()
gg=np.linspace(popt[0]*0.85,popt[0]*1.15,12)
fig, (ax1, ax2) = pl.subplots(1,2, sharex=False,sharey=False,gridspec_kw={ 'width_ratios':[0.95, 0.05]})
#fig.subplots_adjust(left=0.6)
cmap = mpl.cm.cool
norm = mpl.colors.Normalize(vmin=min(gg), vmax=max(gg))
cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=cmap,
                                norm=norm,
                                orientation='vertical')
ggN=(np.max(gg)-gg)/(np.max(gg)-np.min(gg))
ax2.set_title('g')

for i in range(len(gg)):
    ax1.plot(hint, (funcT(hint, gg[i])), c=cmap(ggN[i]),lw=0.85)
ax1.semilogx([],[],alpha=0.0,label='E$^{homo}_{c}$=%.1f MPa' % Ec)
ax1.semilogx(H_total, Da, 'kx', label=labelEcB)
ax1.semilogx(H_total, [1 for d in D],color='grey',ls='--', label='Classical Linear Theory', lw=1)
ax1.semilogx(hint, (funcT(hint, *poptT)), ':',color='black',label='Gradient EB: g$^{fit}$=%.3f' % poptT[0] , lw=1.5)#ax1.xlabel('total height [mm]')
ax1.set_xlabel('total height [mm]')
ax1.set_ylabel('$D/D_{0}$')
ax1.legend()
ax1.autoscale()
fig.show()
pl.savefig('NBD_' +be_dictionary[base_element] +'b.pdf')
# Get back to the working directory
os.chdir(parent_folder)
print('------------------------------END---------------------------------')

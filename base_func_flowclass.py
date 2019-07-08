#!/usr/bin/env python
# -----------------------------------------------------------------------------
# base_func_flowchart.py
# -----------------------------------------------------------------------------
import random 
import numpy as np
import os
import math
from dxfwrite import DXFEngine as dxf
from scipy.spatial import Delaunay , delaunay_plot_2d
from matplotlib import collections as mc
import pylab as pl
import copy
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
#################----> Creating results folder
# Set current direction

class LS:
    def __init__(self, nodes, edges, faces, sym, dim, boundary=[], cros_sec_bound=[]):
        self.nodes = self.data(nodes)
        self.edges = edges
        self.faces = faces
        self.sym = sym
        self.dim = dim
        self.boundary = boundary
#        self.nodeR = nodeR
        self.cros_sec_bound = cros_sec_bound
        self.N = self.length(self)
    class data:
        def __init__(self, nodes):
            self.N=len(nodes)
        
#    def print(self):
#        print("\n mesh: ", self.mesh)
#        print("\n edges: ", self.edges)
#        print("\n faces: ", self.faces )
#        print("\n boundary: ", self.boundary)

def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)
#    print('createFolder( directory F)')

def openFolder(foldername):
    folder = './drawings_/'+ foldername
    createFolder(folder)
    os.chdir(folder)
#    print('openFolder( foldername F)')

def do_modeV( v ):
#    print(' do_modeV( v F) // V**(0.5) G')
    V = 0
    for vi in v:
        V = V + float(vi)**2
    return V**(0.5)

def do_proj_vect(v,w):
#    print(' do_modeV( v F) // proj')
    # projects v over w:
    W = do_modeV(w)**2
    sc = 0
    for i in range(len(v)):
        sc = sc + float(v[i])*float(w[i])
    proj = np.zeros(len(v))
    for i in range(len(v)):
        proj[i] = float((sc/W))*w[i]
    return  proj
    

def get_basic_elem( e ):
    sym = []
    dim = [0,2]
    if e == 0:
        nodes = np.array([[0,0,0],[1,0,0]], dtype="float")
        edges = np.array([[0,1]])
        faces = np.array([[]]) 
        dim = [0]
    elif e == 1:
        nodes = np.array([[1,0],[0,0],[1,1],[0,1]], dtype="float")
        edges = np.array([[0,1],[1,2],[2,3]])
        faces = np.array([[]]) 
        sym = [0,2]    
    elif e == 11:
        nodes = np.array([[0,0],[1,1],[1,0],[0,1]], dtype="float")
        edges = np.array([[0,1],[1,3],[2,0]])
        faces = np.array([[0,1,3],[1,0,2]]) 
        sym = [0,2]
    elif e == 3:
        nodes = np.array([[0,0],[0,0.5],[1,1]], dtype="float")
        edges = np.array([[0,1],[1,2]])
        faces = np.array([[]]) 
        sym = [0,2]
    elif e == 4:
        nodes = np.array([[0,0],[0,1],[1,1],[1,0]], dtype="float")
        edges = np.array([[0,1],[1,2],[2,3],[3,0]])
        faces = np.array([[]]) 
    elif e == 5: # cube
        nodes = np.array([[0,0,0],[1,0,0],[0,1,0],[1,1,0],[0,0,1],[1,0,1],[0,1,1],[1,1,1]], dtype="float")
        edges = np.array([[0,1],[0,2],[2,3],[1,3],[0,4],[1,5],[2,6],[3,7],[4,5],[4,6],[5,7],[6,7]])
        faces = np.array([[0,1,3,2],[4,5,7,6],[0,1,5,4],[0,2,6,4],[1,3,7,5],[2,3,7,6]])
        dim = (0,1,2)
    elif e == 6: # X shape
        nodes = np.array([[0,0,0],[1,0,0],[0,1,0],[1,1,0],[0,0,1],[1,0,1],[0,1,1],[1,1,1]], dtype="float")
        edges = np.array([[0,7],[1,6],[2,5],[3,4]])
        faces = np.array([[]])
        dim = (0,1,2)
    elif e == 7: # hexa-tetrahedrom
        nodes = np.array([[0,0,0],[1,0,0],[0,1,0],[1,1,0],[0,0,1],[1,0,1],[0,1,1],[1,1,1]], dtype="float")
        edges = np.array([[0,1],[0,2],[0,3],[2,3],[1,3],[2,7],[1,7],[1,4],[2,4],[4,6],[4,5],[5,7],[6,7],[4,7]])
        faces = np.array([[0,2,3],[0,1,3],[4,6,7],[4,5,7],[2,4,7],[1,4,7]])
        sym = [0,1,2]
        dim = (0,1,2)
    elif e == 8: # octet-truss
        nodes = np.array([[0,0,0],[1,1,0],[1,0,1],[0,1,1]], dtype="float")
        edges = np.array([[0,1],[0,2],[0,3],[2,3],[1,2],[1,3]])
        faces = np.array([[]])
        sym = [0,1,2]
        dim = (0,1,2)
#    print(' get_basic_elem( e F) // nodes, edges, faces, sym, dim G')
    return nodes, edges, faces, sym, dim

def do_scale(coord, scale):
#    print(' do_scale( coord, scale F) // scaled_coord G')
    if type(coord) == list:
        coord = np.asanyarray(coord)
    max_coord = np.max(coord, axis=0)
    scaled_coord = np.multiply(np.divide(coord,max_coord), scale)
    return scaled_coord

def get_input(coord, edges, faces, sizeXYZ, dim):
    if type(sizeXYZ)!= np.array:
        sizeXYZ = np.asmatrix(sizeXYZ)
    # Amount of nodes
    scale = np.divide(sizeXYZ[:,dim],np.max(coord, axis=0))
    coord = np.multiply(coord,scale)
    N = len(coord)
    # Adding missing dimension in case of node coordinates given in 2D
    if np.shape(coord)[1] != 3:
        coord = np.array([coord[:,0],np.zeros([np.shape(coord)[0],1]),coord[:,1]]).T[0]
        scale = np.array([scale[:,0],np.zeros([np.shape(scale)[0],1]),scale[:,1]]).T[0]
    T_edges = [[a] for a in range(N)]
    for i,j in edges:
        if j not in T_edges[i]:
            T_edges[i].append(j)
#    print('get_input(coord, edges, faces, sizeXYZ, dim F) // coord, T_edges, N, scale G')
    if np.sum(sizeXYZ!=0)==len(dim):
        pass
    else:
        raise ValueError('Element cell size (sizeXYZ) and element dimensions (dim) do not match: sizeXYZ=%s dim=%s'% (sizeXYZ, dim))
        
    return coord, T_edges, faces, N, scale

def check_dimensions(dim, n_G):
    for dimension in [0,1,2]:
        if dimension not in dim:
            n_G[dimension] = 1
#    print('check_dimensions(dim, n_G F) // n_G G')
    return n_G

def get_unique_sorted(column_matrix):
    l_cm = len(column_matrix)
    column_matrix=np.reshape(column_matrix,( l_cm,1))
    index_sort = np.argsort(column_matrix,axis = 0)
    unique_values, uni_sort_ind = np.unique(np.reshape(column_matrix[index_sort]
    ,(l_cm,1)),axis=0,return_inverse = True)
    c_mat_us = uni_sort_ind[np.argsort(index_sort,axis = 0)]
#    print('get_unique_sorted(column_matrix F) // unique_values, c_mat_us G')
    return unique_values, c_mat_us

def check_Dim_Mult(coord, sizeDim):
#    print('check_Dim_Mult(coord, sizeDim F) // coord, delta, idim, n_L G')
    delta, idim = get_unique_sorted(coord)
    if delta[-1] != 0:
        coord = np.matrix(coord*sizeDim/float(delta[-1]))
        delta = float(sizeDim/delta[-1])*delta
    n_L = int(np.shape(delta)[0])
    return coord, delta, idim, n_L
    
def do_translate(tobetranslated, dictionary):
    if type(tobetranslated) == list:
        translation=[]
        word_o = -1
        for word in tobetranslated:
            if word == word_o:
                pass
            else:
                try:
                    translation.append(dictionary[word])
                except:
                    pass
            word_o = word
#    print('do_translate(tobetranslated, dictionary F) // translation G')
    return translation

def check_boundary(num, mesh, p, index, ii, i, j, k, n_L, n_G, Boundary, shape, dim, delmax):
#    print('check_boundary(num_r, mesh_r, ii, i, j, k, n_L, n_G, delmax, Boundary, shape, dim F) // Boundary G')
#    
#    if (i == 0 or n_G[0]) and (j ==0 or n_G[1]) and (k ==0 or n_G[2]):
#        for i in dim:
#            b0 = (num_r[abs(mesh_r[num_r,i] - do_shapeME(np.zeros([1,3]), n_L, n_G, delmax, shape, ii)[0,i])<= 10**-7])
#            for b in b0:
##                flag = True
#                if b not in Boundary[(i*2)]:
#                    Boundary[(i*2)].append(b)
##                    flag = False
#            b1 = (num_r[abs(mesh_r[num_r,i] - np.multiply(n_G, delmax)[i] - do_shapeME(np.array([np.multiply(n_G, delmax)]), 
#                        n_L, n_G, delmax, shape, ii)[0,i])<=10**-7])
##            if b1 
#            for b in b1:
#                if b not in Boundary[(i*2+1)]:
#                    Boundary[(i*2+1)].append(b)
#                   
                        
#    if 1 not in dim:
#        Boundary=[Boundary[0],Boundary[1],[],[],Boundary[4],Boundary[5]]
#    return Boundary
#    num_r = [n for n in num_r[0]]
    dix=[i,j,k]
    flag=True
    ii=np.reshape(np.asarray(ii),(len(ii),3))
    iimax=np.max(ii, axis=0)
    num_r=num[p,:].astype("int")
    for d in dim:
        b0 = num_r[ii[:,d]==0].tolist()
        b1 = num_r[ii[:,d]==iimax[d]].tolist()
        if dix[d] == 0:
            for b in b0:
                if b not in Boundary[(d*2)]:
                    Boundary[2*d].append(b)
        if dix[d] == (n_G[d]-1):
            kk = 1
            symb=Boundary[2*d]
            for b in b1:
                for s in symb:
                    if all(abs(mesh[b,:]-mesh[s,:])**2 <10e-4):
                        index[b]=0
                        num_r[num_r==b]=s
                        b = s
                        kk = kk+1
                if b not in Boundary[(d*2+1)]:               
                    Boundary[2*d+1].append(b)

    return Boundary, num_r, index

def gen_nodegrid(coord, N, sym, n_G, DIM, sizeXYZ, shape, f=[1,1,1]):
#    print('gen_nodegrid(coord, N, sym, n_G, DIM, sizeXYZ, shape,f=[1,1,1] F) // mesh, ind_hash, num, delmax, Boundary G')
    # Get data for the grid generation:
    n_G = check_dimensions( DIM, n_G)
    n_L, ii, delmax, delta = do_pregrid(coord, sizeXYZ)
    nlayer = ((n_L[0]-1)*n_G[0]+1)*(n_L[2]-1) #each row of elements(X direction)
    nplane = ((n_L[0]-1)*n_G[0]+1)*((n_L[2]-1)*n_G[2]+1) #each plane XZ of elements
    # Initialize mesh: matrix with the node coordinates by rows.
    # Coordinates of the grid are ordered WE & SN
    mesh = np.zeros([nplane*((n_L[1]-1)*n_G[1]+1),3])
    # Initialize index: boolean vector to determine the real node rows in the mesh matrix.
    index = np.zeros(nplane*((n_L[1]-1)*n_G[1]+1))
    # Initialize num: Transformation matrix from global to local nodes. 
    # Each row corresponds to one element.
    # Elements are ordered from WE & SN
    num = np.array(np.zeros([n_G[0]*n_G[1]*n_G[2],N]))
    # Initialize Boundary: List of the list of boundary nodes of the basic element. 
    Boundary = [[],[],[],[],[],[]] # [min(x), max(x), min(y), max(y), min(z), max(z)]
    p = -1
    if sym !=-1:
        msize = [1,1,1]
        for rep in sym:
            msize[rep] = 2
        msize.append(N)
        ind_elem_GS = np.zeros(msize)
        msize.append(3)
        iiS = np.zeros(msize)
        coordS = np.zeros(msize)
    for i in range(msize[0]):
        for j in range(msize[1]):
            for k in range(msize[2]):
                indbol = np.array([bool(i),bool(j),bool(k)])
                coords = np.array(coord)
                coords[:,indbol] = np.reshape(np.repeat([delmax[indbol]],len(coord),axis=0),(N,sum(indbol))) - coord[:,indbol]
                coordS[i,j,k] = coords
                iis = np.array(ii)
                iis = [abs(max(ii[0])*indbol[0]-ii[0]),abs(max(ii[1])*indbol[1]-ii[1]),abs(max(ii[2])*indbol[2]-ii[2])]
                ind_elem_Gs = (iis[0] + iis[1]*nplane + iis[2]*nlayer/(n_L[2]-1))
                iiS[i,j,k] = np.reshape(iis, np.shape(iis)[:-1]).T
                ind_elem_GS[i,j,k] = np.reshape(ind_elem_Gs,(len(ind_elem_Gs)))
    del(coords,ind_elem_Gs, iis)
    ind_elem_GS = ind_elem_GS.astype(int)
    p = -1
    for j in range(n_G[1]):
        if 1 in sym:
            jT = j % 2
        else:
            jT = 0
        for k in range(n_G[2]):
#            
            if 2 in sym:
                kT = k % 2
            else:
                kT = 0
            for i in range(n_G[0]):
                if 0 in sym:
                    iT = i % 2
                else:
                    iT = 0
                p = p + 1
                num_r = np.reshape(ind_elem_GS[iT,jT,kT] + i*(n_L[0]-1) + k*nlayer+ j*nplane,(1,N))
                num[p,:] =  num_r
                index[num_r] = np.ones([len(num_r),1])
                [dx,dy,dz] = [float(delmax[0]*i), float(delmax[1]*j), float(delmax[2]*k)]
                mesh[num_r,:] = coordS[iT,jT,kT] + np.repeat([[dx,dy,dz]], N, axis=0) + do_shapeME(
                        coordS[iT,jT,kT]+[dx,dy,dz], n_L,n_G,delmax,shape,iiS[iT,jT,kT], i,j,k)
                Boundary, num[p,:], index = check_boundary(num, mesh, p, index, iiS[iT,jT,kT],i,j,k,n_L,n_G,Boundary,shape,DIM, delmax)
                
    index = index.astype("bool")
    num = num.astype("int")
    mesh = mesh[index] #Take only the nodes of interest from the rectangular grid
    # Hash from index imaginary rectangular grid to node grid.
    # The row number is the global node number and the value is the 
    # correspondant node to the imaginary rectangular grid.
    index = np.arange(0,len(index),1)[index]
    k = -1
    ind_hash = {}
    for i in index:
        k = k + 1
        ind_hash[i] = k
    Boundary = [do_translate(list,  ind_hash) for list in Boundary]
    
    del(indbol)
    return mesh, ind_hash, num, delmax, Boundary

def gen_edges(num, ind_hash, edges):
#    print('gen_edges(num, ind_hash, edges F) // edges_i G')
    # Initializing global adjecency list
    edges_i =  [[a] for a in range(len(ind_hash))]
    #consider only the lower triangle of the symetric matrix
    Nedges = 0
    for row in np.array(num):
        i = -1
        for node_in in row:
            i = i + 1
            node_con = row[edges[i]]
            node_in = ind_hash[node_in]
            for node_out in node_con[1:]:
                node_out = ind_hash[node_out]
                if (node_out not in edges_i[node_in])and(node_in not in edges_i[node_out]):
                    edges_i[node_in].append(node_out)
                    Nedges = Nedges + 1
    return edges_i, Nedges

def gen_CAD(filename, dim, mesh, lim, edges_i, faces_i, t, Export, View, foldername = ''):
#    print('gen_CAD(filename, dim, mesh, lim, edges_i, faces_i, t, Export, View, foldername = '' F)')
    if (Export or View):
        low, high = iter(lim)
#    dim=[0,1,2]
        parent_folder = os.getcwd()
        try:
            openFolder(foldername)
        except:
            pass
            # Initializing file
        if View:
            lines = list()
            try:
                os.remove(filename+'.pdf' )
                print("drawing '%s.pdf' replaced.\n" % filename)
            except: 
                print("drawing '%s.pdf' created.\n" % filename)
            if Export:
                try:
                    os.remove(filename+'.dxf' )
                    print("drawing '%s.dxf' replaced.\n" % filename)
                    drawing = dxf.drawing( filename+'.dxf' )
                except: 
                    print("drawing '%s.dxf' created.\n" % filename)
                    drawing = dxf.drawing( filename+'.dxf' )
                for n0 in edges_i:
                    for n1 in n0[1:]:
                        drawing.add(dxf.polyline([mesh[n0[0],dim], mesh[n1,dim]], layer = 'edges', startwidth =t, endwidth= t))
                drawing.save()
                if np.shape(faces_i)[1] != 0:
                    drawing = dxf.drawing( filename+'f'+'.dxf' )
                    for face in faces_i:
                        f=[]
                        for i in range(len(face)):
                            f.append(tuple(mesh[face[i],dim]))
                        f = dxf.face3d(f, flags=1)
                        f['layer'] = 'faces'
                        f['color'] = 7
                        drawing.add(f)
                        del(f)
                drawing.save()
            for n0 in edges_i:
                    for n1 in n0[1:]:
                        lines.append((mesh[n0[0],dim], mesh[n1,dim]))
            if len(dim) != 3:
                lc = mc.LineCollection(lines, linewidths=2, color='#CCCCCC')#(high[0]-low[0])/len(mesh), color='#CCCCCC')
                fig, aix = pl.subplots()
                for i in range(len(mesh)):
                    aix.annotate(str(i), xy=(mesh[i,dim]), family='Courier New',fontsize=16, color='red' )
                aix.set_xlim(low[0],high[0])
                aix.set_ylim(low[2],high[2])
                aix.add_collection(lc)
                aix.axis('equal')
                aix.axis('off')
                fig.show()
            else:
                fig = pl.figure()
                aix = fig.add_subplot(111, projection='3d', aspect='equal')
                aix.view_init(azim=120)
                lc= Line3DCollection(lines, linewidths=0.1, color='red')
                aix.add_collection3d(lc)
                aix.axis('equal')
                aix.set_xlim3d(low[0]-3,high[0]+3)
                aix.set_ylim3d(low[1]-3,high[1]+3)
                aix.set_zlim3d(low[2]-3,high[2]+3)
                # Hide grid lines
                aix.grid(False)
                # Hide axes ticks
    #            aix.set_xticks([])
    #            aix.set_yticks([])
    #            aix.set_zticks([])
                aix.autoscale_view()
            pl.savefig(filename + '.pdf')
        elif Export:
            try:
                os.remove(filename+'.dxf' )
                print("drawing '%s.dxf' replaced.\n" % filename)
                drawing = dxf.drawing( filename+'.dxf' )
            except: 
                print("drawing '%s.dxf' created.\n" % filename)
                drawing = dxf.drawing( filename+'.dxf' )
            for n0 in edges_i:
                for n1 in n0[1:]:
                    drawing.add(dxf.polyline([mesh[n0[0],dim], mesh[n1,dim]], layer = 'edges', startwidth =t, endwidth= t))
            drawing.save()
            if np.shape(faces_i)[1] != 0:
                drawing = dxf.drawing( filename+'f'+'.dxf' )
                for face in faces_i:
                    f=[]
                    for i in range(len(face)):
                        f.append(tuple(mesh[face[i],dim]))
                    f = dxf.face3d(f, flags=1)
                    f['layer'] = 'faces'
                    f['color'] = 7
                    drawing.add(f)
                    del(f)
            drawing.save()
        os.chdir(parent_folder)
    else:
        pass

def do_intersect(m1,n1,m2,n2):
#    print('do_intersect(m1, n1, m2, n2 F) // np.array([x, 0.0, x*m1+n1]) G')
    x = (n2-n1)/(m1-m2)
    return np.array([x, 0.0, x*m1+n1])

def gen_points_t(c0,c1,c2,T,sizeXYZ):   
#    print('gen_points_t(c0,c1,c2,T,sizeXYZ F) // newp, N G')
    N = 0
    newp =[]
    try:
        m01 = (c1[2]-c0[2])/(c1[0]-c0[0])
        m12 = (c2[2]-c1[2])/(c2[0]-c1[0])
        if m01 == 0 and m12 == 0:
            if bool(c1[2]):
                newp.append(c1-np.array([0,0,T/2]))
                N = N + 1
            if bool(c1[2]-sizeXYZ[2]):
                newp.append(c1+np.array([0,0,T/2]))
                N = N + 1
        else:
            t_v = np.array([[(c1[2]-c0[2]),0.0,-(c1[0]-c0[0])], 
                             [-(c1[2]-c0[2]),0.0,(c1[0]-c0[0])]])
            t_v = c1 + t_v*(T/(2.0*do_modeV(t_v[0]))) # point left and point right
            beta = (m01+m12)/2 # mean slope betweem
            if ((c1[0] == 0.0) or (c1[0] == sizeXYZ[0])) and ((c1[2] == 0.0) or (c1[2] == sizeXYZ[2])):
                for i in range(len(t_v)):
                    t = t_v[i]
                    N = N +1
                    if t[0] < 0.0 or t[0]> sizeXYZ[0]:
                        newp.append(do_intersect(m01, -t[0]*m01+t[2], m01-1, (c1[0]-t[0]*m01+t[2])))
                    if t[2] < 0.0 or t[2]> sizeXYZ[2]:   
                        newp.append(do_intersect(m01, -t[0]*m01+t[2], 0.0, c1[2]))
            elif (c1[0] == 0.0) or (c1[0] == sizeXYZ[0]):
                for i in range(len(t_v)):
                    t = t_v[i]
                    if t[0] < 0.0 or t[0] > sizeXYZ[0]:
                        N = N + 2
                        tt = np.array([[(c2[2]-c1[2]),0.0,-(c2[0]-c1[0])], 
                                         [-(c2[2]-c1[2]),0.0,(c2[0]-c1[0])]])[i]
                        tt = c1 + tt*(T/(2.0*do_modeV(tt)))
                        newp.append([do_intersect(m01, -t[0]*m01+t[2], m01-1, (c1[0]-t[0]*m01+t[2])),
                         do_intersect(m12,-tt[0]*m12+tt[2], m12-1, (c1[0]-tt[0]*m12+tt[2]))])
                    else:
                        N = N +1
                        newp.append(do_intersect(m01, -t[0]*m01+t[2], beta, -c1[0]*beta + c1[2]))
            elif (c1[2] == 0.0) or (c1[2] == sizeXYZ[2]):
                for i in range(len(t_v)):
                    t = t_v[i]
                    if t[2] < 0.0 or t[2] > sizeXYZ[2]:
                        N = N + 2
                        tt = np.array([[(c2[2]-c1[2]),0.0,-(c2[0]-c1[0])], 
                                         [-(c2[2]-c1[2]),0.0,(c2[0]-c1[0])]])[i]
                        tt = c1 + tt*(T/(2.0*do_modeV(tt)))
                        newp.append([do_intersect(m01, -t[0]*m01+t[2], 0.0, c1[2]), 
                                     do_intersect(m12,-tt[0]*m12+tt[2], 0.0, c1[2])])
                    else:
                        N = N + 1
                        newp.append(do_intersect(m01, -t[0]*m01+t[2], beta, -c1[0]*beta + c1[2]))
            else:
                for i in range(len(t_v)):
                    N = N + 1
                    t = t_v[i]
                    newp.append(do_intersect(m01, -t[0]*m01+t[2], beta, -c1[0]*beta + c1[2]))
    except: 
        if bool(c1[0]):
            newp.append(c1-np.array([T/2,0,0]))
            N = N +1
        if bool(c1[0]-sizeXYZ[0]):
            N = N + 1
            newp.append(c1+np.array([T/2,0,0]))
    return newp, N

#def gen_thickness_data(N, coord, t, sizeXYZ, edges):
##    print('gen_thickness_data(N, coord, t, sizeXYZ, edges F) // coco, edges, N + n G')
#    coco = []
#    c, n = gen_points_t(2*coord[0]-coord[1],coord[0],coord[1], t, sizeXYZ)
#    coco.extend(c)
#    edges.extend([[a] for a in range(N, N + n)])
#    edges[0].extend(list(range(N, N + n)))
#    N = N + n
#    for i in range(1,len(coord)-1):
#        c, n = gen_points_t(coord[i-1],coord[i],coord[i+1], t, sizeXYZ)
#        coco.extend(c)
#        edges.extend([[a] for a in range(N, N + n)])
#        edges[i].extend(list(range(N, N + n)))
#        N = n + N
#    c, n = gen_points_t(coord[i],coord[i+1], 2*coord[i+1]-coord[i], t, sizeXYZ)
#    coco.extend(c)
#    edges.extend([[a] for a in range(N, N + n)])
#    edges[i+1].extend(list(range(N, N + n)))
#    return coco, edges, N + n

        
def get_faces(N, t, sizeXYZ, coord, edges, faces=[], show = False):
#    print('get_faces(N, t, sizeXYZ, coord, edges, faces=[], centres=[], show = False F) // newcoord, faces, edges, N, bound_node G')
    sizeXYZ = np.reshape(np.asarray(sizeXYZ), (1,3))
    Nn = int(N)
    if t!=0:
        try:
            newface =[]
            for f in faces:
                coord, edges, newf, N = gen_thicknes_data(N, coord[f], t, sizeXYZ, edges)
                newface.extend(newf)
        except:
            coord, edges, newface, N = gen_thicknes_data(N, coord, t, sizeXYZ, edges)
#
#        # prepare the visualization
#        if show==True:
#            c = np.asarray(newcoord)
#            c = c[faces[l]]
#            polygon = Polygon(c[:,[0,2]], False, hatch ='/',linestyle=':')
#            patches.append(polygon)
#        # Visualization
#        if show == True:
#            fig, aix = plt.subplots()
#            aix.add_collection(PatchCollection(patches, alpha=0.4))
#            aix.axis('equal')
#            plt.xlim(-0,5)
#            plt.ylim(-0,5)
#            plt.show()
        bound_node = set(range(Nn, N ))
#        newcoord = np.asarray(newcoord)
#        N = len(newcoord) 
    else:
        bound_node = set()
        newface=[[]]
    return coord, edges, np.asarray(newface),  N, bound_node

def modeV(v):
#    print('modeV( v F // V**0.5 G')
    V = 0.0
    for vi in v:
        V = V + float(vi)**2
    return V**(0.5)

def proj_vect(v,w):
#    print('proj_vect( v, w F // proj G')
    # projects v over w:
    W = modeV(w)**2
    sc = 0.0
    for i in range(len(v)):
        sc = sc + float(v[i])*float(w[i])
    proj = np.zeros(len(v))
    for i in range(len(v)):
        proj[i] = float((sc/W))*w[i]
    return  proj

def intersect_lines(m1,n1,m2,n2):
#    print('intersect_lines( m1, n1, m2, n2 F // np.array([x, 0.0, x*m1+n1]) G')
    x = (n2-n1)/(m1-m2)
    return np.array([x, 0.0, x*m1+n1])

def check_out_bound( point, sizeXYZ, equal = True):
#    print('check_out_bound( point, sizeXYZ, equal = True F // B G')
    B = []
    if equal:
        for i in [0,2]:
            B.append(point[i] <= 0.0)
            B.append(point[i] >= sizeXYZ[0,i])
    else:
        for i in [0,2]:
            B.append(point[i] < 0.0)
            B.append(point[i] > sizeXYZ[0,i])
    return B

def create_points_t(c0,c1,c2, T, sizeXYZ, N):
#    print('  create_points_t(c0,c1,c2, T, sizeXYZ F // newp, index, N G')
    Nold = N
    index = []
    newp = []
    np.seterr(divide='ignore')
    t_v = np.array([[(c1[2]-c0[2]),0.0,-(c1[0]-c0[0])], 
                     [-(c1[2]-c0[2]),0.0,(c1[0]-c0[0])]])
    t_v = c1 + t_v*(T/(2.0*modeV(t_v[0])))  # point left and point right l01
    
    t_vt = np.array([[(c2[2]-c1[2]),0.0,-(c2[0]-c1[0])], 
                     [-(c2[2]-c1[2]),0.0,(c2[0]-c1[0])]])
    t_vt = c1 + t_vt*(T/(2.0*modeV(t_vt[0]))) #  point left and point right l12   
    deltaX = [(c1[0]-c0[0]) == 0.0 , (c2[0]-c1[0]) == 0.0]
    index = np.ones((2,2),dtype=int)*-1
    if all(deltaX):
        point = np.array([c1 + np.sign((c1[2]-c0[2]))*np.array([T/2,0.0,0.0]), c1+ np.sign(-(c1[2]-c0[2]))*np.array([T/2,0.0,0.0])])
        i = -1
        for p in point:
            i = i+1
            if not(any(check_out_bound(p, sizeXYZ, False))):
                index[:,i]=[N,N]
                newp.append(p)
                N = N + 1
    elif any(deltaX):
        for i in range(len(t_v)):
            if deltaX[0]:
                t = [t_v[i], t_vt[i]]
                mm = np.array([(c2[2]-c1[2])/(c2[0]-c1[0])])
                point = np.array(intersect_lines(mm[0], -t[1][0]*mm[0]+t[1][2], mm[0]-1, (t[0][0]-t[1][0]*mm[0]+t[1][2])))
                if not(any(check_out_bound(point, sizeXYZ, False))):
                    index[:,i]=[N,N]
                    newp.append(point)
                    N = N + 1
                else:
                    point = []
                    point.append(intersect_lines(mm[0], -t[1][0]*mm[0]+t[1][2], mm[0]-1, (c1[0]-t[1][0]*mm[0]+t[1][2])))
                    point.append(intersect_lines(mm[0], -t[1][0]*mm[0]+t[1][2], 0.0, c1[2]))
                    for p in point:
                        if not(any(check_out_bound(p, sizeXYZ, equal= False))):
                            index[1,i] = N
                            newp.append(p)
                            N = N+1
            elif deltaX[1]:
                t = [t_v[i], t_vt[i]]
                mm = np.array([(c1[2]-c0[2])/(c2[0]-c0[0])])
                point = np.array(intersect_lines(mm[0], -t[0][0]*mm[0]+t[0][2], mm[0]-1, (t[1][0]-t[0][0]*mm[0]+t[0][2])))
                if not(any(check_out_bound(point, sizeXYZ, False))):
                    index[:,i]=[N,N]
                    newp.append(point)
                    N = N + 1
                else:
                    point = [intersect_lines(mm[0], -t[0][0]*mm[0]+t[0][2], mm[0]-1, (c1[0]-t[0][0]*mm[0]+t[0][2])),
                            intersect_lines(mm[0], -t[0][0]*mm[0]+t[0][2], 0.0, c1[2])]
                    for p in point:
                        if not(any(check_out_bound(p, sizeXYZ, equal= False))):
                            index[0,i] = N
                            newp.append(p)
                            N = N+1
    elif (np.array([t_v[i] == t_vt[i] for i in range(len(t_v))]).all()):
        if ((c1[2]-c0[2]) == 0.0):
            point = np.array([c1 + np.sign((c1[0]-c0[0]))*np.array([0.0,0.0,-T/2]), c1+ np.sign(-(c1[0]-c0[0]))*np.array([0.0,0.0,-T/2])])
            i = -1
            for p in point:
                i = i+1
                if not(any(check_out_bound(p, sizeXYZ, False))):
                    index[:,i]=[N,N]
                    newp.append(p)
                    N = N + 1
        else:
            mm = (c1[2]-c0[2])/(c1[0]-c0[0])
            for i in range(len(t_v)):
                point = t_v[i]
                t = t_v[i]
                if not(any(check_out_bound(point, sizeXYZ, False))):
                    index[:,i]=[N,N]
                    newp.append(point)
                    N = N + 1
                else:
                    point = []
                    try:
                        point.append(intersect_lines(mm, -t[0]*mm+t[2], mm-1, (c1[0]-t[0]*mm+t[2])))
                    except:
                        pass
                    try:
                        point.append(intersect_lines(mm, -t[0]*mm+t[2], 0.0, c1[2]))
                    except:
                        pass
                    for p in point:
                        if not(any(check_out_bound(p, sizeXYZ, equal= False))):
                            index[:,i]= [N,N]
                            newp.append(p)
                            N = N + 1
                
    else:
        mm = np.array([(c1[2]-c0[2])/(c1[0]-c0[0]), (c2[2]-c1[2])/(c2[0]-c1[0])])
        beta = (np.sum(mm))/2
        for i in range(len(t_v)):
            t = [t_v[i], t_vt[i]]
            point = np.array(intersect_lines(mm[0], -t[0][0]*mm[0]+t[0][2], mm[1], -t[1][0]*mm[1]+t[1][2]))
            if not(any(check_out_bound(point, sizeXYZ, False))):
                index[:,i]=[N,N]
                newp.append(point)
                N = N + 1
            elif mm[0] == mm[1]:
                point = []
                try:
                    point.append(intersect_lines(mm[0], -t[0][0]*mm[0]+t[0][2], mm[0]-1, (c1[0]-t[0][0]*mm[0]+t[0][2])))
                except:
                    pass
                try:
                    point.append(intersect_lines(mm[0], -t[0][0]*mm[0]+t[0][2], 0.0, c1[2]))
                except:
                    pass
                for p in point:
                    if not(any(check_out_bound(p, sizeXYZ, equal= False))):
                        index[:,i]= [N,N]
                        newp.append(p)
                        N = N + 1
            else:
                for j in [0,1]:
                    point = []
                    try:
                        point.append(intersect_lines(mm[j], -t[j][0]*mm[j]+t[j][2], mm[j]-1, (c1[0]-t[j][0]*mm[j]+t[j][2])))
                    except:
                        pass
                    try:
                        point.append(intersect_lines(mm[j], -t[j][0]*mm[j]+t[j][2], 0.0, c1[2]))
                    except:
                        pass
                    for p in point:
                        if not(any(check_out_bound(p, sizeXYZ, equal= False))):
                            if mm[0] == mm[1]:
                                index[:,i]= [N,N]
                                j = 2
                            else:
                                index[j,i]= N
                            newp.append(p)
                            N = N + 1    
    return newp, index, N, Nold

def gen_t_faces(index, edges, faces, N, Nold, i = 1):
#    print('gen_t_faces(index, edges, faces, N, n = 0, i = 1 F // index, edges, faces, N + n G)')
    edges.extend([[a] for a in range(Nold, N)])
    face =  [[i-1],[i-1]]
#    facep = [[i-1],[i-1]]
    for irow in [index[1],index[2]]:
        for ii in [0,1]:
            if irow[ii] != -1:
                face[ii].append(irow[ii])
    for f in face:
        ff = f[:]
        if len(f) > 2:
            if f[1] not in edges[f[0]]:
                edges[f[0]].append(f[1])
            edges[f[1]].append(f[2])
            if f[2] not in edges[i]:
                edges[i].append(f[2])
            ff.append(i)
            faces.append(ff)        
    index = index[-2:]
    return index, edges, faces

def gen_thicknes_data(N, coord, t, sizeXYZ, edges):
#    print('gen_t_coord(N, coord, t, sizeXYZ, edges) F // coco, edges, faces, N G')
    coco = []
    ind = []
    faces=[]
    if (0 in edges[np.shape(coord)[0]-1]) or (int(np.shape(coord)[0]-1) in edges[0]):
         c, index, N, Nold = create_points_t(coord[-1],coord[0],coord[1], t, sizeXYZ, N)
         lasti = index[:]
    else:
         c, index, N, Nold = create_points_t(2*coord[0]-coord[1],coord[0],coord[1], t, sizeXYZ, N)
    coco.extend(c)
    ind.extend(index)
    edges.extend([[a] for a in range(Nold, N)])
    for i in range(1,len(coord)-1):
        c, index, N, Nold = create_points_t(coord[i-1],coord[i],coord[i+1], t, sizeXYZ, N)
        coco.extend(c)
        ind.extend(index)
        ind, edges, faces = gen_t_faces(ind, edges, faces, N, Nold, i)
    if (0 in edges[np.shape(coord)[0]-1]) or (int(np.shape(coord)[0]-1) in edges[0]):
        c, index, N, Nold = create_points_t(coord[i],coord[i+1], coord[0], t, sizeXYZ, N)      
        coco.extend(c)
        ind.extend(index)
        ind, edges, faces = gen_t_faces(ind, edges, faces, N, Nold, i+1)
        ind.extend(lasti)
        ind, edges, faces = gen_t_faces(ind, edges, faces, N, Nold, i+2)
    else:
        c, index, N, Nold = create_points_t(coord[i],coord[i+1], 2*coord[i+1]-coord[i], t, sizeXYZ, N)
        coco.extend(c)
        ind.extend(index)
        ind, edges, faces = gen_t_faces(ind, edges, faces, N, Nold, i+1)
    coord = np.concatenate((coord, np.reshape(np.asarray(coco),(len(coco),3))))
    return coord, edges, faces, N 
   

def gen_faces(num, ind_hash, faces):
#    print('gen_faces(num, ind_hash, faces F) // faces_i G')
    if np.size(faces) != 0: 
        # Initializing global face list
        faces_i =  []
        for num_r in num:
            for face in faces:
                nf = num_r[:,face].tolist()[0]
                nnf = []
                for n in nf:
                    nnf.append(ind_hash[n])
                faces_i.append(nnf)
    else:
        faces_i = faces
    return faces_i

def gen_mesh( N, nodes, edges, faces, bound_node, sizeXYZ, meshSize=[], DoMesh = False):
#    print('gen_mesh( N, nodes, edges, faces, bound_node, sizeXYZ, meshSize=[], DoMesh = False F // nodes, edges, new_faces, N, bound_node G')
    if not(DoMesh):
        new_faces = faces
    else:
        # Generate edge mesh
        minsize, maxsize = [float(m) for m in meshSize]
        old_edges = copy.deepcopy(edges)
        newc = []
        k =-1
        del(edges)
        edges = [[a] for a in range(N)]
        for edge in old_edges:
            for no in edge[1:]:
                nedges = []
                ni=edge[0]
                vect = nodes[no]-nodes[ni]
                dist = np.linalg.norm(vect)
                div = 1
                while ((dist>maxsize) and ((dist/(div+1)>=minsize) and (dist/(div)>maxsize))):
                    div = div + 1
                newcoord = nodes[ni] 
                # Add new edges and new nodes to the data
                for i in range(1,div):
                    k = k + 1
                    newcoord = newcoord + vect*(1./div)
                    newc.append(newcoord)
                    nedges.append(k+N)
                    edges.append([k+N])
                if len(bound_node) !=0:
                    if (ni in bound_node) and (no in bound_node):
                        # Add new bound nodes to the bound node set
                        bound_node = bound_node.union(set(nedges))
                if len(nedges) == 0:
                    edges[ni].append(no)
#                    try:
#                        bound_node.remove(ni)
#                    except:
#                        try:
#                            bound_node.remove(no)
#                        except:
##                            pass
                if len(nedges) != 0:
                    edges[ni].append(nedges[0])
                    edges[no].append(nedges[-1])
                if len(nedges)>=2:
                    for nn in nedges[0:-1]:
                        edges[nn].append(nn+1)  
        # Add new nodes to the node array
        nodes = np.concatenate((nodes, np.reshape(newc,(len(newc),3))))
        N = len(nodes) # new number of nodes
        new_faces = faces

        # TRIANGULATE faces if needed
        if (np.sum(faces) != 0) and (sizeXYZ!=0):

            triang = Delaunay(nodes[:,[0,2]])
            triangg = do_trimming(triang.simplices, bound_node)
#            delaunay_plot_2d(triang)
            plt.figure()
            plt.axis('equal')
            plt.axis('off')
            plt.triplot(nodes[:,0], nodes[:,2], triangg)
            plt.show()
            plt.savefig('tria' + '.pdf')
            new_faces = triangg
    #    if np.size(faces) != 0:
    #        new_nodes, new_edges, new_faces = gen_mesh2D(faces, delta, ii, n_L)
    #        for d in [0,1,2]:
    #            new_nodes[:,d], delt, idim, n_ = check_Dim_Mult(new_nodes[:,d], sizeXYZ[d])
    #            delta.append(delt)
    #            delmax.append(float(max(delt)))
    #            ii.append(idim)
    #            n_L.append(n_)
    #        delmax = np.array(delmax)
    #        del(delt, n_, idim)
    return nodes, edges, new_faces, N, bound_node


def do_trimming(grid, bound_node):
#    print('do_trimming(grid, bound_node F) // gridd G')
    gridd=[]
    for facet in grid:
        if not(all(f in bound_node for f in facet)):
            gridd.append(facet)
        else:
            pass
    gridd = np.asanyarray(gridd)
    return gridd

def gen_mesh2D(faces, delta, ii, n_L):
#    print('gen_mesh2D(faces, delta, ii, n_L F) // tri_faces G')
    # Initialize new_faces:
    tri_faces = []
    x,z=np.meshgrid(delta[0],delta[2])
    x = x.flatten()
    z = z.flatten()
    #define 2D points, asx
    #input data for the Delaunay triangulation of U
    tri = Delaunay(np.vstack([x,z]).T)#triangulate the rectangle U
    
    for face in faces:
        index_L = ii[0][face] + ii[2][face]*n_L[0]
        for row in tri.simplices:
            if any(row) not in index_L:    
                pass
            else:
                tri_faces.append(tri.simplices)
    return tri_faces

def do_pregrid(nodes, sizeXYZ):
#    print('do_pregrid( nodes, sizeXYZ F) // n_L, ii, delmax, delta G')
    # Generating data to prepare the mesh
    ii = []
    n_L = []
    delta = []
    delmax = []
    for d in [0,1,2]:
        nodes[:,d], delt, idim, n_ = check_Dim_Mult(nodes[:,d], sizeXYZ[d])
        delta.append(delt)
        delmax.append(float(max(delt)))
        ii.append(idim)
        n_L.append(n_)
    delmax = np.array(delmax)
    del(delt, n_, idim)
    return n_L, ii, delmax, delta



def do_shapeME(coord, n_L, n_G, delmax ,shape, ii, i=0,j=0, k=0, coef_rand =0.625):
#    print('do_shapeME(coord, n_L, n_G, delmax ,shape, ii, i=0, j=0, k=0, coef_rand =0.625 F // coord G')
    if shape == '':
        coord = np.zeros([1,3])
        
    elif shape == 'random':
        for r in range(len(coord)):
            coord[r,:] = random.uniform(-delmax/coef_rand,delmax/coef_rand)
            
    elif shape == 'progressive':
        for r in range(len(coord)):            
            x = float((i+ii[r,0])**2)
            z = float((k+ii[r,2])**2)
            try:
                y = float((j+ii[r,1])**2)
            except:
                y = 0.0
            coord[r,:] = np.array([x,y,z])
            
    elif shape == 'flower':
        TOT_angleX = 5*math.pi
        TOT_angleY = 1.0
        a = 5.0
        b = 0.3
        r = -1
        for x,y,z in coord:
            r = r+1
            fract_angleX = float((float(i*(n_L[0]-1)+ii[r,0]))/(n_G[0]*(n_L[0]-1)))
            fract_angleY = float((float(j*(n_L[1]-1)+ii[r,1])+1)/(n_G[1]*(n_L[1]-1)))
            xt = (a+z)*math.cos(TOT_angleX*fract_angleX)*math.sin(TOT_angleY*fract_angleY)-x
            yt = (a+z)*math.sin(TOT_angleX*fract_angleX)*math.sin(TOT_angleY*fract_angleY)-y
            zt = (a+z)*(math.cos(TOT_angleY*fract_angleY)+math.log(math.tan(0.5*TOT_angleY*fract_angleY)))+b*TOT_angleX*fract_angleX-z
            coord[r,:] = np.array([xt,yt,zt])
        
            
    elif shape == 'ring':
        TOT_angleX = 2.0*math.pi
        TOT_angleY = 2.0*math.pi
        R = 25.0
        Ri = delmax[2]
        r = -1
        for x,y,z in coord:
            r = r+1
            fract_angleX = float((float(i*(n_L[0]-1)+ii[r,0]))/(n_G[0]*(n_L[0]-1)))
            fract_angleY = float((float(j*(n_L[1]-1)+ii[r,1]))/(n_G[1]*(n_L[1]-1)))
            xt = math.cos(TOT_angleX*fract_angleX)*(R+(Ri+z)*math.cos(TOT_angleY*fract_angleY))-x
            yt = math.sin(TOT_angleX*fract_angleX)*(R+(Ri+z)*math.cos(TOT_angleY*fract_angleY))-y
            zt = (Ri+z)*math.sin(TOT_angleY*fract_angleY)-z
            coord[r,:] = np.array([xt,yt,zt])
            
    elif shape == 'galaxy':
        aa, a = 3, 1
        TOT_angleX = 2*math.pi
        TOT_angleY = 6*math.pi
        r = -1
        for x,y,z in coord:
            r = r+1
            fract_angleX = float((float(i*(n_L[0]-1)+ii[r,0]))/(n_G[0]*(n_L[0]-1)))
            fract_angleY = float((float(j*(n_L[1]-1)+ii[r,1]))/(n_G[1]*(n_L[1]-1)))
            xt = (a+z)*(aa+math.cos(TOT_angleY*fract_angleY/2)*math.sin(TOT_angleX*
                    fract_angleX)-math.sin(TOT_angleY*fract_angleY/2)*math.sin(2*
                                         TOT_angleX*fract_angleX))*math.cos(TOT_angleY*fract_angleY) - x
            yt = (a+z)*(aa+math.cos(TOT_angleY*fract_angleY/2)*math.sin(TOT_angleX*
                    fract_angleX)-math.sin(TOT_angleY*fract_angleY/2)*math.sin(2*
                                         TOT_angleX*fract_angleX))*math.sin(TOT_angleY*fract_angleY) -y
            zt = (a+z)*(math.sin(TOT_angleY*fract_angleY/2)*math.sin(TOT_angleX*
                         fract_angleX)+math.cos(TOT_angleY*fract_angleY/2)*
                         math.sin(2*TOT_angleX*fract_angleX)) -  z
            coord[r,:] = np.array([xt,yt,zt])
            
    elif shape == 'triplefris':
        TOT_angleX = math.pi
        TOT_angleY = math.pi
        r = -1
        for x,y,z in coord:
            r = r+1
            fract_angleX = float((float(i*(n_L[0]-1)+ii[r,0]))/(n_G[0]*(n_L[0]-1)))
            fract_angleY = float((float(j*(n_L[1]-1)+ii[r,1]))/(n_G[1]*(n_L[1]-1)))
            xt = (1+z)*((math.cos(TOT_angleX*fract_angleX)*((1/3)*math.sqrt(2.0)*
                          math.cos(TOT_angleX*fract_angleX)*math.cos(2*TOT_angleY*fract_angleY) + (2/3)*
                          math.sin(TOT_angleX*fract_angleX)*math.cos(TOT_angleY*fract_angleY)))/(1 - math.sqrt(2.0)*
                          math.sin(TOT_angleX*fract_angleX)*math.cos(TOT_angleX*fract_angleX)*
                          math.sin(3*TOT_angleY*fract_angleY))) - x
            yt = (1+z)*((math.cos(TOT_angleX*fract_angleX)*((1/3)*math.sqrt(2.0)*math.cos(TOT_angleX*
                          fract_angleX)*math.sin(2*TOT_angleY*fract_angleY) -(2/3)*
                          math.sin(TOT_angleX*fract_angleX)*math.cos(TOT_angleY*fract_angleY)))/ (1 - math.sqrt(2.0)*
                          math.sin(TOT_angleX*fract_angleX)*math.cos(TOT_angleX*fract_angleX)*
                          math.sin(3*TOT_angleY*fract_angleY))) - y
            zt = (1+z)*((math.cos(TOT_angleX*fract_angleX)*math.cos(TOT_angleX*fract_angleX))/(1 - math.sqrt(2.0)*
                  math.sin(TOT_angleX*fract_angleX)*math.cos(TOT_angleX*fract_angleX)*math.sin(3*
                          TOT_angleY*fract_angleY)) - 1) - z
            coord[r,:] = np.array([xt,yt,zt])
            
    elif shape == 'quadratic':
        r = -1
        for x,y,z in coord:
            r = r +1
            xt = -0.5*x
            yt = 0
            zt = (0.005*(x-5)**2)
            coord[r,:] = np.array([xt,yt,zt])
            
    elif shape == 'shell':
        r = -1
        for x,y,z in coord:
            r = r +1
            xt = 0
            yt = 0
            zt = (0.05*(x-delmax[0]*(n_G[0]*(n_L[0]-1))/2.0)**2) + (0.05*(y-delmax[1]*(n_G[1]*(n_L[1]-1))/2.0)**2)  
            coord[r,:] = np.array([xt,yt,zt])
            
    elif (shape == 'circle') or (shape == 'cylinder'):
        r = -1
        TOT_angle = 2*math.pi
        innerR = 3 # inner radious
        for x,y,z in coord:
            r = r+1
            fract_angle = float((float(i*(n_L[0]-1)+ii[r,0]))/(n_G[0]*(n_L[0]-1)))
            xt = (z+innerR)*math.cos(fract_angle*TOT_angle)-x
            yt = 0
            zt = (z+innerR)*math.sin(fract_angle*TOT_angle)-z
            coord[r,:] = np.array([xt,yt,zt])
                
    elif shape == 'catenoid':
        TOT_angleX = math.pi
        TOT_angleY = 10.0
        c = 5.0 #inner radious
        r = -1
        for x,y,z in coord:
            r = r+1
            fract_angleX = float((float(i*(n_L[0]-1)+ii[r,0]))/(n_G[0]*(n_L[0]-1)))
            fract_angleY = float((float(j*(n_L[1]-1)+ii[r,1]))/(n_G[1]*(n_L[1]-1)))-0.5
            xt = (z+c)*math.cos(TOT_angleX*fract_angleX-math.pi)*math.cosh(TOT_angleY*fract_angleY/c)-x
            yt = (z+c)*math.sin(TOT_angleX*fract_angleX-math.pi)*math.cosh(TOT_angleY*fract_angleY/c)-y
            zt = (TOT_angleY*fract_angleY)-z
            coord[r,:] = np.array([xt,yt,zt])
            
    elif shape == 'sphere':
        TOT_angleX = 2*math.pi
        TOT_angleY = math.pi
        innerR = 5 #inner radious
        r = -1
        for x,y,z in coord:
            r = r+1
            fract_angleX = float((float(i*(n_L[0]-1)+ii[r,0]))/(n_G[0]*(n_L[0]-1)))
            fract_angleY = float((float(j*(n_L[1]-1)+ii[r,1]))/(n_G[1]*(n_L[1]-1)))
            xt = (z+innerR)*math.cos(TOT_angleX*fract_angleX)*math.sin(TOT_angleY*fract_angleY)-x
            yt = (z+innerR)*math.cos(TOT_angleY*fract_angleY)-y
            zt = (z+innerR)*math.sin(TOT_angleX*fract_angleX)*math.sin(TOT_angleY*fract_angleY)-z
            coord[r,:] = np.array([xt,yt,zt])
            
    elif shape == 'NURBS':
        r = -1
        for x,y,z in coord:
            r = r+1
            xt= 0
            yt= 0
            zt = (0.005*(x-delmax[0]*(n_G[0]*(n_L[0]-1))/2.0)*(y-delmax[1]*(n_G[1]*(n_L[1]-1))/2.0))**3
            +(0.005*(x-delmax[0]*(n_G[0]*(n_L[0]-1))/2.0)*(y-delmax[1]*(n_G[1]*(n_L[1]-1))/2.0))**2
            +(0.005*(x-delmax[0]*(n_G[0]*(n_L[0]-1))/2.0)*(y-delmax[1]*(n_G[1]*(n_L[1]-1))/2.0))
            coord[r,:] = np.array([xt,yt,zt])
            
    return coord

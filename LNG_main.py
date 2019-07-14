#!/usr/bin/env python
# -----------------------------------------------------------------------------
# LNG_main.py
# -----------------------------------------------------------------------------

import numpy as np
import os
import math
from dxfwrite import DXFEngine as dxf
from matplotlib import collections as mc
import pylab as pl
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import LNG_engine 
       
class RVE:
    def __init__(self, eletypeID = None, sizeXYZ=[1.0,0.0,1.0], t = 0.0):
        self.sizeXYZ = np.reshape(np.asarray(sizeXYZ), (1,3))
        self.bound_node = {}
        if eletypeID != None:
            
            self._eletype = self._eletype(eletypeID)
        else: #Define your element
            self.eletypeID = eletypeID
            self.nodes = None
            self.edges = None
            self.faces = None
            self.sym = None
            self.dim = None
        self._get_input()
        self.t = t
        if t !=0.0:
            self._get_faces()
    def N(self):
        return len(self.nodes)
    def _eletype(self, eletypeID):
        self.eletypeID = eletypeID
        self.sym = []
        self.dim = [0,2]
        self.faces = np.array([[]])
        def e0(self):
            self.nodes = np.array([[0,0,0],[1,0,0]], dtype="float")
            self.edges = np.array([[0,1]])
            self.dim = [0]
        def e1(self):
            self.nodes = np.array([[1,0],[0,0],[1,1],[0,1]], dtype="float")
            self.edges = np.array([[0,1],[1,2],[2,3]])
            self.sym = [0,2]    
        def e3(self):
            self.nodes = np.array([[0,0],[0,0.5],[1,1]], dtype="float")
            self.edges = np.array([[0,1],[1,2]])
            self.sym = [0,2]
        def e4(self):
            self.nodes = np.array([[0,0],[0,1],[1,1],[1,0]], dtype="float")
            self.edges = np.array([[0,1],[1,2],[2,3],[3,0]])
            self.faces = np.array([[]])
        def e5(self): # cube
            self.nodes = np.array([[0,0,0],[1,0,0],[0,1,0],[1,1,0],[0,0,1],[1,0,1],[0,1,1],[1,1,1]], dtype="float")
            self.edges = np.array([[0,1],[0,2],[2,3],[1,3],[0,4],[1,5],[2,6],[3,7],[4,5],[4,6],[5,7],[6,7]])
            self.faces = np.array([[0,1,3,2],[4,5,7,6],[0,1,5,4],[0,2,6,4],[1,3,7,5],[2,3,7,6]])
            self.dim = (0,1,2)
        def e6(self): # X shape
            self.nodes = np.array([[0,0,0],[1,0,0],[0,1,0],[1,1,0],[0,0,1],[1,0,1],[0,1,1],[1,1,1]], dtype="float")
            self.edges = np.array([[0,7],[1,6],[2,5],[3,4]])
            self.faces = np.array([[]])
            self.dim = (0,1,2)
        def e7(self): # hexa-tetrahedrom
            self.nodes = np.array([[0,0,0],[1,0,0],[0,1,0],[1,1,0],[0,0,1],[1,0,1],[0,1,1],[1,1,1]], dtype="float")
            self.edges = np.array([[0,1],[0,2],[0,3],[2,3],[1,3],[2,7],[1,7],[1,4],[2,4],[4,6],[4,5],[5,7],[6,7],[4,7]])
            self.faces = np.array([[0,2,3],[0,1,3],[4,6,7],[4,5,7],[2,4,7],[1,4,7]])
            self.sym = [0,1,2]
            self.dim = (0,1,2)
        def e8(self): # octet-truss
            self.nodes = np.array([[0,0,0],[1,1,0],[1,0,1],[0,1,1]], dtype="float")
            self.edges = np.array([[0,1],[0,2],[0,3],[2,3],[1,2],[1,3]])
            self.faces = np.array([[]])
            self.sym = [0,1,2]
            self.dim = (0,1,2)
        #    print(' get_basic_elem( e F) // nodes, edges, faces, sym, dim G')
        try:
            eval('e'+str(self.eletypeID)+'(self)')  
        except:
            raise NameError('Element type %s is not defined.' %self.eletypeID)
       
    def _get_input(self):
        if type(self.sizeXYZ)!= np.array:
            self.sizeXYZ = np.asmatrix(self.sizeXYZ)
        # Amount of nodes
        self._scale = np.divide(self.sizeXYZ[:,self.dim],np.max(self.nodes, axis=0))
        self.nodes = np.multiply(self.nodes,self._scale)
        # Adding missing dimension in case of node coordinates given in 2D
        if np.shape(self.nodes)[1] != 3:
            self.nodes = np.array([self.nodes[:,0],np.zeros([np.shape(self.nodes)[0],1]),self.nodes[:,1]]).T[0]
            self._scale = np.array([self._scale[:,0],np.zeros([np.shape(self._scale)[0],1]),self._scale[:,1]]).T[0]
        T_edges = [[a] for a in range(self.N())]
        for i,j in self.edges:
            if j not in T_edges[i]:
                T_edges[i].append(j)
        self.edges = T_edges
    #    print('get_input(coord, edges, faces, sizeXYZ, dim F) // coord, T_edges, N, scale G')
        assert np.sum(self.sizeXYZ!=0)==len(self.dim), ValueError('Element cell size (RVE.sizeXYZ) and element dimensions (self.dim) do not match: sizeXYZ=%s dim=%s'% (self.sizeXYZ, self.dim))
    
    def _get_faces(self):
        Nn = self.N()
        LNG_engine.gen_thickness_data(self)
        self.bound_node = set(range(Nn, self.N() ))
        
    def gen_mesh(self, meshSize):
        if meshSize != []:
            LNG_engine.gen_mesh(self, meshSize)
            print("RVE remeshed with a meshSize of {}, if you want to visualize it type {}.showmesh() ".format(meshSize, type(self).__name__))
        else:
            print("Lattice RVE not remeshed.")
            

    def printDATA(self):
        for x in ['nodes', 'edges', 'faces']:
            try:
                print(str(x)+': {}\n'.format(eval('self.'+x)))
            except:
                pass
    def show(self):
        
        def show2F(self):
            pass
        def show2T(self):
            verts = [self.nodes[f][:,self.dim] for f in self.faces]
            print(verts)
            fig, ax = plt.subplots()
            # Make the collection and add it to the plot.
            coll = PolyCollection(verts, hatch ='/',linestyle=':')
            ax.add_collection(coll)
            ax.autoscale_view()
            plt.show()   
        def show3F(self):
            pass
        
        eval('show'+str(len(self.dim))+str(bool(self.t))[0]+'(self)')

    def showmesh(self):
        plt.figure()
        plt.axis('equal')
        plt.axis('off')
        plt.triplot(self.nodes[:,0], self.nodes[:,2], self.faces)
        plt.show()
                
class LatticeStructure:
    def __init__(self, _RVE, n_G=[1,1,1], shape=''):
        assert type(_RVE) == RVE, 'Please enter an RVE object'
        self.n_G = n_G
        self.shape = shape
        self.RVE = _RVE

    def gen_nodegrid(self):  
        LNG_engine.check_dimensions(self)
        LNG_engine.do_pregrid(self)
        N = self.RVE.N()
        nlayer = ((self.n_L[0]-1)*self.n_G[0]+1)*(self.n_L[2]-1) #each row of elements(X direction)
        nplane = ((self.n_L[0]-1)*self.n_G[0]+1)*((self.n_L[2]-1)*self.n_G[2]+1) #each plane XZ of elements
        # Initialize mesh: matrix with the node coordinates by rows.
        self.nodes = np.zeros([nplane*((self.n_L[1]-1)*self.n_G[1]+1),3]) # Coordinates of the grid are ordered WE & SN
        # Initialize index: boolean vector to determine the real node rows in the mesh matrix.
        index = np.zeros(nplane*((self.n_L[1]-1)*self.n_G[1]+1))
        # Initialize num: Transformation matrix from global to local nodes.
        # Each row corresponds to one element.
        # Elements are ordered from WE & SN
        self.num = np.array(np.zeros([self.n_G[0]*self.n_G[1]*self.n_G[2],N]))
        # Initialize boundary: List of the list of boundary nodes of the basic element.
        self.boundary = [set(),set(),set(),set(),set(),set()] # Boundaries[min(x), max(x), min(y), max(y), min(z), max(z)]
        self.bound, self.TOT_angleX, self.TOT_angleY= [[],[],[]], None, None
        p = -1
        if  self.RVE.sym !=-1:
            msize = [1,1,1]
            for rep in  self.RVE.sym:
                msize[rep] = 2
            msize.append(N)
            ind_elem_GS = np.zeros(msize)
            msize.append(3)
            self.iiS = np.zeros(msize)
            coordS = np.zeros(msize)
        for i in range(msize[0]):
            for j in range(msize[1]):
                for k in range(msize[2]):
                    indbol = np.array([bool(i),bool(j),bool(k)])
                    coords = np.array(self.RVE.nodes)
                    coords[:,indbol] = np.reshape(np.repeat([self.delmax[indbol]],N,axis=0),(N,sum(indbol))) - self.RVE.nodes[:,indbol]
                    coordS[i,j,k] = coords
                    iis = np.array(self.ii)
                    iis = [abs(max(self.ii[0])*indbol[0]-self.ii[0]),abs(max(self.ii[1])*indbol[1]-self.ii[1]),abs(max(self.ii[2])*indbol[2]-self.ii[2])]
                    ind_elem_Gs = (iis[0] + iis[1]*nplane + iis[2]*nlayer/(self.n_L[2]-1))
                    self.iiS[i,j,k] = np.reshape(iis, np.shape(iis)[:-1]).T
                    ind_elem_GS[i,j,k] = np.reshape(ind_elem_Gs,(len(ind_elem_Gs)))
        del(coords,ind_elem_Gs, iis, indbol)
        ind_elem_GS = ind_elem_GS.astype(int)
        p = -1
        self.bound[1]=[]
        for j in range(self.n_G[1]):
            if 1 in  self.RVE.sym:
                jT = j % 2
            else:
                jT = 0
            self.bound[2]=[]
            for k in range(self.n_G[2]):
                if 2 in self.RVE.sym:
                    kT = k % 2
                else:
                    kT = 0
                self.bound[0]=[]
                for i in range(self.n_G[0]):
                    if 0 in  self.RVE.sym:
                        iT = i % 2
                    else:
                        iT = 0
                    p = p + 1
                    self.I = [i,j,k]
                    self.num_r = np.reshape(ind_elem_GS[iT,jT,kT] + i*(self.n_L[0]-1) + k*nlayer+ j*nplane,(1,N))
                    [dx,dy,dz] = [np.float32(self.delmax[0]*i), np.float32(self.delmax[1]*j), np.float32(self.delmax[2]*k)]
                    LNG_engine.do_shapeME(self, coordS[iT,jT,kT]+[dx,dy,dz],self.iiS[iT,jT,kT])
                    self.num[p,:] = self.num_r
                    index[self.num_r] = np.ones([len(self.num_r),1])        
        index = index.astype("bool")
        self.num = self.num.astype("int")
        self.nodes = self.nodes[index].astype("float32") #Take only the nodes of interest from the rectangular grid
        # Hash from index imaginary rectangular grid to node grid.
        # The row number is the global node number and the value is the
        # correspondant node to the imaginary rectangular grid.
        index = np.arange(0,len(index),1)[index]
        self._ind_hash, k = {}, -1
        for i in index:
            k = k + 1
            self._ind_hash[i] = k
        self.boundary = [LNG_engine.do_translate(list(b),  self._ind_hash) for b in self.boundary]
        del (self.bound, self.I)
        
    def gen_edges(self):
        # Initializing global adjecency list
        self.edges =  [[a] for a in range(len(self._ind_hash))]
        #consider only the lower triangle of the symetric matrix
        self.num_edges = 0
        for row in self.num:
            i = -1
            for node_in in row:
                i = i + 1
                node_con = row[self.RVE.edges[i]]
                node_in = self._ind_hash[node_in]
                for node_out in node_con[1:]:
                    node_out = self._ind_hash[node_out]
                    if (node_out not in self.edges[node_in])and(node_in not in self.edges[node_out]):
                        self.edges[node_in].append(node_out)
                        self.num_edges = self.num_edges + 1
                        
    def gen_faces(self):
        
        # Initializing global face list
        try:
            self.faces =  []
            for num_r in self.num:
                for face in self.RVE.faces:
                    nf = num_r[face].tolist()
                    nnf = []
                    for n in nf:
                        nnf.append(self._ind_hash[n])
                    self.faces.append(nnf)
        except:
            pass

    def gen_CAD(self, filename ='', foldername = ''):
        if filename =='': filename = self.shape +'_RVE'+ str(self.RVE.eletypeID
                                                            )+'t'+str(self.RVE.t) +'nG'+ ''.join(map(str, self.n_G[:]))
        if foldername == '': foldername = self.shape
        parent_folder = os.getcwd()
        try:
            LNG_engine.openFolder(foldername)
            print('folder ')+ foldername + ' created inside /drawings_/.'
        except:
            NameError: 'path: '+parent_folder+foldername+'not found'
            # Initializing file
        try:
            os.remove(filename+'.dxf' )
            print("drawing '%s.dxf' replaced.\n" % filename)
            drawing = dxf.drawing( filename+'.dxf' )
        except:
            print("drawing '%s.dxf' created.\n" % filename)
            drawing = dxf.drawing( filename+'.dxf' )
        for n0 in self.edges:
            for n1 in n0[1:]:
                drawing.add(dxf.polyline([self.nodes[n0[0],self.RVE.dim], self.nodes[n1,self.RVE.dim]], layer = 'edges'))
        drawing.save()
        if np.shape(self.faces)[1] != 0:
            drawing = dxf.drawing( filename+'f'+'.dxf' )
            for face in self.faces:
                f=[]
                for i in range(len(face)):
                    f.append(tuple(self.nodes[face[i],self.RVE.dim]))
                f = dxf.face3d(f, flags=1)
                f['layer'] = 'faces'
                f['color'] = 7
                drawing.add(f)
                del(f)
        drawing.save()
        os.chdir(parent_folder)

    def show(self, lim = [], filename ='', foldername = ''):
        lines = list()
        if lim == []: lim = [np.zeros(3),self.delmax*self.n_G]
        lim, low, high = [], lim[0], lim[-1]
        if filename =='': filename = self.shape +'_RVE'+ str(self.RVE.eletypeID
                                                            )+'t'+str(self.RVE.t) +'nG'+ ''.join(map(str, self.n_G[:]))
        if foldername == '': foldername = self.shape

        parent_folder = os.getcwd()
        try:
            LNG_engine.openFolder(foldername)
            print('folder ')+ foldername + ' created inside /drawings_/.'

        except:
            NameError: 'path: '+parent_folder+foldername+'not found'
        try:
            os.remove(filename+'.pdf' )
            print("drawing '%s.pdf' replaced.\n" % filename)
        except:
            print("drawing '%s.pdf' created.\n" % filename)        
        for n0 in self.edges:
                for n1 in n0[1:]:
                    lines.append((self.nodes[n0[0],self.RVE.dim], self.nodes[n1,self.RVE.dim]))

        if len(self.RVE.dim) != 3:
            lc = mc.LineCollection(lines, color='grey')#(high[0]-low[0])/(len(self.nodes)), color='#CCCCCC')
            fig, aix = pl.subplots()
#            for i in range(len(self.nodes)):
#                aix.annotate(str(i), xy=(self.nodes[i,self.RVE.dim]), family='Courier New',fontsize=16, color='red' )
            aix.set_xlim(low[0],high[0])
            aix.set_ylim(low[2],high[2])
            aix.add_collection(lc)
            aix.axis('equal')
            aix.axis('off')
            fig.show()
        else:
            fig = pl.figure()
            aix = fig.add_subplot(111, projection='3d')
            aix.view_init(azim=120)
            lc= Line3DCollection(lines, linewidths=1, color='red')
            aix.add_collection3d(lc)
            aix.set_xlim3d(low[0]-3,high[0]+3)
            aix.set_ylim3d(low[1]-3,high[1]+3)
            aix.set_zlim3d(low[2]-3,high[2]+3)
            # Hide grid lines
            aix.grid(False)
            # Hide axes ticks
            aix.set_xticks([])
            aix.set_yticks([])
            aix.set_zticks([])
            aix.autoscale_view()
        pl.savefig(filename + '.pdf')
        os.chdir(parent_folder)


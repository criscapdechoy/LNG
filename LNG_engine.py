#!/usr/bin/env python
# -----------------------------------------------------------------------------
# LNG_engine.py
# -----------------------------------------------------------------------------

import numpy as np
import os
import math
from scipy.spatial import Delaunay
import copy
import matplotlib.pyplot as plt

def modeV(v):
#    print('modeV( v F // V**0.5 G')
    V = 0.0
    for vi in v:
        V = V + np.float32(vi)**2
    return V**(0.5)

def proj_vect(v,w):
    """ proj_vect(v,w) give the projection of vector v over w"""
    W = modeV(w)**2
    sc = 0.0
    for i in range(len(v)):
        sc = sc + float(v[i])*np.float32(w[i])
    proj = np.zeros(len(v))
    for i in range(len(v)):
        proj[i] = np.float32((sc/W))*w[i]
    return  proj

def intersect_lines(m1,n1,m2,n2):
    np.seterr(all =None)
    x = (n2-n1)/(m1-m2)
    return np.array([x, 0.0, x*m1+n1])

def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)

def openFolder(foldername):
    folder = './drawings_/'+ foldername
    createFolder(folder)
    os.chdir(folder)

def check_out_bound( point, sizeXYZ, equal = True):
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
    Nold, index, newp = N, [], []
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
    faces= []
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

def gen_thickness_data(self):
    N, coord, t, sizeXYZ, edges = self.N(), self.nodes, self.t, self.sizeXYZ, self.edges
    coco, ind, faces= [[],[],[]]
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
        faces[-1][-1] = 0
        edges = edges[:-1]
    else:
        c, index, N, Nold = create_points_t(coord[i],coord[i+1], 2*coord[i+1]-coord[i], t, sizeXYZ, N)
        coco.extend(c)
        ind.extend(index)
        ind, edges, faces = gen_t_faces(ind, edges, faces, N, Nold, i+1)
    self.nodes = np.concatenate((coord, np.reshape(np.asarray(coco),(len(coco),3))))
    self.edges, self.faces = edges, faces

def gen_mesh(self, meshSize):
    def do_trimming(grid, bound_node):
        gridd=[]
        for facet in grid:
            if not(all(f in bound_node for f in facet)):
                gridd.append(facet)
            else:
                pass
        gridd = np.asanyarray(gridd)
        return gridd
    N = self.N()
    minsize, maxsize = [float(m) for m in meshSize]
    old_edges = copy.deepcopy(self.edges)
    newc = []
    k =-1
    self.edges = [[a] for a in range(N)]
    for edge in old_edges:
        for no in edge[1:]:
            nedges = []
            ni=edge[0]
            vect = self.nodes[no]-self.nodes[ni]
            dist = np.linalg.norm(vect)
            div = 1
            while ((dist>maxsize) and ((dist/(div+1)>=minsize) and (dist/(div)>maxsize))):
                div = div + 1
            newcoord = self.nodes[ni]
            # Add new edges and new nodes to the data
            for i in range(1,div):
                k = k + 1
                newcoord = newcoord + vect*(1./div)
                newc.append(newcoord)
                nedges.append(k+N)
                self.edges.append([k+N])
            if len(self.bound_node) !=0:
                if (ni in self.bound_node) and (no in self.bound_node):
                    # Add new bound nodes to the bound node set
                    self.bound_node = self.bound_node.union(set(nedges))
            if len(nedges) == 0:
                self.edges[ni].append(no)
            if len(nedges) != 0:
                self.edges[ni].append(nedges[0])
                self.edges[no].append(nedges[-1])
            if len(nedges)>=2:
                for nn in nedges[0:-1]:
                    self.edges[nn].append(nn+1)
    # Add new nodes to the node array
    self.nodes = np.concatenate((self.nodes, np.reshape(newc,(len(newc),3))))
    # TRIANGULATE faces if needed
    if (np.sum(self.faces) != 0):
        triang = Delaunay(self.nodes[:,[0,2]])
        self.faces = do_trimming(triang.simplices, self.bound_node)


def get_unique_sorted(column_matrix):
    l_cm = len(column_matrix)
    column_matrix=np.reshape(column_matrix,( l_cm,1))
    index_sort = np.argsort(column_matrix,axis = 0)
    unique_values, uni_sort_ind = np.unique(np.reshape(column_matrix[index_sort],(l_cm,1)),axis=0,return_inverse = True)
    c_mat_us = uni_sort_ind[np.argsort(index_sort,axis = 0)]
    return unique_values, c_mat_us

def check_Dim_Mult(coord, sizeDim):
    delta, idim = get_unique_sorted(coord)
    if delta[-1] != 0:
        coord = np.matrix(coord*sizeDim/float(delta[-1]))
        delta = float(sizeDim/delta[-1])*delta
    n_L = int(np.shape(delta)[0])
    return coord, delta, idim, n_L
def check_dimensions(self):
    for dimension in [0,1,2]:
        if dimension not in self.RVE.dim:
            self.n_G[dimension] = 1
def do_pregrid(self):
#    print('do_pregrid( nodes, sizeXYZ F) // n_L, ii, self.delmax, delta G')
# Generating data to prepare the mesh
    ii = []
    n_L = []
    delta = []
    self.delmax = []
    for d in [0,1,2]:
        ## AQUIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
        self.RVE.nodes[:,d], delt, idim, n_ = check_Dim_Mult(self.RVE.nodes[:,d], self.RVE.sizeXYZ[0,d])
        delta.append(delt)
        self.delmax.append(float(max(delt)))
        ii.append(idim)
        n_L.append(n_)
    self.delmax = np.array(self.delmax)
    del(delt, n_, idim)
    self.n_L, self.ii, self.delmax, self.delta = n_L, ii, self.delmax, delta

def do_translate(tobetranslated, dictionary):

    assert (type(tobetranslated) == list), 'Please insert a list object'
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
    return translation


def do_shapeME(self, coord, ii):
    def shape_(self,r):
        return np.zeros([1,3])

    def shape_random(self,r, coef_rand =0.625):
        return np.random.uniform(-self.delmax/coef_rand,self.delmax/self.coef_rand,(1,3))

    def shape_progressive(self,r):
        try:
            return np.array([float((self.I[0]+ii[r,0])**2), float((self.I[1]+ii[r,1])**2),  float((self.I[2]+ii[r,2])**2)])
        except:
            try:
                return np.array([float((self.I[0]+ii[r,0])**2), 0.0,  float((self.I[2]+ii[r,2])**2)])
            except:
                return np.array([float((self.I[0]+ii[r,0])**2),0.0,0.0])

    def shape_flower(self,r, a = 5.0, b = 0.3):
        self.TOT_angleX, self.TOT_angleY = 5*math.pi, 1.0
        xt = (a+z)*math.cos(self.TOT_angleX*self._fract_angleX)*math.sin(self.TOT_angleY*self._fract_angleY)-x
        yt = (a+z)*math.sin(self.TOT_angleX*self._fract_angleX)*math.sin(self.TOT_angleY*self._fract_angleY)-y
        zt = (a+z)*(math.cos(self.TOT_angleY*self._fract_angleY)+math.log(math.tan(0.5*self.TOT_angleY*self._fract_angleY)+0.1))+b*self.TOT_angleX*self._fract_angleX-z
        return np.array([xt,yt,zt])

    def shape_ring(self,r, R = 25.0, Ri = 3.0):
        self.TOT_angleX, self.TOT_angleY = 2.0*math.pi, 2.0*math.pi
        xt = math.cos(self.TOT_angleX*self._fract_angleX)*(R+(Ri+z)*math.cos(self.TOT_angleY*self._fract_angleY))-x
        yt = math.sin(self.TOT_angleX*self._fract_angleX)*(R+(Ri+z)*math.cos(self.TOT_angleY*self._fract_angleY))-y
        zt = (Ri+z)*math.sin(self.TOT_angleY*self._fract_angleY)-z
        return np.array([xt,yt,zt])

    def shape_galaxy(self,r,aa= 3.5, a=5):
        self.TOT_angleX, self.TOT_angleY = 2*math.pi,4*math.pi
        xt = (a+z)*(aa+math.cos(self.TOT_angleY*self._fract_angleY/2)*math.sin(self.TOT_angleX*
                self._fract_angleX)-math.sin(self.TOT_angleY*self._fract_angleY/2)*math.sin(2*
                                     self.TOT_angleX*self._fract_angleX))*math.cos(self.TOT_angleY*self._fract_angleY) - x
        yt = (a+z)*(aa+math.cos(self.TOT_angleY*self._fract_angleY/2)*math.sin(self.TOT_angleX*
                self._fract_angleX)-math.sin(self.TOT_angleY*self._fract_angleY/2)*math.sin(2*
                                     self.TOT_angleX*self._fract_angleX))*math.sin(self.TOT_angleY*self._fract_angleY) -y
        zt = (a+z)*(math.sin(self.TOT_angleY*self._fract_angleY/2)*math.sin(self.TOT_angleX*
                     self._fract_angleX)+math.cos(self.TOT_angleY*self._fract_angleY/2)*
                     math.sin(2*self.TOT_angleX*self._fract_angleX)) -  z
        return np.array([xt,yt,zt])

    def shape_triplefris(self,r):
        self.TOT_angleX, self.TOT_angleY = math.pi, 3*math.pi
        xt = (1+z)*((math.cos(self.TOT_angleX*self._fract_angleX)*((1/3)*math.sqrt(2.0)*
                      math.cos(self.TOT_angleX*self._fract_angleX)*math.cos(2*self.TOT_angleY*self._fract_angleY) + (2/3)*
                      math.sin(self.TOT_angleX*self._fract_angleX)*math.cos(self.TOT_angleY*self._fract_angleY)))/(1 - math.sqrt(2.0)*
                      math.sin(self.TOT_angleX*self._fract_angleX)*math.cos(self.TOT_angleX*self._fract_angleX)*
                      math.sin(3*self.TOT_angleY*self._fract_angleY))) - x
        yt = (1+z)*((math.cos(self.TOT_angleX*self._fract_angleX)*((1/3)*math.sqrt(2.0)*math.cos(self.TOT_angleX*
                      self._fract_angleX)*math.sin(2*self.TOT_angleY*self._fract_angleY) -(2/3)*
                      math.sin(self.TOT_angleX*self._fract_angleX)*math.cos(self.TOT_angleY*self._fract_angleY)))/ (1 - math.sqrt(2.0)*
                      math.sin(self.TOT_angleX*self._fract_angleX)*math.cos(self.TOT_angleX*self._fract_angleX)*
                      math.sin(3*self.TOT_angleY*self._fract_angleY))) - y
        zt = (1+z)*((math.cos(self.TOT_angleX*self._fract_angleX)*math.cos(self.TOT_angleX*self._fract_angleX))/(1 - math.sqrt(2.0)*
              math.sin(self.TOT_angleX*self._fract_angleX)*math.cos(self.TOT_angleX*self._fract_angleX)*math.sin(3*
                      self.TOT_angleY*self._fract_angleY)) - 1) - z
        return np.array([xt,yt,zt])

    def shape_quadratic(self,r):
        xt = -0.5*x
        yt = 0
        zt = (0.005*(x-5)**2)
        return np.array([xt,yt,zt])

    def shape_shell(self,r):
        xt = 0
        yt = 0
        zt = (0.05*(x-self.delmax[0]*(self.n_G[0]*(self.n_L[0]-1))/2.0)**2) + (0.05*(y-self.delmax[1]*(self.n_G[1]*(self.n_L[1]-1))/2.0)**2)
        return np.array([xt,yt,zt])

    def shape_circle(self,r, Ri = 3.0):
        self.TOT_angleX = 2*math.pi
        return shape_cylinder(self,r)

    def shape_cylinder(self,r,Ri = 3.0):

        self.TOT_angleX = 2*math.pi
        xt = (z+Ri)*math.cos(self._fract_angleX*self.TOT_angleX)-x
        yt = 0
        zt = (z+Ri)*math.sin(self._fract_angleX*self.TOT_angleX)-z
        return np.array([xt,yt,zt])

    def shape_catenoid(self,r,Ri = 5.0):
        self.TOT_angleX, self.TOT_angleY = 2*math.pi, math.pi
        xt = (z+Ri)*math.cos(self.TOT_angleX*self._fract_angleX-math.pi)*math.cosh(self.TOT_angleY*self._fract_angleY/Ri)-x
        yt = (z+Ri)*math.sin(self.TOT_angleX*self._fract_angleX-math.pi)*math.cosh(self.TOT_angleY*self._fract_angleY/Ri)-y
        zt = (self.TOT_angleY*self._fract_angleY)-z
        return np.array([xt,yt,zt])
#
    def shape_sphere(self, r, Ri=3.0):
        self.TOT_angleX, self.TOT_angleY = 2*math.pi, math.pi
        xt = (z+Ri)*math.cos(self.TOT_angleX*self._fract_angleX)*math.sin(self.TOT_angleY*self._fract_angleY)-x
        yt = (z+Ri)*math.cos(self.TOT_angleY*self._fract_angleY)-y
        zt = (z+Ri)*math.sin(self.TOT_angleX*self._fract_angleX)*math.sin(self.TOT_angleY*self._fract_angleY)-z
        return np.array([xt,yt,zt])

    def shape_NURBS(self,r):
        xt = 0
        yt = 0
        zt = (0.005*(x-self.delmax[0]*(self.n_G[0]*(self.n_L[0]-1))/2.0)*(y-self.delmax[1]*(self.n_G[1]*(self.n_L[1]-1))/2.0))**3
        +(0.005*(x-self.delmax[0]*(self.n_G[0]*(self.n_L[0]-1))/2.0)*(y-self.delmax[1]*(self.n_G[1]*(self.n_L[1]-1))/2.0))**2
        +(0.005*(x-self.delmax[0]*(self.n_G[0]*(self.n_L[0]-1))/2.0)*(y-self.delmax[1]*(self.n_G[1]*(self.n_L[1]-1))/2.0))
        return np.array([xt,yt,zt])

    r = -1

    for x,y,z in coord:
        r = r+1
        self._fract_angleX = float((float(self.I[0]*(self.n_L[0]-1)+ii[r,0]))/(self.n_G[0]*(self.n_L[0]-1)))
        if self.n_L[1] != 1: self._fract_angleY = float((float(self.I[1]*(self.n_L[1]-1)+ii[r,1]))/(self.n_G[1]*(self.n_L[1]-1)))
        self.nodes[self.num_r[0,r]] = coord[r] + np.around(eval('shape_'+str(self.shape)+'(self,r)'),12)
        flag = False
        for n in self.RVE.dim:
            if (self.I[n] == 0)&(ii[r,n] == 0):
                for symmetric in self.bound[n]:
                    if np.all(np.power(self.nodes[self.num_r[0,r]]-self.nodes[symmetric,:],2)<1e-6)&(flag==False):
                        self.num_r[0,r] = symmetric
                        flag = True
                if (flag==False):
                    self.boundary[2*n].add(self.num_r[0,r])
                    if self.TOT_angleX !=None: self.bound[n].append(self.num_r[0,r])
            if (self.I[n] == (self.n_G[n]-1))&(ii[r,n] == self.n_L[n]-1)&(flag==False):
                for symmetric in reversed(self.bound[n]):
                    if np.all(np.power(self.nodes[self.num_r[0,r]]-self.nodes[symmetric,:],2)<1e-6)&(flag==False):
                        self.num_r[0,r] = symmetric
                        flag = True
                if flag == False:
#                             print(4,r)
                    self.boundary[2*n+1].add(self.num_r[0,r])
                    if self.TOT_angleY !=None: self.bound[n].append(self.num_r[0,r])

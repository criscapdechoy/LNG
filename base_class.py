#!/usr/bin/env python
# -----------------------------------------------------------------------------
# base.py
# -----------------------------------------------------------------------------
import timeit
import base_func_class as base_func
import numpy as np
import math

def main(TypeElem = 8, shape ='cylinder', n_G = [30,30,2], meshSize=[0.08,0.1], sizeXYZ=(1.0,1.0,1.0), 
         nodeR = 0.0, ExtractGeom = True, View = False):
    start_time = timeit.default_timer()
    #t=2*0.7*sizeXYZ[2]/3
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
    
    # Basic element definition: 
    LS = base_func.get_basic_elem( TypeElem )
    # Preparing geometry and topology
    LS, scale = base_func.get_input(LS, sizeXYZ)
    # Generate edge width, "t", if given:(N, t, sizeXYZ, coord, edges, faces=[], show = False):
    LS = base_func.get_faces( LS, sizeXYZ, show = False)
    # Mesh pattern generation 
    LS  = base_func.gen_mesh( LS, sizeXYZ, meshSize, DoMesh = False)
    # Generating global nodes minsize=t/2, maxsize=t*(2**0.5)/2
    LS, ind_hash, num, delmax = base_func.gen_nodegrid(LS, n_G, sizeXYZ, shape)
    # Generating global edges
    LS.edges = base_func.gen_edges(num, ind_hash, LS.edges)
    # Generating global faces
    LS.faces = base_func.gen_faces(num, ind_hash, LS.faces)
    # Delete random edges
    delnod = np.random.randint(LS.N, size=int(0.1*LS.N))
    for con in delnod:
        while len(LS.edges[con])==1:
            con = np.random.randint(LS.N, size=1)[0]
        del(LS.edges[con][np.random.randint(1,len(LS.edges[con]),size=1)[0]])
    # Creating CADfile
    filename = shape +'geom_t'+ str(TypeElem)+ str(LS.nodeR) +'ms'+ ''.join(map(str, n_G[:]))
    base_func.gen_CAD(filename, LS, [np.zeros(3), delmax*n_G], ExtractGeom, View, foldername = shape)  
    elapsed = timeit.default_timer() - start_time
    print("%s nodes" % LS.N )
    print("elapsed time: %s seconds" % elapsed )
    return LS
    
if __name__== "__main__":
    LS = main()
    pass


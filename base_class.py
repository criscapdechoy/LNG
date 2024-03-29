#!/usr/bin/env python
# -----------------------------------------------------------------------------
# base.py
# -----------------------------------------------------------------------------
import timeit
import math
import numpy as np
from LNG_main import RVE, LatticeStructure

def build(TypeElem = 8, shape ='cylinder', n_G = [1,1,2], meshSize=[], sizeXYZ=[1.0,1.0,1.0], Export = True, View = False):
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
    Elem = RVE(5, sizeXYZ=[5.0,5.0,5.0])
    # Mesh generation 
    Elem.gen_mesh(meshSize = [])
    # Initializing lattice structure net
    LS = LatticeStructure( Elem, [10,5,1], shape )

    # Generating global nodes 
    LS.gen_nodegrid()
    # Generating global edges
    LS.gen_edges()
    # Generating global faces
    LS.gen_faces()
    # Delete random edges
    delnod = np.random.randint(len(LS.mesh), size=int(0.1*len(LS.mesh)))
    for con in delnod:
        while len(LS.edges[con])==1:
            con = np.random.randint(len(LS.mesh), size=1)[0]
        del(LS.edges[con][np.random.randint(1,len(LS.edges[con]),size=1)[0]])
    # Creating CADfile
    if Export == True: LS.gen_CAD()
    if View == True: LS.show()
    elapsed = timeit.default_timer() - start_time
    print("%s nodes" % len(LS.mesh) )
    print("elapsed time: %s seconds" % elapsed )
    return LS
    
if __name__== "__main__":
    build()

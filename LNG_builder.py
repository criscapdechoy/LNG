#!/usr/bin/env python
# -----------------------------------------------------------------------------
# LNG_builder.py
# -----------------------------------------------------------------------------
import timeit
import math
import numpy as np
from LNG_main import RVE, LatticeStructure

def main(TypeElem = 8, shape ='galaxy', n_G = [20,60,1], meshSize=[], sizeXYZ=[1,1,1], Export = True, View = False):
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
    start_time = timeit.default_timer()

    Elem = RVE(TypeElem, sizeXYZ)
    # Mesh generation
    Elem.gen_mesh(meshSize)
    # Initializing lattice structure net
    LS = LatticeStructure( Elem, n_G, shape )
    # Generating global nodes
    LS.gen_nodegrid()
    # Generating global edges
    LS.gen_edges()
    # Generating global faces
    LS.gen_faces()
    elapsed = timeit.default_timer() - start_time
    # Delete random edges
#    delnod = np.random.randint(len(LS.nodes), size=int(0.1*len(LS.nodes)))
#    for con in delnod:
#        while len(LS.edges[con])==1:
#            con = np.random.randint(len(LS.nodes), size=1)[0]
#        del(LS.edges[con][np.random.randint(1,len(LS.edges[con]),size=1)[0]])
#    # Creating CADfile
#    if Export == True: LS.gen_CAD()
#    if View == True: LS.show(lim=[[-20,-20,-20],[20,20,20]])

    print("{} nodes, {} edges and {} faces. ".format(len(LS.nodes), LS.num_edges, len(LS.faces)) )
    print("elapsed time: %s seconds" % elapsed )
    return LS

if __name__== "__main__":
    LS = main()

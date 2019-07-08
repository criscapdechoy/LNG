#!/usr/bin/env python
# -----------------------------------------------------------------------------
# base.py
# -----------------------------------------------------------------------------
import timeit
import base_func_flowchart as base_func
import numpy as np
import math

def main(TypeElem = 6, shape ='', n_G = [15,15,15], meshSize=[0.08,0.1], sizeXYZ=[1.0,1.0,1.0], 
         t = 0.0, ExtractGeom = True, View = False):
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
    nodes, edges, faces, sym, dim = base_func.get_basic_elem( TypeElem )
    # Check dimensions coherency
#    n_G = base_func.check_dimensions(dim, n_G)
    # Preparing geometry and topology
    nodes, edges, faces, N , scale = base_func.get_input(nodes, edges, faces, sizeXYZ, dim)
    # Generate edge width, "t", if given:(N, t, sizeXYZ, coord, edges, faces=[], show = False):
    nodes, edges, faces, N, bound_node = base_func.get_faces( N, t,  sizeXYZ, nodes, edges, faces, show = False)
    # Mesh pattern generation 
    nodes, edges, faces, N, bound_node = base_func.gen_mesh( N, nodes, edges, faces, bound_node, sizeXYZ, meshSize, DoMesh = False)
    # Generating global nodes minsize=t/2, maxsize=t*(2**0.5)/2
    mesh, ind_hash, num, delmax, Boundary = base_func.gen_nodegrid(nodes, N, sym, n_G, dim, sizeXYZ, shape)
    # Generating global edges
    edges, Nedges = base_func.gen_edges(num, ind_hash, edges)
    total_length = 0.0
    for edge in edges:
        for eo in edge[1:]:
            ei = edge[0]
            ledge = base_func.do_modeV(mesh[ei]-mesh[eo])
            total_length = total_length + ledge
    total_volume = (total_length-len(mesh))*0.01**2
    # Generating global faces
    faces = base_func.gen_faces(num, ind_hash, faces)

    # Creating CADfile
    filename = shape +'geom_t'+ str(TypeElem)+ str(t) +'ms'+ ''.join(map(str, n_G[:]))
    base_func.gen_CAD(filename, dim, mesh, [np.zeros(3),delmax*n_G], edges, faces, t, ExtractGeom, View, foldername = shape)  
    elapsed = timeit.default_timer() - start_time
    print("%s nodes and %s edges" % (len(mesh),Nedges) )
    print("elapsed time: %s seconds" % elapsed )
    return mesh, edges, faces, Boundary
    
if __name__== "__main__":
    mesh, edges, faces, Boundary = main()
    pass


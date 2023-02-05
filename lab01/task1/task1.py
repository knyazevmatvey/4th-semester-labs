import gmsh
import sys

gmsh.initialize()

gmsh.model.add("task1")


# defining parameters for geometry
s = 0.9     # big torus radius
r = 0.3     # inner radius of the tube
R = 0.5     # outer radius of the tube

# creates torus with center (0,0,0), with radii (s,R), (s,r)
tor2 = gmsh.model.occ.addTorus(0, 0, 0, s, r)
tor1 = gmsh.model.occ.addTorus(0, 0, 0, s, R)

# finds difference between volumes = tor1 - tor2 (argument '3' corresponds to dims)
gmsh.model.occ.cut([(3, tor1)], [(3, tor2)], 3)

# changing mesh size
gmsh.option.setNumber("Mesh.MeshSizeFactor", 0.2)

gmsh.model.occ.synchronize()
gmsh.model.mesh.generate(3)
gmsh.write("task1.msh")

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
gmsh.finalize()

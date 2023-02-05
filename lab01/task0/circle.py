import gmsh
import sys

gmsh.initialize()

gmsh.model.add("circle")

lc = 0.01 # mesh parameter ??? (idk)

# defining points
x = 0
y = 0
z = 0
r = 0.1
p0 = gmsh.model.geo.addPoint(x, y, z, lc)
p1 = gmsh.model.geo.addPoint(x+r, y, z, lc)
p2 = gmsh.model.geo.addPoint(x, y+r, z, lc)
p3 = gmsh.model.geo.addPoint(x-r, y, z, lc)
p4 = gmsh.model.geo.addPoint(x, y-r, z, lc)

# defining circle arcs
c1 = gmsh.model.geo.addCircleArc(p1, p0, p2)
c2 = gmsh.model.geo.addCircleArc(p2, p0, p3)
c3 = gmsh.model.geo.addCircleArc(p3, p0, p4)
c4 = gmsh.model.geo.addCircleArc(p4, p0, p1)

# defining surface
gmsh.model.geo.addCurveLoop([c1, c2, c3, c4], 1)
gmsh.model.geo.addPlaneSurface([1], 1)

gmsh.model.geo.synchronize()
gmsh.model.mesh.generate(2)
gmsh.write("circle.msh")

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
gmsh.finalize()

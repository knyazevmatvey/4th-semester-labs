import gmsh
import sys

gmsh.initialize()

gmsh.model.add("cyl")

lc = 0.01 # mesh parameter ??? (idk)

# defining points for lower circle
x = 0
y = 0
z = 0
r = 0.1
p0 = gmsh.model.geo.addPoint(x, y, z, lc)
p1 = gmsh.model.geo.addPoint(x+r, y, z, lc)
p2 = gmsh.model.geo.addPoint(x, y+r, z, lc)
p3 = gmsh.model.geo.addPoint(x-r, y, z, lc)
p4 = gmsh.model.geo.addPoint(x, y-r, z, lc)

# defining points for upper circle
x = 0
y = 0
z = 0.3
r = 0.1
p5 = gmsh.model.geo.addPoint(x, y, z, lc)
p6 = gmsh.model.geo.addPoint(x+r, y, z, lc)
p7 = gmsh.model.geo.addPoint(x, y+r, z, lc)
p8 = gmsh.model.geo.addPoint(x-r, y, z, lc)
p9 = gmsh.model.geo.addPoint(x, y-r, z, lc)

# defining lower circle arcs
c1 = gmsh.model.geo.addCircleArc(p1, p0, p2)
c2 = gmsh.model.geo.addCircleArc(p2, p0, p3)
c3 = gmsh.model.geo.addCircleArc(p3, p0, p4)
c4 = gmsh.model.geo.addCircleArc(p4, p0, p1)

# defining upper circle arcs
c5 = gmsh.model.geo.addCircleArc(p6, p5, p7)
c6 = gmsh.model.geo.addCircleArc(p7, p5, p8)
c7 = gmsh.model.geo.addCircleArc(p8, p5, p9)
c8 = gmsh.model.geo.addCircleArc(p9, p5, p6)

# connecting lower and upper points
c9 = gmsh.model.geo.addLine(p1, p6)
c10 = gmsh.model.geo.addLine(p3, p8)

# defining lower circle surface
gmsh.model.geo.addCurveLoop([c1, c2, c3, c4], 1)
gmsh.model.geo.addPlaneSurface([1], 1)

# defining upper circle surface
gmsh.model.geo.addCurveLoop([c5, c6, c7, c8], 2)
gmsh.model.geo.addPlaneSurface([2], 2)

# here i had to use .occ (instead of .geo) because idk how to do
# non-plane surface otherwise
circ1 = gmsh.model.occ.addCircle(0, 0, 0, 0.1)
circ2 = gmsh.model.occ.addCircle(0, 0, 0.3, 0.1)
gmsh.model.occ.addCurveLoop([circ1], 3)
gmsh.model.occ.addCurveLoop([circ2], 4)
gmsh.model.occ.addThruSections([3, 4], 3)
# defining cylinder walls
#gmsh.model.geo.addCurveLoop([9, 5, 6, -10, -2, -1], 3)
#gmsh.model.geo.addCurveLoop([-10, 3, 4, 9, -8, -7], 4)
#gmsh.model.geo.addRuledSurface([3], 3)
#gmsh.model.geo.addRuledSurface([4], 4)

gmsh.model.occ.synchronize()
gmsh.model.geo.synchronize()
gmsh.model.mesh.generate(2)
gmsh.write("cyl.msh")

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
gmsh.finalize()











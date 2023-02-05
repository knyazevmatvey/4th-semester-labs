import gmsh
import sys

gmsh.initialize()

gmsh.model.add("task0")

# defining points
lc = 1e-2
gmsh.model.geo.addPoint(0, 0, 0, lc, 1)
gmsh.model.geo.addPoint(0.1, 0, 0, lc, 2)
gmsh.model.geo.addPoint(0.1, 0.1, 0, lc, 3)
gmsh.model.geo.addPoint(0, 0.1, 0, lc, 4)
gmsh.model.geo.addPoint(0, 0, 0.1, lc, 5)
gmsh.model.geo.addPoint(0.1, 0, 0.1, lc, 6)
gmsh.model.geo.addPoint(0.1, 0.1, 0.1, lc, 7)
gmsh.model.geo.addPoint(0, 0.1, 0.1, lc, 8)

# defining lines
for i in range(3):
    gmsh.model.geo.addLine(i+1, i+2, i+1)
gmsh.model.geo.addLine(4, 1, 4)
for i in range(3):
    gmsh.model.geo.addLine(i+5, i+6, i+5)
gmsh.model.geo.addLine(8, 5, 8)
for i in range(4):
    gmsh.model.geo.addLine(i+1, i+5, i+9)

# defining surfaces
gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 1)
gmsh.model.geo.addCurveLoop([5, 6, 7, 8], 2)
gmsh.model.geo.addCurveLoop([2, 11, -6, -10], 3)
gmsh.model.geo.addCurveLoop([1, 10, -5, -9], 4)
gmsh.model.geo.addCurveLoop([9, -8, -12, 4], 5)
gmsh.model.geo.addCurveLoop([3, 12, -7, -11], 6)
for i in range(6):
    gmsh.model.geo.addPlaneSurface([i+1], i+1)

loop = gmsh.model.geo.addSurfaceLoop([i+1 for i in range(6)])
gmsh.model.geo.addVolume([loop])

gmsh.model.geo.synchronize()

gmsh.model.mesh.generate(3)

gmsh.write("cube.msh")

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()

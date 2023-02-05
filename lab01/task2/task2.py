import gmsh
import math
import os
import sys

gmsh.initialize()

# loading stl mesh
path = os.path.dirname(os.path.abspath(__file__))
gmsh.merge(os.path.join(path, 'cat.stl'))
'''
# classifying surfaces (so that angles of one surface are less than 'angle')
angle = 40
forceParametrizablePatches = False
includeBoundary = True
curveAngle = 180

gmsh.model.mesh.classifySurfaces(angle * math.pi / 180., includeBoundary,
                                 forceParametrizablePatches,
                                 curveAngle * math.pi / 180.)
'''
# this line ignores errors in createGeometry
#gmsh.option.setNumber("General.AbortOnError", 0)
#gmsh.model.mesh.createGeometry()

# creating a volume from all surfaces
s = gmsh.model.getEntities(2)
l = gmsh.model.geo.addSurfaceLoop([s[i][1] for i in range(len(s))])
gmsh.model.geo.addVolume([l])

gmsh.model.geo.synchronize()

print('before generating')

# changing mesh size
gmsh.option.setNumber("Mesh.MeshSizeFactor", 0.2)

gmsh.model.mesh.generate(3)
gmsh.write('task2.msh')

print('after generating')

# Launch the GUI to see the results:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()


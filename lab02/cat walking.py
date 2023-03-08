import vtkmodules.all as vtk
import numpy as np
import gmsh
import os
import math

# tail parameters
tail_separator_x = -0.31
tail_separator_y = -0.02
tail_omega = 0.2
type_tail = 1

# leg parameters
leg_separator_y = -0.14
right_left_leg_separator_z = 0.20
forward_back_leg_separator = -0.15
leg_period = 3
leg_ampl = 0.40
type_leg = 2
type_forward_right_leg = 3
type_forward_left_leg = 4
type_back_right_leg = 5
type_back_left_leg = 6

forward_leg_axis = [-0.06, -0.12, 0]
back_leg_axis = [-0.25, -0.08, 0]


# time step
tau = 0.05

# overall velocity
v_average = 0.1 * leg_ampl

# rotation axis parameters (tail)
axis_y = 0.00
axis_z = 0.20


# Класс расчётной сетки
class CalcMesh:

    # Конструктор сетки, полученной из stl-файла
    def __init__(self, nodes_coords, tetrs_points):
        # defining time
        self.time = 0

        # 3D-сетка из расчётных точек
        # Пройдём по узлам в модели gmsh и заберём из них координаты
        self.nodes = np.array([nodes_coords[0::3], nodes_coords[1::3], nodes_coords[2::3]])
        self.number_of_nodes = len(self.nodes[0])

        # Модельная скалярная величина распределена как-то вот так
        self.smth = np.power(self.nodes[0, :], 2) + np.power(self.nodes[1, :], 2)

        # separating nodes into types
        self.node_type = np.zeros(shape=(int(len(nodes_coords) / 3)))
        for i in range(self.number_of_nodes):
            if (self.nodes[0, i] < tail_separator_x) and (self.nodes[1, i] > tail_separator_y):
                self.node_type[i] = type_tail

            if self.nodes[1, i] < leg_separator_y:
                if self.nodes[0, i] > forward_back_leg_separator:
                    if self.nodes[2, i] > right_left_leg_separator_z:
                        self.node_type[i] = type_forward_right_leg
                    else:
                        self.node_type[i] = type_forward_left_leg
                else:
                    if self.nodes[2, i] > right_left_leg_separator_z:
                        self.node_type[i] = type_back_right_leg
                    else:
                        self.node_type[i] = type_back_left_leg


        # Тут может быть скорость, но сейчас здесь нули
        self.velocity = np.zeros(shape=(3, int(len(nodes_coords) / 3)), dtype=np.double)

        # Пройдём по элементам в модели gmsh
        self.tetrs = np.array([tetrs_points[0::4], tetrs_points[1::4], tetrs_points[2::4], tetrs_points[3::4]])
        self.tetrs -= 1

    def rotate(self, direction, point, omega, moving_type):
        """
        the method rotates points of type 'moving_type' around axis going through 'point' parallel to
        'dir' with angular velocity 'omega'
        direction -- an array with axis coordinates (presumably a unit vector)
        point -- an array with coordinates of a point on the axis
        moving_type -- type of points to be moved
        omega -- the angular velocity
        """
        for i in range(self.number_of_nodes):
            if self.node_type[i] == moving_type:
                self.velocity[0, i] = omega * (direction[1] * (self.nodes[2, i] - point[2])
                                               - direction[2] * (self.nodes[1, i] - point[1]))
                self.velocity[1, i] = omega * (direction[2] * (self.nodes[0, i] - point[0])
                                               - direction[0] * (self.nodes[2, i] - point[2]))
                self.velocity[2, i] = omega * (direction[0] * (self.nodes[1, i] - point[1])
                                               - direction[1] * (self.nodes[0, i] - point[0]))

    # the method recalculates velocity field
    def calc_velocities(self):
        # Обнуляем все скорости перед вычислением
        for i in range(self.number_of_nodes):
            self.velocity[0, i] = 0
            self.velocity[1, i] = 0
            self.velocity[2, i] = 0

        # Двигаем лапы
        k = math.cos(2 * math.pi * self.time / leg_period) * leg_ampl
        self.rotate([0, 0, -1], np.add(forward_leg_axis, [self.time * v_average, 0, 0]), k, type_forward_right_leg)
        self.rotate([0, 0, -1], np.add(forward_leg_axis, [self.time * v_average, 0, 0]), -k, type_forward_left_leg)
        self.rotate([0, 0, -1], np.add(back_leg_axis, [self.time * v_average, 0, 0]), k, type_back_right_leg)
        self.rotate([0, 0, -1], np.add(back_leg_axis, [self.time * v_average, 0, 0]), -k, type_back_left_leg)

        # for i in range(self.number_of_nodes):
        #    if self.node_type[i] == type_tail:
        #       self.velocity[1, i] = -tail_omega * (self.nodes[2, i] - axis_z)
        #        self.velocity[2, i] = tail_omega * (self.nodes[1, i] - axis_y)

        # if self.node_type[i] == type_leg:
        #    self.velocity[0, i] = (leg_separator_y - self.nodes[1, i]) * math.cos(self.time / leg_period) * leg_ampl

        # if self.node_type[i] == type_forward_right_leg:
        #    k = math.cos(self.time / leg_period) * leg_ampl * leg_omega
        #    self.velocity[0, i] = (self.nodes[1, i] - forward_leg_axis_y) * k
        #    self.velocity[1, i] = -(self.nodes[0, i] - forward_leg_axis_x) * k

    # Метод отвечает за выполнение для всей сетки шага по времени величиной tau
    def move(self, tau):
        # По сути метод просто двигает все точки c их текущими скоростями
        self.calc_velocities()
        self.nodes += self.velocity * tau
        self.time += tau

        # Добавляем равномерное движение
        for i in range(self.number_of_nodes):
            self.nodes[0, i] += tau * v_average

    # Метод отвечает за запись текущего состояния сетки в снапшот в формате VTK
    def snapshot(self, snap_number):
        # Сетка в терминах VTK
        unstructuredGrid = vtk.vtkUnstructuredGrid()
        # Точки сетки в терминах VTK
        points = vtk.vtkPoints()

        # Скалярное поле на точках сетки
        smth = vtk.vtkDoubleArray()
        smth.SetName("smth")

        # Векторное поле на точках сетки
        vel = vtk.vtkDoubleArray()
        vel.SetNumberOfComponents(3)
        vel.SetName("vel")

        # Обходим все точки нашей расчётной сетки
        # Делаем это максимально неэффективным, зато наглядным образом
        for i in range(0, len(self.nodes[0])):
            # Вставляем новую точку в сетку VTK-снапшота
            points.InsertNextPoint(self.nodes[0, i], self.nodes[1, i], self.nodes[2, i])
            # Добавляем значение скалярного поля в этой точке
            smth.InsertNextValue(self.smth[i])
            # Добавляем значение векторного поля в этой точке
            vel.InsertNextTuple((self.velocity[0, i], self.velocity[1, i], self.velocity[2, i]))

        # Грузим точки в сетку
        unstructuredGrid.SetPoints(points)

        # Присоединяем векторное и скалярное поля к точкам
        unstructuredGrid.GetPointData().AddArray(smth)
        unstructuredGrid.GetPointData().AddArray(vel)

        # А теперь пишем, как наши точки объединены в тетраэдры
        # Делаем это максимально неэффективным, зато наглядным образом
        for i in range(0, len(self.tetrs[0])):
            tetr = vtk.vtkTetra()
            for j in range(0, 4):
                tetr.GetPointIds().SetId(j, self.tetrs[j, i])
            unstructuredGrid.InsertNextCell(tetr.GetCellType(), tetr.GetPointIds())

        # Создаём снапшот в файле с заданным именем
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetInputDataObject(unstructuredGrid)
        writer.SetFileName("lab02-step-" + str(snap_number) + ".vtu")
        writer.Write()


# Теперь придётся немного упороться:
# (а) построением сетки средствами gmsh,
# (б) извлечением данных этой сетки в свой код.
gmsh.initialize()

# Считаем STL
try:
    path = os.path.dirname(os.path.abspath(__file__))
    gmsh.merge(os.path.join(path, 'cat.stl'))
except:
    print("Could not load STL mesh: bye!")
    gmsh.finalize()
    exit(-1)

# no geometry

# Восстановим геометрию
'''
angle = 40
forceParametrizablePatches = False
includeBoundary = True
curveAngle = 180
gmsh.model.mesh.classifySurfaces(angle * math.pi / 180., includeBoundary, forceParametrizablePatches,
                                 curveAngle * math.pi / 180.)
gmsh.model.mesh.createGeometry()
'''

# Зададим объём по считанной поверхности
surfaces = gmsh.model.getEntities(2)
loop = gmsh.model.geo.addSurfaceLoop([surfaces[i][1] for i in range(len(surfaces))])
gmsh.model.geo.addVolume([loop])

gmsh.model.geo.synchronize()

# Зададим мелкость желаемой сетки
f = gmsh.model.mesh.field.add("MathEval")
gmsh.model.mesh.field.setString(f, "F", "4")
gmsh.model.mesh.field.setAsBackgroundMesh(f)

# Построим сетку
gmsh.model.mesh.generate(3)

# Теперь извлечём из gmsh данные об узлах сетки
nodeTags, nodesCoord, parametricCoord = gmsh.model.mesh.getNodes()

# И данные об элементах сетки тоже извлечём, нам среди них нужны только тетраэдры, которыми залит объём
GMSH_TETR_CODE = 4
tetrsNodesTags = None
elementTypes, elementTags, elementNodeTags = gmsh.model.mesh.getElements()
for i in range(0, len(elementTypes)):
    if elementTypes[i] != GMSH_TETR_CODE:
        continue
    tetrsNodesTags = elementNodeTags[i]

if tetrsNodesTags is None:
    print("Can not find tetra data. Exiting.")
    gmsh.finalize()
    exit(-2)

print("The model has %d nodes and %d tetrs" % (len(nodeTags), len(tetrsNodesTags) / 4))

# На всякий случай проверим, что номера узлов идут подряд и без пробелов
for i in range(0, len(nodeTags)):
    # Индексация в gmsh начинается с 1, а не с нуля. Ну штош, значит так.
    assert (i == nodeTags[i] - 1)
# И ещё проверим, что в тетраэдрах что-то похожее на правду лежит.
assert (len(tetrsNodesTags) % 4 == 0)

gmsh.finalize()

mesh = CalcMesh(nodesCoord, tetrsNodesTags)
mesh.snapshot(0)

# Делаем шаги по времени,
# на каждом шаге считаем новое состояние и пишем его в VTK
for i in range(1, 60):
    mesh.move(tau)
    mesh.snapshot(i)

print("finished")

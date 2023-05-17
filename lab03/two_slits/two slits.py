from fenics import *
import mshr

# parameters
b = 0.10
d = 0.30
wavelength = 0.1
k = 2*3.14/wavelength

# Create mesh and define function space
n=500
domain = mshr.Circle(Point(0.,0.),1.0,n)
mesh = mshr.generate_mesh(domain, n, "cgal")
print("generated mesh")

# Defining mixed Function space
V_real = FiniteElement("P", mesh.ufl_cell(), 1)
V_im = FiniteElement("P", mesh.ufl_cell(), 1)
V = FunctionSpace(mesh, V_real * V_im)

x = SpatialCoordinate(mesh)



# Define variational problem
(u_real, u_im) = TrialFunctions(V)
(v, tau) = TestFunctions(V)

# Defining different walls
mf = MeshFunction("size_t", mesh, 2)
tol = 1E-10

# Dirichlet
def left_wall(x, on_boundary):
    return x[0] <= tol
#class LeftWall(SubDomain):
#    def inside(self, x, on_boundary):
#        return x[0] <= tol

class K(UserExpression):
    def eval(self, value, x):
        "Set value[0] to value at point x"
        if (x[1] >= d/2-b/2) and (x[1] <= d/2+b/2):
            value[0] = 1
        elif (x[1] >= -d/2-b/2) and (x[1] <= -d/2+b/2):
            value[0]= 1
        else:
            value[0] = 0
src = K(degree=0)

class BoundaryOtherWalls(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (x[0] >= 0 - tol)

b_other = BoundaryOtherWalls()
b_other.mark(mf, 1)

zero = Constant("0")
one = Constant("1")

# Dirichlet condition
bc_1 = DirichletBC(V.sub(0), src, left_wall)
bx_2 = DirichletBC(V.sub(1), zero, left_wall)
bcs = [bc_1]


ds = Measure('ds', domain=mesh, subdomain_data=mf)

one = Constant("1")

# Main equation
F = (inner(grad(u_real), grad(v))*dx - k**2*v*u_real*dx + 
     inner(grad(tau), grad(u_im))*dx - k**2*tau*u_im*dx +
     v*k*u_im*ds(1) - tau*k*u_real*ds(1))
a = lhs(F)
L = rhs(F)

# Compute solution
u = Function(V)
solve(a == L, u, bcs)

(ans_re, ans_im) = u.split()
# Save solution to file in VTK format
vtkfile_real = File('two_slits/real.pvd')
vtkfile_real << ans_re

vtkfile_im = File('two_slits/im.pvd')
vtkfile_im << ans_im
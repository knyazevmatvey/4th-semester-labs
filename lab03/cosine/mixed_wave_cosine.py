from fenics import *

# Create mesh and define function space
mesh = UnitSquareMesh(100, 100)

# Defining mixed Function space
V_real = FiniteElement("P", mesh.ufl_cell(), 2)
V_im = FiniteElement("P", mesh.ufl_cell(), 2)
V = FunctionSpace(mesh, V_real * V_im)

x = SpatialCoordinate(mesh)

# Defining wave number (corresponding to n_refr = 1)
wavelength = 1.0
k = 2*3.14/wavelength

# Define variational problem
(u_real, u_im) = TrialFunctions(V)
(v, tau) = TestFunctions(V)



zero = Constant("0")
one = Constant("0")


# Dirichlet condition
def boundary(x, on_boundary):
    return on_boundary

u_D_re = Expression("cos(k*x[0])", k=k, degree=2)
u_D_im = Expression("sin(k*x[0])", k=k, degree=2)

bc_re = DirichletBC(V.sub(0), u_D_re, boundary)
bc_im = DirichletBC(V.sub(1), u_D_im, boundary)
bcs = [bc_re, bc_im]



src = Constant("0")
# impedance of outside
#z = k
z = 0
# Main equation
F = (inner(grad(u_real), grad(v))*dx - k**2*v*u_real*dx + 
     inner(grad(tau), grad(u_im))*dx - k**2*tau*u_im*dx)
a = lhs(F)
L = rhs(F)


# Compute solution
u = Function(V)
solve(a == L, u, bcs)
(ans_re, ans_im) = u.split()

# Save solution to file in VTK format
vtkfile_real = File('test_cosine/real.pvd')
vtkfile_real << ans_re

vtkfile_im = File('test_cosine/im.pvd')
vtkfile_im << ans_im

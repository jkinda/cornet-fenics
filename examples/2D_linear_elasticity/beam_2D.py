# -*- coding: utf-8 -*-
"""
Cantilever beam/plate problem in 2D elasticity (plane strain and plane stress)
uing a traditional displacement-based formulation

Created on Wed Feb 10 10:06:12 2016
@author: Jeremy Bleyer  (jeremy.bleyer@enpc.fr)
"""

from dolfin import *

# Geometry : rectangular beam of length L and height H
L, H = 10., 1.

# Discretization : (Nx x Ny) elements
Nx = 500
Ny = int(H/L*Nx)+1

# Elasticity parameters
E,nu = 100, 0.
# Weight density
rho_g = 1e-2

# 2D model : plane stress or plane strain
plane_stress = True

# Lamé coefficient for constitutive relation
mu = E/2./(1+nu)
lmbda_plane_strain = E*nu/(1+nu)/(1-2*nu)
lmbda_plane_stress = E*nu/(1-nu**2)

# Define epsilon and sigma operator
def eps(v):
    return sym(grad(v))
def sigma(v):
    dim = v.geometric_dimension()
    if dim == 2 and plane_stress:
        lmbda = lmbda_plane_stress
    else:
        lmbda = lmbda_plane_strain
    return lmbda*tr(eps(v))*Identity(dim) + 2.0*mu*eps(v)

# Create mesh
mesh = RectangleMesh(Point(0.,0.),Point(L,H), Nx, Ny, 'crossed')

# Unform vertical loading due to self-weight
f = Constant((0.,-rho_g))

# Define function space
V = VectorFunctionSpace(mesh, 'Lagrange', 2)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
a = inner(sigma(u), eps(v))*dx
l = inner(f, v)*dx

# Define left boundary x=0
def left(x, on_boundary):
    return near(x[0],0.)

# Fully clamped boundary condition
bc = DirichletBC(V, Constant((0.,0.)), left)


u = Function(V, name='Displacement')
#A,b = assemble_system(a,l,bc)
solver = KrylovSolver("cg","hypre_amg")
solver.parameters["monitor_convergence"] = True
#solver.solve(A,u.vector(),b)
solve(a == l, u, bc,solver_parameters={"linear_solver": "cg","preconditioner":"hypre_amg"})


plot(u, mode ='displacement')
plot(sigma(u)[0,0], mode='color')

# Verfication of deflection against beam-theory
w_max = -u(L,H/2.)[1]

if plane_stress: 
    Eb = E
else:
    Eb = E/(1-nu**2)
w_EB = 12*rho_g*L**4/8/Eb/H**2
w_Tim = w_EB*(1 + 2/5.*(Eb/mu)*(H/L)**2) 

print("Deflection : %g" % w_max)
print("Deflection error against Euler-Bernoulli theory : %g" % (w_max/w_EB-1))
print("Deflection error against Timoshenko theory : %g" % (w_max/w_Tim-1))
interactive()
# -*- coding: utf-8 -*-
"""
Biaxial traction of a plate with hole in 2D orthotropic elasticity

Created on Fri Jun 10 09:47:00 2016
@author: Jeremy Bleyer  (jeremy.bleyer@enpc.fr)
"""

from dolfin import *
from mshr import *

# Geometry : rectangular beam of length L and height H
L, R = 1., 0.1

# Discretization : (Nx x Ny) elements
Nx = 50
Nr = int(2*R/L*Nx)+1

# Orthotropic elasticity parameters
Ex, Ey, nux = 100, 10., 0.3
nuy = Ex*nux/Ey
mu = 20.
def C_orthotropic():
    return as_matrix([[Ex,Ex*nux,0.],[Ey*nuy,Ey,0.],[0.,0.,mu]])

# biaxial pressure
p = 1e-2

# Define epsilon and sigma operator
def eps(v):
    return sym(grad(v))
def strain2voigt(e):
    return as_vector([e[0,0],e[1,1],2*e[0,1]])
def voigt2stress(s):
    return as_tensor([[s[0],s[2]],[s[2],s[1]]])
    
def sigma(v):
    epsV = strain2voigt(eps(v))
    C = C_orthotropic()
    return voigt2stress(dot(C,epsV))

# Create mesh
domain = Rectangle(Point(0.,0.), Point(L, L)) - Circle(Point(0., 0.), R)
mesh = generate_mesh(domain, 100)

# Define function space
V = VectorFunctionSpace(mesh, 'Lagrange', 2)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
a = inner(sigma(u), eps(v))*dx


# Define right boundary x=L
class Right(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0],L) and on_boundary
class Top(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1],L) and on_boundary
class Left(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0],0) and on_boundary
class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1],0) and on_boundary
        
boundaries = MeshFunction("size_t", mesh, 1)
boundaries.set_all(0)
Right().mark(boundaries,1)
Top().mark(boundaries,2)
Left().mark(boundaries,3)
Bottom().mark(boundaries,4)
ds = Measure('ds')[boundaries]

def left(x, on_boundary):
    return near(x[0],0) and on_boundary
def bottom(x, on_boundary):
    return near(x[1],0) and on_boundary

l = p*v[0]*ds(1) + p*v[1]*ds(2)
#l = p*v[0]*dx

# symmetry boundary conditions
bc1 = DirichletBC(V.sub(0), Constant(0.), left)
bc2 = DirichletBC(V.sub(1), Constant(0.), bottom)
bc = [bc1,bc2]

u = Function(V, name='Displacement')
solve(a == l, u, bc)


plot(u, mode ='displacement')
plot(sigma(u)[0,0], mode='color')
plot(sigma(u)[1,1], mode='color')

interactive()
# -*- coding: utf-8 -*-
"""
Limit analysis solution of plate with hole in uniaxial traction
using the Linear Matching Method

Created on Fri Jun 10 12:18:00 2016
@author: Jeremy Bleyer  (jeremy.bleyer@enpc.fr)
"""

from dolfin import *
from mshr import *

# Geometry : rectangular beam of length L and height H
L, R = 1., 0.1

# Discretization : (Nx x Ny) elements
Nx = 50
Nr = int(2*R/L*Nx)+1

# Virtual elasticity parameters
E0, lmbda = 1., 1e3

# Yield stress
sig0 = 1.

# biaxial pressure
p1, p2 = Expression("t",t=1.), Constant(0.)

# Define epsilon and sigma operator
def eps(v):
    return sym(grad(v))
    
# Create mesh
domain = Rectangle(Point(0.,0.), Point(L, L)) - Circle(Point(0., 0.), R)
mesh = generate_mesh(domain, 20)

# Define function space
V = VectorFunctionSpace(mesh, 'Lagrange', 2)
V0 = FunctionSpace(mesh,'DG', 0)

E = Function(V0)
E.interpolate(Constant(E0))

def sigma(v):
    return lmbda*tr(eps(v))*Identity(2) + 2*E*eps(v)
def sig_eq(sigma):
    s = sigma - 1/2.*tr(sigma)*Identity(2)
    return sqrt(3/2.*inner(s,s))

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

#l = p1*v[0]*ds(1) + p2*v[1]*ds(2)
l = p1*v[0]*dx

# symmetry boundary conditions
bc1 = DirichletBC(V.sub(0), Constant(0.), left)
bc2 = DirichletBC(V.sub(1), Constant(0.), bottom)
bc = [bc1,bc2]

u = Function(V, name='Displacement')

for i in range(100):

    solve(a == l, u, bc)
    Sig_eq = project(sig_eq(sigma(u)),V0)
    sign = (max(Sig_eq.vector().array()) + min(Sig_eq.vector().array()))/2. 
    E.assign(project(sign/Sig_eq,V0))
    p1.t = p1.t*sig0/max(Sig_eq.vector().array())
    print("Current load mutiplier : %g" % p1.t)
    plot(u, mode ='displacement', key='displ')
    plot(Sig_eq, mode='color', key='norm_sig')

interactive()
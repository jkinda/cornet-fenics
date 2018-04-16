# -*- coding: utf-8 -*-
"""
Modal analysis of a cantilever beam in 2D plane strain elasticity

Created on Wed Jun 15 16:57:23 2016
@author: Jeremy Bleyer  (jeremy.bleyer@enpc.fr)
"""
from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
# Geometry : rectangular beam of length L and height H
L, H = 50., 1.

# Discretization : (Nx x Ny) elements
Nx = 500
Ny = int(H/L*Nx)+1

# Elasticity parameters
E, nu = 100., 0.
# Material density
rho = 1e-3

# Lamé coefficient for constitutive relation
mu = E/2./(1+nu)
lmbda = E*nu/(1+nu)/(1-2*nu)

# Define epsilon and sigma operator
def eps(v):
    return sym(grad(v))
def sigma(v):
    dim = v.geometric_dimension()
    return 2.0*mu*eps(v) + lmbda*tr(eps(v))*Identity(dim)

# Create mesh
mesh = RectangleMesh(Point(0.,0.),Point(L,H), Nx, Ny, 'crossed')

# Define function space
V = VectorFunctionSpace(mesh, 'Lagrange', degree=2)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)


def left(x, on_boundary):
    return near(x[0],0.)

bc = DirichletBC(V, Constant((0.,0.)), left)

k_form = inner(sigma(u),eps(v))*dx
K = assemble(k_form)
bc.apply(K)

m_form = rho*dot(u,v)*dx
M = assemble(m_form)


eigensolver = SLEPcEigenSolver(as_backend_type(K),as_backend_type(M))
eigensolver.parameters['problem_type'] = 'gen_hermitian'
eigensolver.parameters["spectrum"] = "smallest real"
eigensolver.parameters['spectral_transform'] = 'shift-and-invert'
eigensolver.parameters['spectral_shift'] = 0.

N_eig = 5   # number of eigenvalues
print "Computing %i first eigenvalues..." % N_eig
eigensolver.solve(N_eig)

alpha_n = (2*np.array(range(N_eig))+1)*pi/2.
# correction for the first 3 terms, correction is negligible for n>3
e_n = np.array([0.3042, -0.018 , 0.001])
alpha_n[:3] += e_n

# Extraction
for i in range(N_eig):
    # Extract eigenpair
    r, c, rx, cx = eigensolver.get_eigenpair(i)
    
    # eigenfrequency
    omega = sqrt(r)
    omega_beam = alpha_n[i]**2*sqrt(E*H**2/12./rho/L**4)
    
    print "Smallest eigenfrequency: %.5e    Beam theory: %.5e" % (omega,omega_beam)

    # Initialize function and assign eigenvector (renormalize by stiffness matrix)
    eigenmode = Function(V)
    eigenmode.vector()[:] = rx/omega
    
    # Plot eigenmode
    plot(eigenmode, mode="displacement")
    plt.show()

interactive()


#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 08:43:09 2018

@author: Jeremy Bleyer, Ecole des Ponts ParisTech,
Laboratoire Navier (ENPC,IFSTTAR,CNRS UMR 8205)
@email: jeremy.bleyer@enpc.f
"""

from dolfin import *
import numpy as np


# Generate mesh
LX, WY, WZ = 100.0, 20.0, 4.0
mesh = BoxMesh(Point(0, 0, 0), Point(LX, WY, WZ), 10, 5, 4)

# Define functions
degree=1
Ve = VectorElement("CG", mesh.ufl_cell(), 1)
Re = VectorElement("R", mesh.ufl_cell(), 0, dim=6)
Re = VectorElement("R", mesh.ufl_cell(), 0)
W = FunctionSpace(mesh, MixedElement(Ve, Ve))
w  = Function(W) # Incremental displacement
(u, lamb) = split(w)
u_, lamb_ = TestFunctions(W)
du, dlamb = TrialFunctions(W)

# Define boundary conditions
def left(x, on_boundary):
  return near(x[0], 0) and on_boundary
def right(x, on_boundary):
  return near(x[0], LX) and on_boundary
def everywhere(x, on_boundary):
  return x[0] < LX-1e-8

facets = MeshFunction("size_t", mesh, 2)
AutoSubDomain(right).mark(facets, 1)
ds = Measure("ds", subdomain_data=facets, domain=mesh)
x = SpatialCoordinate(mesh)
xG, yG, zG = (assemble(x[i]*ds(1)) for i in range(3))

# BCs
bcs = [DirichletBC(W.sub(0), Constant((0.0, 0.0, 0.0)), left),
       DirichletBC(W.sub(1), Constant((0.0,)*3), everywhere, method="pointwise")]

# Material parameters
E, nu = 2.0, 0.2
mu, lmbda = Constant(E/(2.0*(1.0 + nu))), Constant(E*nu/((1.0 + nu)*(1.0 - 2.0*nu)))

def eps(v):
    return sym(grad(v))
def sigma(v):
    dim = v.geometric_dimension()
    return 2.0*mu*eps(v) + lmbda*tr(eps(v))*Identity(dim)

U = Constant((-1, 0, 0))
a = inner(sigma(du), eps(u_))*dx + dot(lamb_, du-U)*ds(1) + dot(dlamb, u_)*ds(1)

K = PETScMatrix()
b = PETScVector()
assemble_system(lhs(a), rhs(a), bcs, A_tensor=K, b_tensor=b)
solve(K, w.vector(), b)

Vsig = TensorFunctionSpace(mesh, "CG", 1)
sig = Function(Vsig, name="Stress")
sig.assign(project(sigma(u), Vsig))

uu = w.sub(0, deepcopy=True)

file_results = XDMFFile("bucking_3D.xdmf")
file_results.parameters["flush_output"] = True
file_results.parameters["functions_share_mesh"] = True

file_results.write(uu, 0)
file_results.write(sig, 0)


kg_form = inner(grad(du),dot(sigma(u),grad(u_)))*dx

KG = PETScMatrix()
assemble(kg_form, tensor=KG)
for bci in bcs:
    bci.zero(KG)
    

eigensolver = SLEPcEigenSolver(K)
#eigensolver.parameters['problem_type'] = 'gen_hermitian'
eigensolver.parameters['solver'] = 'arnoldi'
eigensolver.parameters['spectral_transform'] = 'shift-and-invert'
eigensolver.parameters['spectral_shift'] = 1e-3
eigensolver.parameters['tolerance'] = 1e-4

N_eig = 3   # number of eigenvalues
print("Computing %i first eigenvalues..." % N_eig)
eigensolver.solve(N_eig)


#plt.figure()
# Extraction
print("Critical buckling loads:")
for i in range(N_eig):
    # Extract eigenpair
    r, c, rx, cx = eigensolver.get_eigenpair(i)
    wx = Function(W, name="Eigenvector "+str(i+1))
    wx.vector().set_local(rx.get_local())
    file_results.write(wx.sub(0, True), 0)
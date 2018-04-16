# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 15:04:26 2017

@author: bleyerj
"""

from dolfin import *
mesh = UnitSquareMesh(100,100)
U = FiniteElement("CG", mesh.ufl_cell(), 1)
Z = FunctionSpace(mesh,MixedElement([U,U,U]))

z = Function(Z)
(y, p, u) = split(z)

F = inner(grad(y) , grad(p))*dx - u*p*dx

yd = Constant(1.)

J = inner(y-yd, y-yd)*dx + Constant(1.)*inner(u,u)*dx

bc_y = DirichletBC(Z.sub(0), 0.0, "on_boundary")
bc_p = DirichletBC(Z.sub(1), 0.0, "on_boundary")
bcs = [bc_y, bc_p]

L = J + F
kkt = derivative(L, z, TestFunction(Z))

solve(kkt == 0, z, bcs)

plot(y)
plot(u)
interactive()
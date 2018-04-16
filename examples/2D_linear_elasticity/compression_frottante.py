# -*- coding: utf-8 -*-
"""
Cantilever beam/plate problem in 2D elasticity (plane strain and plane stress)
uing a traditional displacement-based formulation

Created on Wed Feb 10 10:06:12 2016
@author: Jeremy Bleyer  (jeremy.bleyer@enpc.fr)
"""

from dolfin import *

def compression_frottante(N,order):
    # Geometry : rectangular beam of length L and height H
    L, H = 1., 1.
    
    # Discretization : (Nx x Ny) elements
    Nx = N
    Ny = Nx
    
    # Elasticity parameters
    E,nu = 100, 0.3
    
    delta = 0.1
    
    # 2D model : plane stress or plane strain
    plane_stress = False
    
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
    mesh = RectangleMesh(Point(0.,0.),Point(L,H), Nx, Ny,"crossed")
    
    
    # Define function space
    V = VectorFunctionSpace(mesh, 'Lagrange', order)
    
    # Define variational problem
    u = TrialFunction(V)
    v = TestFunction(V)
    a = inner(sigma(u), eps(v))*dx
    l = inner(Constant((0,0)), v)*dx
    
    # Define left boundary x=0
    def left(x, on_boundary):
        return near(x[0],0.)
    def bottom(x, on_boundary):
        return near(x[1],0.)
    def top(x, on_boundary):
        return near(x[1],H)
    
    # Fully clamped boundary condition
    bc = [DirichletBC(V, Constant((0.,0.)), bottom),DirichletBC(V, Constant((0.,-delta)), top)]
    
    
    u = Function(V, name='Displacement')
    A,b = assemble_system(a,l,bc)
#    solver = KrylovSolver("cg","hypre_amg")
#    solver.parameters["monitor_convergence"] = True
    solver = LUSolver("mumps")
    solver.solve(A,u.vector(),b)
    #solve(a == l, u, bc) #,solver_parameters={"linear_solver": "cg","preconditioner":"hypre_amg"})
    
    Epot = assemble(0.5*inner(sigma(u),eps(u))*dx)
    print "Potential energy:",Epot
    E_a = 2*Epot/delta**2*H/L
    print "Apparent Young's modulus:",E_a
    
    nu_a = -assemble(eps(u)[0,0]*dx)/assemble(eps(u)[1,1]*dx)
    print "Apparent Poisson coefficient:",nu_a

    return E_a,nu_a,u
    
import matplotlib.pyplot as plt
E_ex = 115.5936
N_list = [10,20,50,100,200]
E_ex,_,u_ex = compression_frottante(200,3)
for n in N_list:
    E_a1,nu_a1,u1 = compression_frottante(n,1)
    E_a2,nu_a2,u2 = compression_frottante(n,2)
    plt.loglog(1./n,E_a1/E_ex-1,"ob")
    plt.loglog(1./n,E_a2/E_ex-1,"sr")
#    plt.loglog(1./n,errornorm(u_ex,u1,degree_rise=2),"ob")
#    plt.loglog(1./n,errornorm(u_ex,u2,degree_rise=1),"sr")
    
plt.loglog([1./n for n in N_list],[1e-3/n for n in N_list],"--k")
plt.loglog([1./n for n in N_list],[1e-2/n**2 for n in N_list],"-k")


plt.show()
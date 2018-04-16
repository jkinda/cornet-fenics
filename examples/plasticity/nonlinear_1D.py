from dolfin import *
from mshr import *
from time import time

#parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"]["representation"] = 'uflacs'

set_log_level(INFO)

# Problem parameters
sig0 = Constant(3**0.5)
E = Constant(1.)

mesh = UnitIntervalMesh(50)


# Create function space
V = FunctionSpace(mesh, "CG", 1)
We = FiniteElement("Quadrature", mesh.ufl_cell(), degree=0,quad_scheme='default')
W = FunctionSpace(mesh,We)
P0 = FunctionSpace(mesh,"DG",0)
    
sig = Function(W)
sig_old = Function(W)
n_elas = Function(W)
beta = Function(W)
p = Function(W)
u = Function(V)
du = Function(V)
Du = Function(V)
v = TrialFunction(V)
u_ = TestFunction(V)

# Define boundary conditions
def left(x, on_boundary):
    return on_boundary and near(x[0],0)
def right(x, on_boundary):
    return on_boundary and near(x[0],1)
def eps(u):
    return grad(u)[0]
    
markers = MeshFunction("size_t", mesh,1)
markers.set_all(0)

bc = [DirichletBC(V,0,left),DirichletBC(V,0.1,right)]
bc0 = [DirichletBC(V,0,left),DirichletBC(V,0,right)]

sig_ = Constant(0)
n_elas_ = as_vector([0,0,0,0])
Beta = Constant(0)

Heav = lambda x: (1+x/abs(x))/2.
def sig_el(eps_el):
    return E*eps_el
def sig_tang(eps_el):
    return sig_el(eps_el) - E*beta*eps_el
def proj_sig(deps,old_sig):
    sig_elas = old_sig + sig_el(deps)
    sig_eq = abs(sig_elas)
    f_elas = sig_eq - sig0
    dp = f_elas/E*Heav(f_elas)
    beta = dp/sig_eq
    deps_p = dp/sig_eq*sig
    new_sig = sig_elas-E*deps_p
    return new_sig,beta,dp
def Fext(u):
    return -Constant(0)*u*dx


a_Newton = dot(eps(v),sig_tang(eps(u_)))*dx
res = -inner(eps(u_),sig)*dx + Fext(u_)

tic = time()
niter=0
nitermax = 1e3
for i in range(10):
    A,Res = assemble_system(a_Newton,res,bc)
    nRes = Res.norm("l2")
#    nRes=1.
    Du.interpolate(Constant(0.))
    print "Increment:",str(i+1)
    while nRes > 1e-8 and niter < nitermax:
        solve(A, du.vector(), Res,"mumps")
        Du.assign(Du+du)
        deps = eps(Du)
        sig_,beta_,dp_ = proj_sig(deps,sig_old)
#        local_project(sig_,W,sig)
#        local_project(beta_,W,beta)
        sig.assign(project(sig_,W))
        beta.assign(project(beta_,W))
#        Fint_ex = inner(voigt_strain(eps(u_)),sig)*dx(metadata = {"quadrature_degree":2}) 
#        res = -Fint_ex+Fext(u_)
##        Res = assemble(res)
        A,Res = assemble_system(a_Newton,res,bc0)
        nRes = Res.norm("l2")
        print "    Residual:",nRes
        niter += 1
    u.assign(u+Du)
    sig_old.assign(sig)
    p.assign(p+project(dp_,W))
#    plot(project(p,P0),mode="color",key="p")
#    plot(project(sig[3],P0),mode="color",key="s")
#    interactive()


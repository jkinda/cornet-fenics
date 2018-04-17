from dolfin import *
from mshr import *
from time import time

parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"]["representation"] = 'uflacs'

set_log_level(INFO)

# Problem parameters
sig0 = Constant(3**0.5)
E = Constant(100.)
nu = Constant(0.3)
lmbda = E*nu/(1+nu)/(1-2*nu)
mu = E/2./(1+nu)

B = 0.5
L = 5.
H = 3.
N = 1

mesh = RectangleMesh(Point(0,0),Point(L,-H),50,30)


# Create function space
deg_vel = 2
ns = deg_vel-1
V = VectorFunctionSpace(mesh, "CG", deg_vel)
We = VectorElement("Quadrature", mesh.ufl_cell(), degree=ns,dim=4,quad_scheme='default')
W = FunctionSpace(mesh,We)
W0e = FiniteElement("Quadrature", mesh.ufl_cell(), degree=ns,quad_scheme='default')
W0 = FunctionSpace(mesh, W0e)
P0 = FunctionSpace(mesh,"DG",0)
    
sig = Function(W)
sig_old = Function(W)
n_elas = Function(W)
beta = Function(W0)
p = Function(W0)
u = Function(V)
du = Function(V)
Du = Function(V)
v = TrialFunction(V)
u_ = TestFunction(V)

# Define boundary conditions
def left(x, on_boundary):
    return on_boundary and near(x[0],0)
def bottom(x, on_boundary):
    return on_boundary and near(x[1],-H)
def right(x, on_boundary):
    return on_boundary and near(x[0],L)
class Footing(SubDomain):
    def inside(self,x, on_boundary):
        return near(x[1],0) and x[0] <= B and on_boundary
def eps(u):
    return sym(grad(u))
    
markers = MeshFunction("size_t", mesh,1)
markers.set_all(0)
Footing().mark(markers,1)
ds = Measure('ds')[markers]

bc = [DirichletBC(V.sub(1),0,bottom),DirichletBC(V.sub(0),0,left),DirichletBC(V.sub(1),-0.02,markers,1)]
bc0 = [DirichletBC(V.sub(1),0,bottom),DirichletBC(V.sub(0),0,left),DirichletBC(V.sub(1),0,markers,1)]
#bc = [DirichletBC(V,Constant((0.,0.)),left),DirichletBC(V.sub(1),Constant(-50.),right)]
#bc0 = [DirichletBC(V,Constant((0.,0.)),left),DirichletBC(V.sub(1),Constant(0.),right)]


sig_ = Constant(as_vector([0,0,0,0]))
n_elas_ = as_vector([0,0,0,0])
Beta = Constant(0)

Heav = lambda x: (1+x/abs(x))/2.
def sig3D(eps_el):
    return lmbda*tr(eps_el)*Identity(3) + 2*mu*eps_3D(eps_el)
def sig2D(eps_el):
    return lmbda*tr(eps_el)*Identity(2) + 2*mu*eps_el
def sig3D_tang(eps_el):
    N_elas = as_3D_tensor(n_elas)
    return sig3D(eps_el) - 3*mu*(1-beta)*inner(N_elas,eps_el)*N_elas-2*mu*beta*dev(eps_el)
def as_3D_tensor(X,Voigt=False):
    if Voigt:
        a = 2.
    else:
        a = 1.
    return as_tensor([[X[0],X[3]/a,0],[X[3]/a,X[1],0],[0,0,X[2]]])
def dev(a):
    return a -1/3.*tr(a)*Identity(3)
def proj_sig(deps,old_sig):
    sig_n = as_3D_tensor(old_sig)
    sig_elas = sig_n + sig3D(deps)
    s = dev(sig_elas)
    sig_eq = sqrt(3/2.*inner(s,s))
    f_elas = sig_eq - sig0
    n_elas = s/sig_eq*Heav(f_elas)
    dp = f_elas/3./mu*Heav(f_elas)
    beta = 3*mu*dp/sig_eq
    deps_p = 3*dp/2./sig_eq*s
    new_sig = sig_elas-2*mu*deps_p
    return as_vector([new_sig[0,0],new_sig[1,1],new_sig[2,2],new_sig[0,1]]),as_vector([n_elas[0,0],n_elas[1,1],n_elas[2,2],n_elas[0,1]]),beta,dp
def voigt_strain(e):
    return as_vector([e[0,0],e[1,1],0,2*e[0,1]])    
def eps_3D(e):
    return as_tensor([[e[0,0],e[0,1],0],[e[0,1],e[1,1],0],[0,0,0]])
def Fext(u):
    return -Constant(0)*u[1]*ds(1)

metadata = {"quadrature_degree":ns,"quadrature_scheme":"default"}

def local_project(v,V,u=None):
    dv = TrialFunction(V)
    v_ = TestFunction(V)
    a_proj = inner(dv,v_)*dx(metadata=metadata)
    b_proj = inner(v,v_)*dx(metadata=metadata)
    solver = LocalSolver(a_proj,b_proj)
    solver.factorize()
    if u is None:
        u = Function(V)
        solver.solve_local_rhs(u)
        return u
    else:
        solver.solve_local_rhs(u)
        return

a_Newton = inner(eps_3D(eps(v)),sig3D_tang(eps_3D(eps(u_))))*dx(metadata = metadata)
res = -inner(voigt_strain(eps(u_)),sig)*dx(metadata = metadata) + Fext(u_)

tic = time()
niter=0
nitermax = 1e3
for i in range(10):
    A,Res = assemble_system(a_Newton,res,bc)
    nRes = Res.norm("l2")
    Du.interpolate(Constant((0,0)))
    print "Increment:",str(i+1)
    while nRes > 1e-8 and niter < nitermax:
        solve(A, du.vector(), Res,"mumps")
        Du.assign(Du+du)
        deps = eps(Du)
        sig_,n_elas_,beta_,dp_ = proj_sig(deps,sig_old)
        local_project(sig_,W,sig)
        local_project(n_elas_,W,n_elas)
        local_project(beta_,W0,beta)
#        Fint_ex = inner(voigt_strain(eps(u_)),sig)*dx(metadata = {"quadrature_degree":2}) 
#        res = -Fint_ex+Fext(u_)
##        Res = assemble(res)
        A,Res = assemble_system(a_Newton,res,bc0)
        nRes = Res.norm("l2")
        print "    Residual:",nRes
        niter += 1
    u.assign(u+Du)
    sig_old.assign(sig)
    p.assign(p+local_project(dp_,W0))
#    plot(project(p,P0),mode="color",key="p")
#    plot(project(sig[3],P0),mode="color",key="s")
#    interactive()


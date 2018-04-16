from dolfin import *
set_log_level(ERROR)


mu, sig0 = 1., 10.
N = 10
R = 2*mu

ppos = lambda x: (x+abs(x))/2.

def eps(v):
    return sym(grad(v))

# Read mesh
mesh = UnitSquareMesh(N,N,"crossed")

# Create function space
V1 = VectorFunctionSpace(mesh, "CG", 2)
V0 = TensorFunctionSpace(mesh,"DG",1)
Vp = FunctionSpace(mesh,"CG",1)
V = MixedFunctionSpace([V1,Vp])
# Define functions
d = Function(V0)
s = Function(V0)

du = TrialFunction(V)
u_ = TestFunction(V)
u = Function(V)

(dv,dp) = split(du)
(v_,p_) = split(u_)
(v,p) = split(u)

a = (inner(eps(dv),eps(v_)) - dp*tr(eps(v_)) -p_*tr(eps(dv)))*dx
l = inner((d-R*s),eps(v_))*dx

def cavity_border(x, on_boundary):
    return (near(x[0],0) or near(x[1],0) or near(x[0],1.)) and on_boundary
def surface(x, on_boundary):
    TOL = 0.
    return near(x[1],1.) and (x[0] > TOL) and (x[0] < 1.-TOL) and on_boundary
    
bc_u1 = DirichletBC(V.sub(0),Constant((0.,0.)), cavity_border)
bc_u2 = DirichletBC(V.sub(0),Constant((1.,0.)), surface)
bc = [bc_u1,bc_u2]

niter = 1000
tol = 1e-6
for i in range(niter):
    
    solve(a == l, u, bc)
    sig = s + R*eps(v)
    
    nsig = sqrt(inner(sig,sig))    
    
    plot(v,title="Velocity",key="vel")
  
    d.assign(project(ppos(1-sig0/nsig)*sig/(2*mu+R),V0))
    
    residual = project(eps(v) - d,V0)
    s.assign(s + 1/R*residual)
    
    norm = assemble(inner(residual,residual)**2*dx)
    print("Residual L2 norm : %g" % norm)
    if norm < tol:
        print("Converged in %i iterations" % i)        
        break

interactive()

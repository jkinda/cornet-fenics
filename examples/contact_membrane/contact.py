"""This demo uses PETSc's TAO solver for nonlinear (bound-constrained)
optimisation problems to solve a contact mechanics problem in FEniCS.
"""
from dolfin import *
import numpy as np

set_log_level(INFO)

if not has_tao():
    print("DOLFIN must be compiled with PETSc to run this demo.")
    exit(0)
    
k = Constant(1.)

f = Expression("2")
adh = Constant(0.1)
rho = Constant(0.05)

k1 = 2
k2 = 16
#obstacle = Expression("0.1*sin(k2*pi*x[0])*sin(k2*pi*x[1])*sin(k1*pi*x[0])*sin(k1*pi*x[1])+0.25",k1=k1,k2=k2)
obstacle = Expression("0.1")

def eps(u):
    return grad(u)

# Read mesh
N = 100
mesh = RectangleMesh(Point(0.,0.),Point(1.,1.),N,N,"crossed")
#mesh = UnitIntervalMesh(N)

# Create function space
V = FunctionSpace(mesh, "Lagrange", 1)
# Define functions
du = TrialFunction(V)            # Incremental displacement
v  = TestFunction(V)             # Test function
u  = Function(V)                 # Displacement

d = Function(V,name="Obstacle")
d.interpolate(obstacle)

# Stored strain energy density (elastic membrane + adhesion)
psi = 0.5*inner(eps(u),k*eps(u)) - adh*exp(-(d-u)/rho) - f*u

# Total potential energy
Pi = psi*dx

# Compute first variation of Pi (directional derivative about u in the
# direction of v)
F = derivative(Pi, u, v)

# Compute Jacobian of F
J = derivative(F, u, du)


def border(x, on_boundary):
    return on_boundary
    
bc = DirichletBC(V, Constant(0.), border)

# Define the minimisation problem by using OptimisationProblem class
class ContactProblem(OptimisationProblem):

    def __init__(self):
        OptimisationProblem.__init__(self)

    # Objective function
    def f(self, x):
        u.vector()[:] = x
        return assemble(Pi)

    # Gradient of the objective function
    def F(self, b, x):
        u.vector()[:] = x
        assemble(F, tensor=b)

    # Hessian of the objective function
    def J(self, A, x):
        u.vector()[:] = x
        assemble(J, tensor=A)

        
u.assign(-d)
bc.apply(u.vector())
u_min = interpolate(Constant(-1e8), V)
u_max = d
bc.apply(u_min.vector())
bc.apply(u_max.vector())

# Create the PETScTAOSolver
solver = PETScTAOSolver()
info(solver.parameters,True)
# Set some parameters
solver.parameters["method"] = "gpcg"  
#solver.parameters["linear_solver"] = "nash"
solver.parameters["line_search"] = "gpcg"
solver.parameters["preconditioner"] = "ml_amg"
solver.parameters["maximum_iterations"] = 1000
solver.parameters["monitor_convergence"] = True
solver.parameters["report"] = True



# Solve the problem
solver.solve(ContactProblem(), u.vector(), u_min.vector(), u_max.vector())
plot(u,key="u1")
print np.max(abs(u.vector().array())),assemble(psi*dx)


u.interpolate(Constant(0.))
solver.solve(ContactProblem(), u.vector(), u_min.vector(), u_max.vector())
plot(u,key="u2")
print np.max(abs(u.vector().array())),assemble(psi*dx)


# Plot the current configuration

File('contact_tao_displacement.pvd') << u
print("Finished!")

#plot(adh/rho*exp(-(d-u)/rho))
interactive()

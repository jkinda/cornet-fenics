
..    # gedit: set fileencoding=utf8 :

.. _TipsAndTricks:

Efficient assignement of vector-valued functions
------------------------------------------------

Suppose that you have a vector-valued (dimension 2) function ``v`` defined on space ``VV`` and two scalar valued functions ``u0`` and ``u1`` defined on space ``V`` that you want to assign as the components of ``v``.  
Projection of ``as_vector([u0,u1])`` onto ``VV`` or interpolation of Expression is too slow for big meshes. An efficient way to perform this is to use the ``FunctionAssigner`` between space ``VV`` and spaces ``[V,V]`` as follows :

assigner = FunctionAssigner(VV, [V, V])
assigner.assign(vv, [u0, u1]) 


Efficient projection on DG or Quadrature spaces
------------------------------------------------

For projecting a Function on a DG or Quadrature space, that is a space with no coupling between elements, the projection can be performed element-wise. For this purpose, using the LocalSolver is much more faster than performing a global projection :


def local_project(v,V):
    dv = TrialFunction(V)
    v_ = TestFunction(V)
    a_proj = inner(dv,v_)*dx(metadata=metadata)
    b_proj = inner(v,v_)*dx(metadata=metadata)
    solver = LocalSolver(a_proj,b_proj)
    solver.factorize()
    u = Function(V)
    solver.solve_local_rhs(u)
    return u

Local factorizations can be cached if projection is performed many times.


Displacement-controlled for nonlinear problems
----------------------------------------------

Real Lagrange multiplier of loading

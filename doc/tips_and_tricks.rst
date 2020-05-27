
Tips and Tricks
================

Contents:

.. toctree::
   :maxdepth: 1

   demo/tips_and_tricks/mass_lumping.ipynb
   demo/tips_and_tricks/computing_reactions.ipynb

In construction...

.. _TipsTricksProjection:

------------------------------------------------
Efficient projection on DG or Quadrature spaces
------------------------------------------------


For projecting a Function on a DG or Quadrature space, that is a space with no coupling between elements, the projection can be performed element-wise. For this purpose, using the LocalSolver is much  faster than performing a global projection::

 metadata={"quadrature_degree": deg}
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

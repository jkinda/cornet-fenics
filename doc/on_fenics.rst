------------------------
A few words about FEniCS
------------------------

Contrary to most finite element softwares in the computational mechanics community, FEniCS is a multi-purpose software that allows for the discretization and numerical resolution of a wide range of problems governed by partial differential equations.

The formulation of a problem in FEniCS relies on the definition of discrete functional spaces depending on a mesh. For instance, we can define :

.. code-block:: python

	# A continuous Galerkin ("CG") scalar function space interpolated with polynomials of degree 2
	Vh = FunctionSpace(mesh,"CG",2)
	# A discontinous Galerkin ("DG") vectorial function space (the dimension is fixed by the mesh) interpolated with polynomials of degree 1
	Wh = VectorFunctionSpace(mesh,"DG",1)

The next step is to define test and trial functions which will then be used to define forms (linear, bilinear or non-linear) representing the weak variational formulation of the problem which we aim at solving.




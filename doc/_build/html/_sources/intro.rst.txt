=========================
Introduction
=========================

These numerical tours will introduce you to a wide variety of topics in Computational Continuum Mechanics using the finite element software FEniCS.

You can find instructions on how to install FEniCS on the FEniCS project website http://fenicsproject.org. In the following numerical tours, we will use the Python interface for the different FEniCS scripts.

FEniCS is also distributed along with an important number of documented or undocumented examples, some of them will be revisited in these tours but do not hesitate to look for other interesting examples.

In particular, we advise you to go through the documentation and tutorials to get started with FEniCS.

--------------------
A few words about FEniCS
--------------------

Contrary to most finite element softwares in the computational mechanics community, FEniCS is a multi-purpose software that allows for the discretization and numerical resolution of a wide range of problems governed by partial differential equations.

The formulation of a problem in FEniCS relies on the definition of discrete functional spaces depending on a mesh. For instance, we can define :

.. code-block:: python

	# A continuous Galerkin ("CG") scalar function space interpolated with polynomials of degree 2
	Vh = FunctionSpace(mesh,"CG",2)
	# A discontinous Galerkin ("DG") vectorial function space (the dimension is fixed by the mesh) interpolated with polynomials of degree 1
	Wh = VectorFunctionSpace(mesh,"DG",1)

The next step is to define test and trial functions which will then be used to define forms (linear, bilinear or non-linear) representing the weak variational formulation of the problem which we aim at solving.








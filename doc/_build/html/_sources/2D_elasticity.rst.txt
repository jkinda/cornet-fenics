=========================
 2D linear elasticity
=========================

In this first numerical tour, we will show how to compute a static solution for a 2D isotropic linear elastic medium, either in plane stress or in plane strain, using FEniCS in a tradtional displacement-based finite element formulation.

We will illustrate this on the case of a cantilever beam modeled as a 2D medium of dimensions L x H. We first define the geometrical parameters and mesh density which will be used. The mesh is generated using the RectangleMesh function and we choose a crossed configuration for the mesh structure.

We now define the material parameters which are here given in terms of a Young's modulus and a Poisson coefficient. In the following, we will need to define the constitutive relation relating \sigma as a function of \varepsilon. Let us recall that the general expression for a 3D medium of the linear elastic isotropic constitutive relation is given by :

here we will consider a 2D model either in plane strain or in plane stress. Irrespective of this choice, we will work only with a 2D displacement vector u and will subsequently define the strain operator \epsilon as follows 

.. code-block:: python

   def eps(v):
	return sym(grad(v))

which computes the 2x2 plane components of the symmetrized gradient tensor of any 2D vectorial field.
In the plane strain case, the full 3D strain tensor is defined as follows, so that the 2x2 plane part of the stress tensor is defined in the same way as for the 3D case.

In the plane stress case, an out-of-plane \epsilon_{zz}  

Hence, the 2D constitutive relation can be defined as follows, by changing only the value of the Lamé coefficient \lambda. 


.. code-block:: python

    def sigma(v):
	dim = 2
	if plane_stress:
	    lmbda = lmbda_plane_stress
	else:
	    lmbda = lmbda_plane_strain
	return lmbda*tr(eps(v))*Identity(dim) + 2.0*mu*eps(v)

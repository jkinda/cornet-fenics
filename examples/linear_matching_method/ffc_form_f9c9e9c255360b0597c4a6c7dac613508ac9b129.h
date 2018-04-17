// This code conforms with the UFC specification version 1.6.0
// and was automatically generated by FFC version 1.6.0.
// 
// This code was generated with the following parameters:
// 
//   convert_exceptions_to_warnings: False
//   cpp_optimize:                   True
//   cpp_optimize_flags:             '-O2'
//   epsilon:                        1e-14
//   error_control:                  False
//   form_postfix:                   False
//   format:                         'ufc'
//   name:                           'ffc'
//   no-evaluate_basis_derivatives:  True
//   optimize:                       False
//   precision:                      15
//   quadrature_degree:              -1
//   quadrature_rule:                'auto'
//   representation:                 'auto'
//   restrict_keyword:               ''
//   split:                          False

#ifndef __FFC_FORM_F9C9E9C255360B0597C4A6C7DAC613508AC9B129_H
#define __FFC_FORM_F9C9E9C255360B0597C4A6C7DAC613508AC9B129_H

#include <cmath>
#include <stdexcept>
#include <fstream>
#include <ufc.h>

/// This class defines the interface for a finite element.

class ffc_form_f9c9e9c255360b0597c4a6c7dac613508ac9b129_finite_element_0: public ufc::finite_element
{
public:

  /// Constructor
  ffc_form_f9c9e9c255360b0597c4a6c7dac613508ac9b129_finite_element_0() : ufc::finite_element()
  {
    // Do nothing
  }

  /// Destructor
  virtual ~ffc_form_f9c9e9c255360b0597c4a6c7dac613508ac9b129_finite_element_0()
  {
    // Do nothing
  }

  /// Return a string identifying the finite element
  virtual const char* signature() const
  {
    return "FiniteElement('Discontinuous Lagrange', Domain(Cell('triangle', 2)), 0, None)";
  }

  /// Return the cell shape
  virtual ufc::shape cell_shape() const
  {
    return ufc::triangle;
  }

  /// Return the topological dimension of the cell shape
  virtual std::size_t topological_dimension() const
  {
    return 2;
  }

  /// Return the geometric dimension of the cell shape
  virtual std::size_t geometric_dimension() const
  {
    return 2;
  }

  /// Return the dimension of the finite element function space
  virtual std::size_t space_dimension() const
  {
    return 1;
  }

  /// Return the rank of the value space
  virtual std::size_t value_rank() const
  {
    return 0;
  }

  /// Return the dimension of the value space for axis i
  virtual std::size_t value_dimension(std::size_t i) const
  {
    return 1;
  }

  /// Evaluate basis function i at given point x in cell (actual implementation)
  static void _evaluate_basis(std::size_t i,
                              double* values,
                              const double* x,
                              const double* vertex_coordinates,
                              int cell_orientation)
  {
    // Compute Jacobian
    double J[4];
    compute_jacobian_triangle_2d(J, vertex_coordinates);
    
    // Compute Jacobian inverse and determinant
    double K[4];
    double detJ;
    compute_jacobian_inverse_triangle_2d(K, detJ, J);
    
    
    // Compute constants
    
    // Get coordinates and map to the reference (FIAT) element
    
    // Reset values
    *values = 0.0;
    
    // Array of basisvalues
    double basisvalues[1] = {0.0};
    
    // Declare helper variables
    
    // Compute basisvalues
    basisvalues[0] = 1.0;
    
    // Table(s) of coefficients
    static const double coefficients0[1] = \
    {1.0};
    
    // Compute value(s)
    for (unsigned int r = 0; r < 1; r++)
    {
      *values += coefficients0[r]*basisvalues[r];
    } // end loop over 'r'
  }

  /// Evaluate basis function i at given point x in cell (non-static member function)
  virtual void evaluate_basis(std::size_t i,
                              double* values,
                              const double* x,
                              const double* vertex_coordinates,
                              int cell_orientation) const
  {
    _evaluate_basis(i, values, x, vertex_coordinates, cell_orientation);
  }

  /// Evaluate all basis functions at given point x in cell (actual implementation)
  static void _evaluate_basis_all(double* values,
                                  const double* x,
                                  const double* vertex_coordinates,
                                  int cell_orientation)
  {
    // Element is constant, calling evaluate_basis.
    _evaluate_basis(0, values, x, vertex_coordinates, cell_orientation);
  }

  /// Evaluate all basis functions at given point x in cell (non-static member function)
  virtual void evaluate_basis_all(double* values,
                                  const double* x,
                                  const double* vertex_coordinates,
                                  int cell_orientation) const
  {
    _evaluate_basis_all(values, x, vertex_coordinates, cell_orientation);
  }

  /// Evaluate order n derivatives of basis function i at given point x in cell (actual implementation)
  static void _evaluate_basis_derivatives(std::size_t i,
                                          std::size_t n,
                                          double* values,
                                          const double* x,
                                          const double* vertex_coordinates,
                                          int cell_orientation)
  {
throw std::runtime_error("// Function evaluate_basis_derivatives not generated (compiled with -fno-evaluate_basis_derivatives)");
  }

  /// Evaluate order n derivatives of basis function i at given point x in cell (non-static member function)
  virtual void evaluate_basis_derivatives(std::size_t i,
                                          std::size_t n,
                                          double* values,
                                          const double* x,
                                          const double* vertex_coordinates,
                                          int cell_orientation) const
  {
    _evaluate_basis_derivatives(i, n, values, x, vertex_coordinates, cell_orientation);
  }

  /// Evaluate order n derivatives of all basis functions at given point x in cell (actual implementation)
  static void _evaluate_basis_derivatives_all(std::size_t n,
                                              double* values,
                                              const double* x,
                                              const double* vertex_coordinates,
                                              int cell_orientation)
  {
    // Element is constant, calling evaluate_basis_derivatives.
    _evaluate_basis_derivatives(0, n, values, x, vertex_coordinates, cell_orientation);
  }

  /// Evaluate order n derivatives of all basis functions at given point x in cell (non-static member function)
  virtual void evaluate_basis_derivatives_all(std::size_t n,
                                              double* values,
                                              const double* x,
                                              const double* vertex_coordinates,
                                              int cell_orientation) const
  {
    _evaluate_basis_derivatives_all(n, values, x, vertex_coordinates, cell_orientation);
  }

  /// Evaluate linear functional for dof i on the function f
  virtual double evaluate_dof(std::size_t i,
                              const ufc::function& f,
                              const double* vertex_coordinates,
                              int cell_orientation,
                              const ufc::cell& c) const
  {
    // Declare variables for result of evaluation
    double vals[1];
    
    // Declare variable for physical coordinates
    double y[2];
    switch (i)
    {
    case 0:
      {
        y[0] = 0.333333333333333*vertex_coordinates[0] + 0.333333333333333*vertex_coordinates[2] + 0.333333333333333*vertex_coordinates[4];
      y[1] = 0.333333333333333*vertex_coordinates[1] + 0.333333333333333*vertex_coordinates[3] + 0.333333333333333*vertex_coordinates[5];
      f.evaluate(vals, y, c);
      return vals[0];
        break;
      }
    }
    
    return 0.0;
  }

  /// Evaluate linear functionals for all dofs on the function f
  virtual void evaluate_dofs(double* values,
                             const ufc::function& f,
                             const double* vertex_coordinates,
                             int cell_orientation,
                             const ufc::cell& c) const
  {
    // Declare variables for result of evaluation
    double vals[1];
    
    // Declare variable for physical coordinates
    double y[2];
    y[0] = 0.333333333333333*vertex_coordinates[0] + 0.333333333333333*vertex_coordinates[2] + 0.333333333333333*vertex_coordinates[4];
    y[1] = 0.333333333333333*vertex_coordinates[1] + 0.333333333333333*vertex_coordinates[3] + 0.333333333333333*vertex_coordinates[5];
    f.evaluate(vals, y, c);
    values[0] = vals[0];
  }

  /// Interpolate vertex values from dof values
  virtual void interpolate_vertex_values(double* vertex_values,
                                         const double* dof_values,
                                         const double* vertex_coordinates,
                                         int cell_orientation,
                                         const ufc::cell& c) const
  {
    // Evaluate function and change variables
    vertex_values[0] = dof_values[0];
    vertex_values[1] = dof_values[0];
    vertex_values[2] = dof_values[0];
  }

  /// Map coordinate xhat from reference cell to coordinate x in cell
  virtual void map_from_reference_cell(double* x,
                                       const double* xhat,
                                       const ufc::cell& c) const
  {
    throw std::runtime_error("map_from_reference_cell not yet implemented.");
  }

  /// Map from coordinate x in cell to coordinate xhat in reference cell
  virtual void map_to_reference_cell(double* xhat,
                                     const double* x,
                                     const ufc::cell& c) const
  {
    throw std::runtime_error("map_to_reference_cell not yet implemented.");
  }

  /// Return the number of sub elements (for a mixed element)
  virtual std::size_t num_sub_elements() const
  {
    return 0;
  }

  /// Create a new finite element for sub element i (for a mixed element)
  virtual ufc::finite_element* create_sub_element(std::size_t i) const
  {
    return 0;
  }

  /// Create a new class instance
  virtual ufc::finite_element* create() const
  {
    return new ffc_form_f9c9e9c255360b0597c4a6c7dac613508ac9b129_finite_element_0();
  }

};

/// This class defines the interface for a local-to-global mapping of
/// degrees of freedom (dofs).

class ffc_form_f9c9e9c255360b0597c4a6c7dac613508ac9b129_dofmap_0: public ufc::dofmap
{
public:

  /// Constructor
  ffc_form_f9c9e9c255360b0597c4a6c7dac613508ac9b129_dofmap_0() : ufc::dofmap()
  {
    // Do nothing
  }

  /// Destructor
  virtual ~ffc_form_f9c9e9c255360b0597c4a6c7dac613508ac9b129_dofmap_0()
  {
    // Do nothing
  }

  /// Return a string identifying the dofmap
  virtual const char* signature() const
  {
    return "FFC dofmap for FiniteElement('Discontinuous Lagrange', Domain(Cell('triangle', 2)), 0, None)";
  }

  /// Return true iff mesh entities of topological dimension d are needed
  virtual bool needs_mesh_entities(std::size_t d) const
  {
    switch (d)
    {
    case 0:
      {
        return false;
        break;
      }
    case 1:
      {
        return false;
        break;
      }
    case 2:
      {
        return true;
        break;
      }
    }
    
    return false;
  }

  /// Return the topological dimension of the associated cell shape
  virtual std::size_t topological_dimension() const
  {
    return 2;
  }

  /// Return the geometric dimension of the associated cell shape
  virtual std::size_t geometric_dimension() const
  {
    return 2;
  }

  /// Return the dimension of the global finite element function space
  virtual std::size_t global_dimension(const std::vector<std::size_t>&
                                       num_global_entities) const
  {
    return num_global_entities[2];
  }

  /// Return the dimension of the local finite element function space for a cell
  virtual std::size_t num_element_dofs() const
  {
    return 1;
  }

  /// Return the number of dofs on each cell facet
  virtual std::size_t num_facet_dofs() const
  {
    return 0;
  }

  /// Return the number of dofs associated with each cell entity of dimension d
  virtual std::size_t num_entity_dofs(std::size_t d) const
  {
    switch (d)
    {
    case 0:
      {
        return 0;
        break;
      }
    case 1:
      {
        return 0;
        break;
      }
    case 2:
      {
        return 1;
        break;
      }
    }
    
    return 0;
  }

  /// Tabulate the local-to-global mapping of dofs on a cell
  virtual void tabulate_dofs(std::size_t* dofs,
                             const std::vector<std::size_t>& num_global_entities,
                             const ufc::cell& c) const
  {
    dofs[0] = c.entity_indices[2][0];
  }

  /// Tabulate the local-to-local mapping from facet dofs to cell dofs
  virtual void tabulate_facet_dofs(std::size_t* dofs,
                                   std::size_t facet) const
  {
    switch (facet)
    {
    case 0:
      {
        
        break;
      }
    case 1:
      {
        
        break;
      }
    case 2:
      {
        
        break;
      }
    }
    
  }

  /// Tabulate the local-to-local mapping of dofs on entity (d, i)
  virtual void tabulate_entity_dofs(std::size_t* dofs,
                                    std::size_t d, std::size_t i) const
  {
    if (d > 2)
    {
    throw std::runtime_error("d is larger than dimension (2)");
    }
    
    switch (d)
    {
    case 0:
      {
        
        break;
      }
    case 1:
      {
        
        break;
      }
    case 2:
      {
        if (i > 0)
      {
      throw std::runtime_error("i is larger than number of entities (0)");
      }
      
      dofs[0] = 0;
        break;
      }
    }
    
  }

  /// Tabulate the coordinates of all dofs on a cell
  virtual void tabulate_coordinates(double* dof_coordinates,
                                    const double* vertex_coordinates) const
  {
    dof_coordinates[0] = 0.333333333333333*vertex_coordinates[0] + 0.333333333333333*vertex_coordinates[2] + 0.333333333333333*vertex_coordinates[4];
    dof_coordinates[1] = 0.333333333333333*vertex_coordinates[1] + 0.333333333333333*vertex_coordinates[3] + 0.333333333333333*vertex_coordinates[5];
  }

  /// Return the number of sub dofmaps (for a mixed element)
  virtual std::size_t num_sub_dofmaps() const
  {
    return 0;
  }

  /// Create a new dofmap for sub dofmap i (for a mixed element)
  virtual ufc::dofmap* create_sub_dofmap(std::size_t i) const
  {
    return 0;
  }

  /// Create a new class instance
  virtual ufc::dofmap* create() const
  {
    return new ffc_form_f9c9e9c255360b0597c4a6c7dac613508ac9b129_dofmap_0();
  }

};

/// This class defines the interface for the tabulation of the cell
/// tensor corresponding to the local contribution to a form from
/// the integral over a cell.

class ffc_form_f9c9e9c255360b0597c4a6c7dac613508ac9b129_cell_integral_0_otherwise: public ufc::cell_integral
{
public:

  /// Constructor
  ffc_form_f9c9e9c255360b0597c4a6c7dac613508ac9b129_cell_integral_0_otherwise() : ufc::cell_integral()
  {
    // Do nothing
  }

  /// Destructor
  virtual ~ffc_form_f9c9e9c255360b0597c4a6c7dac613508ac9b129_cell_integral_0_otherwise()
  {
    // Do nothing
  }

  /// Tabulate which form coefficients are used by this integral
  virtual const std::vector<bool> & enabled_coefficients() const
  {
    static const std::vector<bool> enabled({true});
    return enabled;
  }

  /// Tabulate the tensor for the contribution from a local cell
  virtual void tabulate_tensor(double*  A,
                               const double * const *  w,
                               const double*  vertex_coordinates,
                               int cell_orientation) const
  {
    // Compute Jacobian
    double J[4];
    compute_jacobian_triangle_2d(J, vertex_coordinates);
    
    // Compute Jacobian inverse and determinant
    double K[4];
    double detJ;
    compute_jacobian_inverse_triangle_2d(K, detJ, J);
    
    // Set scale factor
    const double det = std::abs(detJ);
    
    // Compute cell volume
    
    
    // Compute circumradius of triangle in 2D
    
    
    // Array of quadrature weights.
    static const double W1 = 0.5;
    // Quadrature points on the UFC reference element: (0.333333333333333, 0.333333333333333)
    
    // Values of basis functions at quadrature points.
    static const double FE0[1][1] = \
    {{1.0}};
    
    // Reset values in the element tensor.
    for (unsigned int r = 0; r < 1; r++)
    {
      A[r] = 0.0;
    } // end loop over 'r'
    
    // Compute element tensor using UFL quadrature representation
    // Optimisations: ('eliminate zeros', False), ('ignore ones', False), ('ignore zero tables', False), ('optimisation', False), ('remove zero terms', False)
    
    // Loop quadrature points for integral.
    // Number of operations to compute element tensor for following IP loop = 6
    for (unsigned int ip = 0; ip < 1; ip++)
    {
      
      // Coefficient declarations.
      double F0 = 0.0;
      
      // Total number of operations to compute function values = 2
      for (unsigned int r = 0; r < 1; r++)
      {
        F0 += FE0[0][0]*w[0][0];
      } // end loop over 'r'
      
      // Number of operations for primary indices: 4
      for (unsigned int j = 0; j < 1; j++)
      {
        // Number of operations to compute entry: 4
        A[j] += FE0[0][j]*0.506291612414141/(F0)*W1*det;
      } // end loop over 'j'
    } // end loop over 'ip'
  }

};

/// This class defines the interface for the assembly of the global
/// tensor corresponding to a form with r + n arguments, that is, a
/// mapping
///
///     a : V1 x V2 x ... Vr x W1 x W2 x ... x Wn -> R
///
/// with arguments v1, v2, ..., vr, w1, w2, ..., wn. The rank r
/// global tensor A is defined by
///
///     A = a(V1, V2, ..., Vr, w1, w2, ..., wn),
///
/// where each argument Vj represents the application to the
/// sequence of basis functions of Vj and w1, w2, ..., wn are given
/// fixed functions (coefficients).

class ffc_form_f9c9e9c255360b0597c4a6c7dac613508ac9b129_form_0: public ufc::form
{
public:

  /// Constructor
  ffc_form_f9c9e9c255360b0597c4a6c7dac613508ac9b129_form_0() : ufc::form()
  {
    // Do nothing
  }

  /// Destructor
  virtual ~ffc_form_f9c9e9c255360b0597c4a6c7dac613508ac9b129_form_0()
  {
    // Do nothing
  }

  /// Return a string identifying the form
  virtual const char* signature() const
  {
    return "2985146edeefe65ee0498fff797d4783aab99180e594ec6a295315a04c58fb3b7dca3788b8fc6155ddefa2365e88e7120ba40cba9f2fd67913026075c54ad148";
  }


  /// Return the rank of the global tensor (r)
  virtual std::size_t rank() const
  {
    return 1;
  }

  /// Return the number of coefficients (n)
  virtual std::size_t num_coefficients() const
  {
    return 1;
  }

  /// Return original coefficient position for each coefficient (0 <= i < n)
  virtual std::size_t original_coefficient_position(std::size_t i) const
  {
    static const std::vector<std::size_t> position({0});
    return position[i];
  }


  /// Create a new finite element for argument function i
  virtual ufc::finite_element* create_finite_element(std::size_t i) const
  {
    switch (i)
    {
    case 0:
      {
        return new ffc_form_f9c9e9c255360b0597c4a6c7dac613508ac9b129_finite_element_0();
        break;
      }
    case 1:
      {
        return new ffc_form_f9c9e9c255360b0597c4a6c7dac613508ac9b129_finite_element_0();
        break;
      }
    }
    
    return 0;
  }

  /// Create a new dofmap for argument function i
  virtual ufc::dofmap* create_dofmap(std::size_t i) const
  {
    switch (i)
    {
    case 0:
      {
        return new ffc_form_f9c9e9c255360b0597c4a6c7dac613508ac9b129_dofmap_0();
        break;
      }
    case 1:
      {
        return new ffc_form_f9c9e9c255360b0597c4a6c7dac613508ac9b129_dofmap_0();
        break;
      }
    }
    
    return 0;
  }


  /// Return the number of cell domains
  virtual std::size_t max_cell_subdomain_id() const
  {
    return 0;
  }

  /// Return the number of exterior facet domains
  virtual std::size_t max_exterior_facet_subdomain_id() const
  {
    return 0;
  }

  /// Return the number of interior facet domains
  virtual std::size_t max_interior_facet_subdomain_id() const
  {
    return 0;
  }

  /// Return the number of vertex domains
  virtual std::size_t max_vertex_subdomain_id() const
  {
    return 0;
  }

  /// Return the number of custom domains
  virtual std::size_t max_custom_subdomain_id() const
  {
    return 0;
  }


  /// Return whether the form has any cell integrals
  virtual bool has_cell_integrals() const
  {
    return true;
  }

  /// Return whether the form has any exterior facet integrals
  virtual bool has_exterior_facet_integrals() const
  {
    return false;
  }

  /// Return whether the form has any interior facet integrals
  virtual bool has_interior_facet_integrals() const
  {
    return false;
  }

  /// Return whether the form has any vertex integrals
  virtual bool has_vertex_integrals() const
  {
    return false;
  }

  /// Return whether the form has any custom integrals
  virtual bool has_custom_integrals() const
  {
    return false;
  }


  /// Create a new cell integral on sub domain subdomain_id
  virtual ufc::cell_integral* create_cell_integral(std::size_t subdomain_id) const
  {
    return 0;
  }

  /// Create a new exterior facet integral on sub domain subdomain_id
  virtual ufc::exterior_facet_integral* create_exterior_facet_integral(std::size_t subdomain_id) const
  {
    return 0;
  }

  /// Create a new interior facet integral on sub domain subdomain_id
  virtual ufc::interior_facet_integral* create_interior_facet_integral(std::size_t subdomain_id) const
  {
    return 0;
  }

  /// Create a new vertex integral on sub domain subdomain_id
  virtual ufc::vertex_integral* create_vertex_integral(std::size_t subdomain_id) const
  {
    return 0;
  }

  /// Create a new custom integral on sub domain subdomain_id
  virtual ufc::custom_integral* create_custom_integral(std::size_t subdomain_id) const
  {
    return 0;
  }


  /// Create a new cell integral on everywhere else
  virtual ufc::cell_integral* create_default_cell_integral() const
  {
    return new ffc_form_f9c9e9c255360b0597c4a6c7dac613508ac9b129_cell_integral_0_otherwise();
  }

  /// Create a new exterior facet integral on everywhere else
  virtual ufc::exterior_facet_integral* create_default_exterior_facet_integral() const
  {
    return 0;
  }

  /// Create a new interior facet integral on everywhere else
  virtual ufc::interior_facet_integral* create_default_interior_facet_integral() const
  {
    return 0;
  }

  /// Create a new vertex integral on everywhere else
  virtual ufc::vertex_integral* create_default_vertex_integral() const
  {
    return 0;
  }

  /// Create a new custom integral on everywhere else
  virtual ufc::custom_integral* create_default_custom_integral() const
  {
    return 0;
  }

};

#endif
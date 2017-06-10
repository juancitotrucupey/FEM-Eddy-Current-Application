/*
 *    An MIT forward solver code based on the deal.II (www.dealii.org) library.
 *    Copyright (C) 2013-2015 Ross Kynch & Paul Ledger, Swansea Unversity.
 * 
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.    
*/

#include <all_data.h>

using namespace dealii;

namespace EquationData
{
    
  // Electromagnetic constants:
  const double constant_epsilon0 = 8.85418782e-12; // electric constant (permittivity)
  const double constant_mu0 = 4.0*numbers::PI*1e-7;//1.25663706e-6; // magnetic constant (permeability)
  const double constant_sigma0 = 0.0; // background conductivity, set to regularise system (can't solve curlcurlE = f).
  
  // Material Parameters:
  double param_omega = 1.0; // angular frequency, rads/sec.
  double param_regularisation = 0.01;
  
  // Material parameters for non-conducting region:
  double param_epsilon_background = 0.0;
  double param_sigma_background = 0.0;
  double param_mur_background = 1.0;
  
  // 
  // Material parameter of the conducting object
  
  double param_epsilon_conducting = 0.0;
  double param_sigma_conducting = 0.0;
  double param_mur_conducting = 1.0;  
  
  /* Vectors holding equation parameters
   * intially assume 1 objects.
   * This can be adjusted via the parameter file and then
   * use vector.reinit(number_of objects) later on.
   */
  Vector<double> param_mur(1);
  Vector<double> param_sigma(1);
  Vector<double> param_epsilon(1);
  // Kappa = Kappa_re + i*Kappa_im = -omega.^2*epr + i*omega*sigma
  Vector<double> param_kappa_re(1); 
  Vector<double> param_kappa_im(1);
  
  // Factor for the RHS of the equation.
  double rhs_factor = 1.0;
  
  bool neumann_flag = false;
}

namespace PolarizationTensor
{
  bool enable = false;
  bool takeComplexConj = false;
  std::vector< FullMatrix<double> > polarizationTensor(2, FullMatrix<double> (3,3));
  std::vector< Vector<double> > H0(2, Vector<double> (3));
}

namespace MeshData
{
  // Flag for external mesh
  // (can be used within code to make use of deal's generator)
  bool external_mesh = false;
  
  std::string mesh_filename = "unspecified";
  
  // Set default to unspecified
  // if not changed, then don't use any of the other info.
  std::string boundary_shape = "unspecified";
  
  // For cubes/cuboids/boxes/etc:
  double xmax;
  double ymax;
  double zmax;
  double xmin;
  double ymin;
  double zmin;
  
  // For cylinders & spheres
  double height;
  double radius;
  double inner_radius;
  
  // For Toruses
  double major_radius;
  double minor_radius;
  
  Point<3> centre;
  
}

namespace IO_Data
{
  // Input/output data for use with the input/output classes
  std::string mesh_filename = "mesh.ucd";
  std::string parameter_filename = "input.prm";
  std::string output_filename = "solution";
  std::string output_filetype = "vtk";
  unsigned int n_subdivisions = 2;
}

namespace PreconditionerData
{
  bool use_direct = false;
  double solver_tolerance = 1e-6;
  double strengthen_diagonal = 0.0;
  unsigned int extra_off_diagonals = 0;
  bool right_preconditioning = true;
  
  bool constrain_gradients = true;
}


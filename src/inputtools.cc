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

#include <inputtools.h>

using namespace dealii;

namespace InputTools
{
  template <int dim>
  void read_in_mesh (const std::string mesh_name,
                     Triangulation<dim> &triangulation)
                                              
  {
    // Intended for reading in .ucd files generated by Cubit/Trelis
    // These meshes have material ids (blocks) starting at 1
    // so we take 1 away from all material_id's in the mesh.
    
    GridIn<dim> gridin;
    gridin.attach_triangulation(triangulation);
    std::ifstream mesh_file(mesh_name);
    gridin.read_ucd(mesh_file);

    // Adjust material_id to start at 0 instead of 1.
    for (typename Triangulation<dim>::active_cell_iterator
      cell = triangulation.begin_active();
    cell != triangulation.end();
    ++cell)
    {
      cell->set_material_id(cell->material_id()-1);
    }
  }
  // Template instantiation
  template void read_in_mesh<3>(const std::string,
                                Triangulation<3> &);
  
  
  // CLASS PARAMETERREADER MEMBERS
  ParameterReader::ParameterReader(ParameterHandler &paramhandler)
  :
  prm(paramhandler)
  {}
  
  ParameterReader::~ParameterReader ()
  {}
  
  void ParameterReader::declare_parameters()
  {
    // Declare subsections with:
    //      prm.enter_subsection("subsection_name")
    // declare entries within subsection with:
    //       prm.declare("entry_name", "default_value", Pattern::type, "Description")
    
    prm.enter_subsection ("Mesh Data");
    {
      prm.declare_entry("external mesh", "true",
                        Patterns::Bool(),
                        "Use an external mesh file (specified by mesh file).");      
      
      prm.declare_entry("mesh file", "mesh.ucd",
                        Patterns::Anything(),
                        "Name of the mesh file (with extension)");
      
      prm.declare_entry("recovery region id", "1",
                        Patterns::Integer(0,1000),
                        "Recovery region material ID");
      
      prm.declare_entry("background region id", "0",
                        Patterns::Integer(0,1000),
                        "Background region material ID (no recovery performed here)");
      
      // Geometry info about mesh:
      // Can pass info to the code so we can use particular shapes & geometry
      // information below:
      prm.declare_entry("boundary shape", "unspecified",
                        Patterns::Anything(),
                        "Shape of boundary (cube/cylinder_z)");
      
      // general & box shapes
      // Typically assume centred at (0,0,0), but allow it to change.
      prm.declare_entry("xmax", "0.0",
                        Patterns::Double(),
                        "Maximum value of x");
      
      prm.declare_entry("xmin", "0.0",
                        Patterns::Double(),
                        "Minimum value of x");
      
      prm.declare_entry("ymax", "0.0",
                        Patterns::Double(),
                        "Maximum value of y");
      
      prm.declare_entry("ymin", "0.0",
                        Patterns::Double(),
                        "Minimum value of y");
      
      prm.declare_entry("zmax", "0.0",
                        Patterns::Double(),
                        "Maximum value of z");
      
      prm.declare_entry("zmin", "0.0",
                        Patterns::Double(),
                        "Minimum value of z");
      
      prm.declare_entry("centre","0.0, 0.0, 0.0",
                        Patterns::List(Patterns::Double(),3,3,","),
                        "Centre of mesh (if applicable)");
      
      // Cylinder shapes:
      // Can also use zmax for cylinder
      prm.declare_entry("radius", "0.0",
                        Patterns::Double(0),
                        "Radius of mesh (if applicable)");
      
      prm.declare_entry("height", "0.0",
                        Patterns::Double(0),
                        "Height of mesh (if applicable)");
      
      // Spheres:
      prm.declare_entry("inner radius", "0.0",
                        Patterns::Double(0),
                        "Radius of inner mesh (if applicable)");

      // Torus shapes:
      prm.declare_entry("major radius", "0.0",
                        Patterns::Double(0),
                        "Major radius of a torus");
      
      prm.declare_entry("minor radius", "0.0",
                        Patterns::Double(0),
                        "Minor radius of a torus");
    }
    prm.leave_subsection ();
    
    prm.enter_subsection ("Output Parameters");
    {
      prm.declare_entry("Output filename", "solution",
                        Patterns::Anything(),
                        "Name of the output file (without extension)");
      
      prm.declare_entry("Output filetype", "vtk",
                        Patterns::Anything(),
                        "Output file extension");
      
      prm.declare_entry("n_subdivisions", "2",
                        Patterns::Double(0),
                        "Number of subdivisions in VTK file");
      
    }
    prm.leave_subsection ();
    
    prm.enter_subsection ("Material Parameters");
    {
      prm.declare_entry("omega", "1.0",
                        Patterns::Double(0),
                        "Angular frequency");
      
      prm.declare_entry("regularisation parameter", "0.01",
                        Patterns::Double(0),
                        "Regularisation Parameter");
      
      prm.declare_entry("background epsilon", "0.0",
                        Patterns::Double(0),
                        "Background permittivity");
      
      prm.declare_entry("background mur", "1.0",
                        Patterns::Double(0),
                        "Background (relative) permeability");
      
      prm.declare_entry("background sigma", "0.0",
                        Patterns::Double(0),
                        "Background conductivity");
      
      prm.declare_entry("object epsilon", "0.0",
                        Patterns::Double(0),
                        "Object permittivity");
      
      prm.declare_entry("object mur", "1.0",
                        Patterns::Double(0),
                        "Object (relative) permability");
      
      prm.declare_entry("object sigma", "0.0",
                        Patterns::Double(0),
                        "Object conductivity");
      
      prm.declare_entry("neumann flag", "false",
                        Patterns::Bool(),
                        "Boundary condition flag. True = Neumann, False = Dirichlet");
      
    }
    prm.leave_subsection ();
  
    prm.enter_subsection("Preconditioner Data");
    {
      // Disable the gmres solver:
      prm.declare_entry("use direct", "false",
                        Patterns::Bool(),
                        "Enable the sparse direct solver (disables the GMRES solver)");
      
      prm.declare_entry("solver tolerance", "1e-6",
                        Patterns::Double(),
                        "Iterative solver tolerance");
      
      // Data for the preconditioner
      prm.declare_entry("diagonal strengthening", "0",
                        Patterns::Double(0),
                        "SparseILU diagonal strengthening");
      
      prm.declare_entry("extra off diagonals", "0",
                        Patterns::Integer(0),
                        "SparseILU extra off diagonals");
      
      // Options for the solver:
      prm.declare_entry("right preconditioning", "true",
                        Patterns::Bool(),
                        "Left or right preconditioning");
      
      prm.declare_entry("constrain gradients", "true",
                        Patterns::Bool(),
                        "Flag to constrain gradients (to zero) in the non-conducting region");
    }
    prm.leave_subsection ();
    
    prm.enter_subsection("Polarization Tensor");
    {
      // Enable use
      prm.declare_entry("use polarization tensor", "false",
                        Patterns::Bool(),
                        "Enable/disable use of Polarization Tensor");
      
      prm.declare_entry("complex conjugate", "false",
                        Patterns::Bool(),
                        "Take the complex conjugate of the tensor before use");
      
      // Read in via list
      prm.declare_entry("polarization tensor real",
                        "1.0, 0.0, 0.0; 0.0, 1.0, 0.0; 0.0, 0.0, 1.0",
                        Patterns::List(Patterns::List(Patterns::Double(),3,3,","),3,3,";"),
                        "Real part of Polarization Tensor");
            
      prm.declare_entry("polarization tensor imaginary",
                        "1.0, 0.0, 0.0; 0.0, 1.0, 0.0; 0.0, 0.0, 1.0",
                        Patterns::List(Patterns::List(Patterns::Double(),3,3,","),3,3,";"),
                        "Imaginary part of Polarization Tensor");
      
      prm.declare_entry("uniform field real","0.0, 0.0, 1.0",
                        Patterns::List(Patterns::Double(),3,3,","),
                        "Real part of uniform background field");
      
      prm.declare_entry("uniform field imaginary","0.0, 0.0, 0.0",
                        Patterns::List(Patterns::Double(),3,3,","),
                        "Imaginary part of uniform background field");
    }
    prm.leave_subsection ();
   

  }
  
  void ParameterReader::get_matrix_from_list(const std::string entry,
                                             FullMatrix<double> &matrix_out,
                                             const unsigned int mRows,
                                             const unsigned int nCols)
  {
    /* Outputs a matrix read in by ParameterHandler as a List (Patterns::List(Patterns::Double())
     * - Assumes a square matrix but could be extended to handle non-square.
     * - Also assumes the separator is a comma - could be extended to handle any.
     * 
     * Requires that the matrix is entered using ; to signal a new row
     * and , to signal a new column
     */
    
    std::stringstream wholeline(prm.get(entry));
    std::string rowStr;
    std::string val;
    for (unsigned int i=0; i<mRows; ++i)
    {
      std::getline(wholeline,rowStr,';');
      std::stringstream row(rowStr);
      for (unsigned int j=0; j<nCols; ++j)
      {
        std::getline(row, val, ',');
        std::stringstream converter(val);
        converter >> matrix_out(i,j);
      }
    }
  }
  
  void ParameterReader::get_vector_from_list(const std::string entry,
                                             Vector<double> &vector_out,
                                             const unsigned int vector_length)
  {
    // Outputs a vector read in by ParameterHandler as a List (Patterns::List(Patterns::Double())
    // Note: this will resize the vector if the vector_length is not given explicitly.
    
    
    if (vector_length == 0)
    {
      const unsigned int actual_length = get_vector_length_from_list(entry, ',');
      vector_out.reinit(actual_length);
    }
    
    std::stringstream wholeline(prm.get(entry));
    std::string val;
    for (unsigned int i=0; i<vector_out.size(); ++i)
    {
      std::getline(wholeline,val,',');
      std::stringstream converter(val);
      converter >> vector_out(i);
    }
  }
  
  unsigned int ParameterReader::get_vector_length_from_list(const std::string entry,
                                                            const char delimiter)
  {
    // Returns the length of a vector read in as a list.

    std::string s = prm.get(entry);
    // Work out number of entries by the number of delimiters in the string + 1
    // e.g. if we have a vector 1, 2, 3 then there are 2 ","'s there so we have length 3.
    size_t n = std::count(s.begin(), s.end(), delimiter);
    
    // check the result is sane:
    Assert(n < UINT_MAX, ExcIndexRange (n, 0, UINT_MAX));
    
    unsigned int length = static_cast<unsigned int>(n);
    return length+1;
  }
  
  
  void ParameterReader::copy_to_equation_data()
  {
    // temp vector data:
    Vector<double> temp_vector1(3);
    
    // Input data:
    prm.enter_subsection("Mesh Data");
    
    MeshData::external_mesh = prm.get_bool("external mesh");
    
    IO_Data::mesh_filename = prm.get("mesh file");
    MeshData::mesh_filename = prm.get("mesh file");
    
    //InverseProblemData::recovery_region_id = prm.get_integer("recovery region id");
    //InverseProblemData::background_region_id = prm.get_integer("background region id");
    
    // Geometry info about mesh:
    MeshData::boundary_shape = prm.get("boundary shape");
    MeshData::xmax = prm.get_double("xmax");
    MeshData::xmin = prm.get_double("xmin");

    MeshData::ymax = prm.get_double("ymax");
    MeshData::ymin = prm.get_double("ymin");
      
    MeshData::zmax = prm.get_double("zmax");
    MeshData::zmin = prm.get_double("zmin");
      
    MeshData::radius = prm.get_double("radius");
    MeshData::height = prm.get_double("height");
    
    MeshData::inner_radius = prm.get_double("inner radius");
    
    MeshData::major_radius = prm.get_double("major radius");
    MeshData::minor_radius = prm.get_double("minor radius");
    
    get_vector_from_list("centre", temp_vector1, 3);
    for (unsigned int d=0; d<3; ++d)
    {
      MeshData::centre(d) = temp_vector1(d);
    }
    prm.leave_subsection();
    
    // Output data:
    prm.enter_subsection("Output Parameters");
    
    IO_Data::output_filename = prm.get("Output filename");
    IO_Data::output_filetype = prm.get("Output filetype");
    IO_Data::n_subdivisions = prm.get_double("n_subdivisions");
    
    prm.leave_subsection();
    
    // Material parameters:
    prm.enter_subsection("Material Parameters");
    
    EquationData::param_epsilon_background = prm.get_double("background epsilon");
    EquationData::param_mur_background = prm.get_double("background mur");
    EquationData::param_sigma_background = prm.get_double("background sigma");
    
    EquationData::param_omega = prm.get_double("omega");
    
    EquationData::param_regularisation = prm.get_double("regularisation parameter");
    
    EquationData::param_epsilon_conducting = prm.get_double("object epsilon");
    EquationData::param_mur_conducting = prm.get_double("object mur");
    EquationData::param_sigma_conducting = prm.get_double("object sigma");
    EquationData::neumann_flag = prm.get_bool("neumann flag");
    
    prm.leave_subsection();
    
    // Excitation and Sensor Coil Data:
    // Data:
    //      number of coils
    //      array centre
    //      array radius
    //      array start angle
    //      coil radius
    
    prm.enter_subsection("Preconditioner Data");
    PreconditionerData::use_direct = prm.get_bool("use direct");
    PreconditionerData::solver_tolerance = prm.get_double("solver tolerance");
    PreconditionerData::strengthen_diagonal = prm.get_double("diagonal strengthening");
    PreconditionerData::extra_off_diagonals = prm.get_integer("extra off diagonals");
    PreconditionerData::right_preconditioning = prm.get_bool("right preconditioning");
    PreconditionerData::constrain_gradients = prm.get_bool("constrain gradients");
    prm.leave_subsection();
    
    prm.enter_subsection("Polarization Tensor");
    PolarizationTensor::enable = prm.get_bool("use polarization tensor");
    PolarizationTensor::takeComplexConj = prm.get_bool("complex conjugate");
    get_matrix_from_list("polarization tensor real", PolarizationTensor::polarizationTensor[0], 3 , 3);
    get_matrix_from_list("polarization tensor imaginary", PolarizationTensor::polarizationTensor[1], 3, 3);
    // Take complex conjugate if required.
    if (PolarizationTensor::takeComplexConj)
    {
      FullMatrix<double> tempmat = PolarizationTensor::polarizationTensor[1];
      for (unsigned int d1=0;d1<3;++d1)
      {
        for (unsigned int d2=0;d2<3;++d2)
        {
          PolarizationTensor::polarizationTensor[1](d1,d2) = -tempmat(d1,d2);
        }
      }
    }
    get_vector_from_list("uniform field real", PolarizationTensor::H0[0],3);
    get_vector_from_list("uniform field imaginary", PolarizationTensor::H0[1],3);
    prm.leave_subsection();
   
  }
  
  void ParameterReader::read_parameters (const std::string parameter_file)
  {
    declare_parameters();
    prm.read_input (parameter_file);
  }
  void ParameterReader::read_and_copy_parameters (const std::string parameter_file)
  {
    declare_parameters();
    prm.read_input (parameter_file);
    copy_to_equation_data();
  }
  void ParameterReader::copy_parameters ()
  {
    copy_to_equation_data();
  }
  // template instantiation not needed (derived class) ???
  // END CLASS PARAMETERREADER
}

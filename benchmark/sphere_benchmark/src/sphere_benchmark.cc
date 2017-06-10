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

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>

#include <deal.II/grid/manifold_lib.h>

#include <deal.II/fe/fe_nedelec.h>

// Archivos que componen este solucionador
#include <all_data.h>
#include <backgroundfield.h>
#include <eddycurrentfunction.h>
#include <forwardsolver.h>
#include <inputtools.h>
#include <mydofrenumbering.h>
#include <mypreconditioner.h>
#include <myvectortools.h>
#include <outputtools.h>
#include <myfe_nedelec.h>


using namespace dealii;

//*****shereBenchmark*****************************************************
namespace sphereBenchmark
{
  //__________sphereBenchmark______________________________________-
  template <int dim>
  class sphereBenchmark
  {
  public:
    sphereBenchmark (const unsigned int poly_order,
                     const unsigned int mapping_order=2);
    ~sphereBenchmark ();
    void run(std::string input_filename,
             std::string output_filename,
             const unsigned int href);
  private:
    Triangulation<dim> tria; //Clase de dela.ii, triangulación
    FESystem<dim> fe; //clase de deal.ii, objeto de elementos finitos
    DoFHandler<dim> dof_handler; //Clase de deal.ii, administrar los grados de libertad

    const unsigned int poly_order;
    const unsigned int mapping_order;
    
    void initialise_materials();
    void process_mesh(bool neuman_flag);
  };
  //________sphereBenchmark,end____________________________
  //----Constructor---------
  template <int dim>
  sphereBenchmark<dim>::sphereBenchmark(const unsigned int poly_order,
                                        const unsigned int mapping_order)
  :
  fe (MyFE_Nedelec<dim>(poly_order), mapping_order), //MyFE_Nedelec<> es una clase definida en myfe_nedelec.h
  dof_handler (tria),
  poly_order(poly_order),
  mapping_order(mapping_order)
  {
  }
  
  //-----Destructor------
  template <int dim>
  sphereBenchmark<dim>::~sphereBenchmark ()
  {
    dof_handler.clear ();  
  }
  
  //------initialise_materials-------------------
  template <int dim>
  void sphereBenchmark<dim>::initialise_materials()
  {
    EquationData::param_mur.reinit(2);
    EquationData::param_mur(0) = EquationData::param_mur_background;       
    EquationData::param_mur(1) = EquationData::param_mur_conducting;
    
    EquationData::param_sigma.reinit(2);
    EquationData::param_sigma(0) = EquationData::param_sigma_background;
    EquationData::param_sigma(1) = EquationData::param_sigma_conducting;
    
    EquationData::param_epsilon.reinit(2);
    EquationData::param_epsilon(0) = EquationData::param_epsilon_background;
    EquationData::param_epsilon(1) = EquationData::param_epsilon_conducting;
    
    EquationData::rhs_factor = 0.0; // No source term for this problem.
    
    
    // kappa = -omega^2*epr + i*omega*sigma;
    // i.e. kappa_re = -omega^2
    //      kappa_im = omega*sigma
    EquationData::param_kappa_re.reinit(EquationData::param_mur.size());
    EquationData::param_kappa_im.reinit(EquationData::param_mur.size());
    for (unsigned int i=0;i<EquationData::param_mur.size();i++) // note mur and kappa must have same size:
    {
      // Unregularised - will be handled in forward solver.
      EquationData::param_kappa_re(i) = 0.0;
      EquationData::param_kappa_im(i) = EquationData::param_sigma(i)*EquationData::param_omega*EquationData::constant_mu0;
     
    }
  }
  
  //----------process_mesh---------
  template <int dim>
  void sphereBenchmark<dim>::process_mesh(bool neumann_flag)
  {
    // Routine to process the read in mesh
    typename Triangulation<dim>::cell_iterator cell;
    const typename Triangulation<dim>::cell_iterator endc = tria.end();
    
    // Make the interior sphere's boundary a spherical boundary:    
    // First set all manifold_ids to 0 (default).
    cell = tria.begin ();
    for (; cell!=endc; ++cell)
    {
      cell->set_all_manifold_ids(numbers::invalid_manifold_id); //numbers::invalid_manifold_id es un marcador de Deal.ii
    }
    // Now find those on the surface of the sphere.
    // Do this by looking through all cells with material_id
    // of the conductor, then finding any faces of those cells which
    // are shared with a cell with material_id of the non-conductor.
    cell = tria.begin ();
    for (; cell!=endc; ++cell)
    {
      // First if cell lies on boundary of the inner sphere, then flag as spherical (manifold 100).
      if (cell->material_id() == 1)
      {
        for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
        {
          if (cell->neighbor(face)->material_id() == 0)
          {
            cell->face(face)->set_all_manifold_ids(100);
          }
        }
      }
      // Also if cell lies on outer boundary, then flag as spherical (manifold 100).
      else if (cell->at_boundary())
      {
        for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
        {
          if (cell->face(face)->at_boundary())
          {
            cell->face(face)->set_all_manifold_ids(100);
          }
        }
      }
      //Es necesario este ultimo ciclo? Me parece que está poniendo como esférico el manifold de todas las celdas del no conductor
      if (cell->material_id() == 0)
      {
        cell-> set_all_manifold_ids(100);
      }
    }
    
    // Here we set the boundary ids for the boundary conditions
    // and may choose to mark the boundary as curved, etc.
    // Can also perform mesh refinement here.
    
    // Set boundaries to neumann (boundary_id = 10)
    if (neumann_flag)
    {
      cell = tria.begin ();
      for (; cell!=endc; ++cell)
      {
        for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
        {
          if (cell->face(face)->at_boundary())
          {
            cell->face(face)->set_all_boundary_ids (10);
          }
        }
      }
    }
    else
    {
      cell = tria.begin ();
      for (; cell!=endc; ++cell)
      {
        for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
        {
          if (cell->face(face)->at_boundary())
          {
            cell->face(face)->set_all_boundary_ids (0);
          }
        }
      }
    }
   }
  
  //---------run---------------------------
  template <int dim>
  void sphereBenchmark<dim>::run(std::string input_filename, 
                                 std::string output_filename,
                                 const unsigned int href)
  {
    
    ParameterHandler prm; //ParameterHandler clase perteneciente a deal.ii
    InputTools::ParameterReader param(prm); //Clase definida en InputTools   
    param.read_and_copy_parameters(input_filename); //Leer los parámetros del archivo input_filename
    
    const double sphere_radius=MeshData::inner_radius;
    static const SphericalManifold<dim> sph_boundary; //SphericalManifold clase perteneciente a deal.ii
    if (!MeshData::external_mesh)//Si no hay un archivo con el mallado
    {
      
      const double ball_radius = sphere_radius; //sphere_radious=inner_radious
      // Make the central shell account for 1/4 of the shell part of the mesh:
      const double middle_radius = 0.25*(MeshData::radius - MeshData::inner_radius) + MeshData::inner_radius;
      const double outer_radius = MeshData::radius;
      Triangulation<dim> inner_ball;
      Triangulation<dim> outer_ball;
      Triangulation<dim> whole_ball;
      Triangulation<dim> middle_shell;
      Triangulation<dim> outer_shell;
      
      GridGenerator::hyper_ball(whole_ball,
                                Point<dim> (0.0,0.0,0.0),
                                sphere_radius); //Clase perteneciente a deal.ii, se crea la esfera interior (sphere_radious=inner_radious)
      typename Triangulation<dim>::cell_iterator cell = whole_ball.begin();
      typename Triangulation<dim>::cell_iterator endc = whole_ball.end();
      for (;cell!=endc; ++cell)
      {
        cell->set_material_id(1); //Tipo de material_id en la forntera interna, debe ser conductor.
      }   
           
      GridGenerator::hyper_shell(middle_shell,
                                 Point<dim> (0.0,0.0,0.0),
                                 ball_radius,
                                 middle_radius,
                                 6); //Clase perteneciente a deal.ii, se crea el cascaron hasta middle_radious (ball_radious=inner_radious)
      cell = middle_shell.begin();
      endc = middle_shell.end();
      for (;cell!=endc; ++cell)
      {
        cell->set_material_id(0);//Tipo de material_id en el cascaron interno, debe ser no conductor
      }
      Triangulation<dim> middle_ball;
      GridGenerator::merge_triangulations(whole_ball,
                                          middle_shell,
                                          middle_ball); //Clase perteneciente a deal.ii, se crean la triangulación middle_ball como la unión de while_ball y middle_shell
      
      GridGenerator::hyper_shell(outer_shell,
                                 Point<dim> (0.0,0.0,0.0),
                                 middle_radius,
                                 outer_radius,
                                 6); //Clase perteneciente a deal.ii, se crea el cascarón externo
      cell = outer_shell.begin();
      endc = outer_shell.end();
      for (;cell!=endc; ++cell)
      {
        cell->set_material_id(0);//Tipo de material_id en el cascarón externo, debe ser no conductor
      }
      GridGenerator::merge_triangulations(middle_ball,
                                          outer_shell,
                                          tria); //Clase perteneciente a dela.ii, se unen las triangulaciones midle_ball y outer_ball para generar tria
      process_mesh(EquationData::neumann_flag); //Función interna que asigna condiciones de frontera 
      
      tria.set_manifold (100, sph_boundary);//Función de Triangulation de deal.ii, asigna un manifold a una parte de la triangulación(no se a qué parte) asigna un manifold esférico a todas las lartes que están marcadas con 100
      
      if (href>0)
      {
        tria.refine_global(href); //Función de Triangulation de deal.ii, refina href veces el mallado de la triangulación tria
      }
      
    }
    else //Si hay un archivo con el mallado
    {
      InputTools::read_in_mesh<dim>(IO_Data::mesh_filename,
                                    tria); //Leer el mallado del archivo externo
      
      
      
      process_mesh(EquationData::neumann_flag); //Función interna que asigna condiciones de frontera 
      // Set the marked boundary to be spherical:
      tria.set_manifold (100, sph_boundary); //Función de Triangulation de deal.ii, asigna una variedad esférica a las celdas con id 100 
      if (href>0)
      {
        tria.refine_global(href); //Función de Triangulation de deal.ii, refina href veces el mallado de la triangulación tria
      }
    }
    
    initialise_materials(); //función de la clase sphereBenchmark
    
    deallog << "Number of active cells:       "
    << tria.n_active_cells() //Función de Triangulation de dela.ii. Cantidad total de celdas activas en la triangualción tria
    << std::endl;
    
    // Now setup the forward problem:
    dof_handler.distribute_dofs (fe); //Clase perteneciente a deal.ii, .distribute.dofs(fe) indexa los grados de libertas del objeto de elementos finitos fe
    
    const MappingQ<dim> mapping(mapping_order, (mapping_order>1 ? true : false)); //Clase de deal.ii. Se utilizan polinomios de orden mapping_order en todas las celdas; si mapping_order >1 entonces se utilizan polinomios de orden mapping_order, en caso contrario se utilizan las bases de menos orden posbile según la triangulación.

//     const MappingQ1<dim> mapping;
    ForwardSolver::EddyCurrent<dim, DoFHandler<dim>> eddy(mapping,
                                                          dof_handler,
                                                          fe,
                                                          PreconditionerData::use_direct); //Clase definida por el usuario
    
    deallog << "Number of degrees of freedom: "
    << dof_handler.n_dofs() //función de DoFHandler de dela.ii. Cantidad de grados de libertad.  
    << std::endl;
    
    // assemble the matrix for the eddy current problem:
    deallog << "Assembling System Matrix...." << std::endl;
    eddy.assemble_matrices(dof_handler); //Funcion de la clase EddyCurrent definida por el ususario
    deallog << "Matrix Assembly complete. " << std::endl;
    
    // initialise the linear solver - precomputes any inverses for the preconditioner, etc:
    deallog << "Initialising Solver..." << std::endl;
    eddy.initialise_solver(); //Funcion de la clase EddyCurrent definida por el ususario
    deallog << "Solver initialisation complete. " << std::endl;
    
    // construct RHS for this field:
    Vector<double> uniform_field(dim+dim); //La clase Vector<double> pertenece a deal.ii
    for (unsigned int d=0;d<dim;++d)
    {
      uniform_field(d) = PolarizationTensor::H0[0](d);//Parte real del campo incidente
    }
    
    EddyCurrentFunction<dim>* boundary_conditions;//Se crea una función de deal.ii que calcula las varaibles relacionadas al campo incidente
    if (PolarizationTensor::enable)
    {
      boundary_conditions
      = new backgroundField::conductingObject_polarization_tensor<dim>
          (PolarizationTensor::H0,
           PolarizationTensor::polarizationTensor); //usar condiciones de frontera para un vector de polarización conocido     
    }
    else
    {
      boundary_conditions
      = new backgroundField::conductingSphere<dim> (sphere_radius, PolarizationTensor::H0); //usar condiciones de frontera para un campo magnético incidente conocido
    }
    
    // assemble rhs
    deallog << "Assembling System RHS...." << std::endl;
    eddy.assemble_rhs(dof_handler,
                      *boundary_conditions); //funcion de la clase ForwardSolver::EddyCurrent<dim, DoFHandler<dim>>
    deallog << "Matrix RHS complete. " << std::endl;
    
    deallog << "Solving... " << std::endl;
    Vector<double> solution; //La clase Vector<double> pertenece a deal.ii
    // solve system & storage in the vector of solutions:
    unsigned int n_gmres_iterations;
    eddy.solve(solution,n_gmres_iterations); //funcion de la clase ForwardSolver::EddyCurrent<dim, DoFHandler<dim>>, n_gmres_iterations de termina la cantidad de iteraciones usadas al solucionar el problema

    deallog << "Computed solution. " << std::endl;
    
    // Output error to screen:
    
    double l2err_wholeDomain;
    double l2err_conductor;
    double hcurlerr_wholeDomain;
    double hcurlerr_conductor;
    
    MyVectorTools::calcErrorMeasures(mapping,
                                     dof_handler,
                                     solution,
                                     *boundary_conditions,
                                     l2err_wholeDomain,
                                     l2err_conductor,
                                     hcurlerr_wholeDomain,
                                     hcurlerr_conductor);
    
    deallog << "L2 Errors | Whole Domain: " << l2err_wholeDomain << "   Conductor Only:" << l2err_conductor << std::endl;
    deallog << "HCurl Error | Whole Domain: " << hcurlerr_wholeDomain << "   Conductor Only:" << hcurlerr_conductor << std::endl;
    // Short version:
    std::cout << tria.n_active_cells()
    << " " << dof_handler.n_dofs()
    << " " << n_gmres_iterations
    << " " << l2err_wholeDomain
    << " " << l2err_conductor
    << " " << hcurlerr_wholeDomain
    << " " << hcurlerr_conductor
    << std::endl;
    
    {
      std::ostringstream tmp;
      tmp << output_filename;    
      OutputTools::output_to_vtk<dim, DoFHandler<dim>>(mapping,
                                                       dof_handler,
                                                       solution,
                                                       tmp.str(),
                                                       *boundary_conditions);
    }
    
    // Output the solution fields along a line to a text file:
    Point <dim> xaxis;
    Point <dim> yaxis;
    Point <dim> zaxis;
    Point <dim> diagaxis;
    if (MeshData::boundary_shape=="cube")
    { 
      xaxis = Point<dim> (MeshData::xmax, 0.0, 0.0);
      yaxis = Point<dim> (0.0, MeshData::ymax, 0.0);
      zaxis = Point<dim> (0.0, 0.0, MeshData::zmax);
      diagaxis = Point<dim> (MeshData::xmax, MeshData::ymax, MeshData::zmax);
    }
    else if (MeshData::boundary_shape=="cylinder_z")
    {
      xaxis = Point<dim> (MeshData::radius, 0.0, 0.0);
      yaxis = Point<dim> (0.0, MeshData::radius, 0.0);
      zaxis = Point<dim> (0.0, 0.0, MeshData::zmax);
      diagaxis = Point<dim> (MeshData::radius*cos(numbers::PI/4.0), MeshData::radius*sin(numbers::PI/4.0), MeshData::zmax);
    }
    else if (MeshData::boundary_shape=="sphere")
    {
      xaxis = Point<dim> (MeshData::radius, 0.0, 0.0);
      yaxis = Point<dim> (0.0, MeshData::radius, 0.0);
      zaxis = Point<dim> (0.0, 0.0, MeshData::radius);
      diagaxis = Point<dim> (MeshData::radius/sqrt(3), MeshData::radius/sqrt(3), MeshData::radius/sqrt(3));
    }
    
    {
    std::stringstream tmp;
    tmp << output_filename << "_xaxis";
    std::string xaxis_str = tmp.str();
    OutputTools::output_radial_values<dim> (mapping,
                                            dof_handler,
                                            solution,
                                            *boundary_conditions,
                                            uniform_field,
                                            xaxis,
                                            xaxis_str);
    }
    {
      std::stringstream tmp;
      tmp << output_filename << "_yaxis";
      std::string yaxis_str = tmp.str();
      OutputTools::output_radial_values<dim> (mapping,
                                              dof_handler,
                                              solution,
                                              *boundary_conditions,
                                              uniform_field,
                                              yaxis,
                                              yaxis_str);
    }
    {
      std::stringstream tmp;
      tmp << output_filename << "_zaxis";
      std::string zaxis_str = tmp.str();
      OutputTools::output_radial_values<dim> (mapping,
                                              dof_handler,
                                              solution,
                                              *boundary_conditions,
                                              uniform_field,
                                              zaxis,
                                              zaxis_str);
    }
    {
      std::stringstream tmp;
      tmp << output_filename << "_diagaxis";
      std::string diagaxis_str = tmp.str();
      OutputTools::output_radial_values<dim> (mapping,
                                              dof_handler,
                                              solution,
                                              *boundary_conditions,
                                              uniform_field,
                                              diagaxis,
                                              diagaxis_str);
    }
        
    delete boundary_conditions;
  }
}
//*****************sphereBenchmark,end********************************************************

int main (int argc, char* argv[])
{
  using namespace dealii;
  
  const int dim = 3;
  // Set default input:
  unsigned int p_order = 0;
  unsigned int href = 0;
  unsigned int mapping_order = 1;
  std::string output_filename = "sphere";
  std::string input_filename = "../input_files/sphere_benchmark.prm";
  
  // Allow for input from command line:
  if (argc > 0)
  {
    for (int i=1;i<argc;i++)
    {
      if (i+1 != argc)
      {
        std::string input = argv[i];
        if (input == "-p")
        {
          std::stringstream strValue;
          strValue << argv[i+1];
          strValue >> p_order;
        }
        if (input == "-m")
        {
          std::stringstream strValue;
          strValue << argv[i+1];
          strValue >> mapping_order;
          if (mapping_order < 1)
          {
            std::cout << "ERROR: mapping order must be > 0" << std::endl;
            return 1;
          }
        }
        if (input == "-i")
        {
          input_filename = argv[i+1];
        }
        if (input == "-o")
        {
          output_filename = argv[i+1];
        }
        if (input == "-h") // h refinement
        {
          std::stringstream strValue;
          strValue << argv[i+1];
          strValue >> href;
        }
      }
    }
  }
  
  // Only output to logfile, not console:
  deallog.depth_console(0);
  std::ostringstream deallog_filename;
  deallog_filename << output_filename << "_p" << p_order << ".deallog";
  std::ofstream deallog_file(deallog_filename.str());
  deallog.attach(deallog_file);
  
  sphereBenchmark::sphereBenchmark<dim> eddy_voltages(p_order,
                                                    mapping_order); //Clase definida en este archivo
  eddy_voltages.run(input_filename,
                    output_filename,
                    href); //Función de la clase spehreBenchmark<dim>
  
  deallog_file.close();
  return 0;
}

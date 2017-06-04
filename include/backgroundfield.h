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

#ifndef BACKGROUNDFIELD_H
#define BACKGROUNDFIELD_H

// deal.II includes:
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

// std includes:
#include <complex>

// My includes:
#include <eddycurrentfunction.h>
#include <all_data.h>

using namespace dealii;


// TODO: Make all of the functions of type EddyCurrentFunction.
//       Will need to add a few extra functions (e.g. the zero_xxx boolean functions).
//Clases derivadas de EddyCurrentFunction<dim> definida en curl_function.h, que a su vez es una clase derivada de Function<dim> de Deal.ii
namespace backgroundField
{
//_______________________________conductingSphere_______________________________________________________________________________ 
  // Field for a conducting sphere in a uniform background field.
  // includes a function to compute the perturbed field.
  template<int dim>
  class conductingSphere : public EddyCurrentFunction<dim>
  {
  public:
    conductingSphere(double sphere_radius,
                     const std::vector<Vector<double>> &uniform_field);
        
    void vector_value_list (const std::vector<Point<dim> > &points,
                            std::vector<Vector<double> >   &value_list,
                            const types::material_id &mat_id) const;
        
    void curl_value_list (const std::vector<Point<dim> > &points,
                          std::vector<Vector<double> >   &value_list,
                          const types::material_id &mat_id) const;
                          
    void perturbed_field_value_list (const std::vector< Point<dim> > &points,
                                     std::vector<Vector<double> > &value_list,
                                     const types::material_id &mat_id) const;
    
    void check_spherical_coordinates(const std::vector< Point<dim> > &points,
                                     std::vector<Vector<double> > &value_list) const;
                                     
    bool zero_vector() const { return false;}
    bool zero_curl() const { return false;}
    bool zero_perturbed() const { return false;}
                                        
  private:
    double sphere_radius;
    double constant_p;
    std::complex<double> constant_C;
    std::complex<double> constant_D;
    double constant_B_magnitude;
    const std::vector<Vector<double>> uniform_field;
  };
//_______________________________conductingObject_polarization_tensor_____________________________________________________________________
  // Field for a conducting object in a uniform background field.
  // Specifically, for an object where the polarization tensor is known.
  // This allows us to output values via the perturbed_field_value_list, 
  // using the pertubation tensor which must be read in as input.
  // Takes as input:
  // - a uniform field (3 real, 3 imag values), so [0][:] = real part.
  //                                               [1][:] = imag part.
  // - a polarization tensor [2][3][3], so [0][:][:] = real part
  //                                       [1][:][:] = imag part.
  template<int dim>
  class conductingObject_polarization_tensor : public EddyCurrentFunction<dim>
  {
  public:
    conductingObject_polarization_tensor(const std::vector<Vector<double> > &uniform_field,
                                         const std::vector<FullMatrix<double> > &polarizationTensor);
    
    void curl_value_list (const std::vector<Point<dim> > &points,
                          std::vector<Vector<double> >   &value_list,
                          const types::material_id &mat_id) const;

    void perturbed_field_value_list (const std::vector< Point<dim> > &points,
                                     std::vector<Vector<double> > &value_list,
                                     const types::material_id &mat_id) const;

    bool zero_curl() const { return false;}
    bool zero_perturbed() const { return false;}
                                        
  private:
    const std::vector< Vector<double> > uniform_field;
    const std::vector<FullMatrix<double> > polarizationTensor;
  };
//_____________________________WavePropagation______________________________________________________________________
  // Wave propagation
  // E = p*exp(i*k*x), x in R^3
  //
  // with p orthogonal to k, omega = |k| & |p|=1.
  template<int dim>
  class WavePropagation : public EddyCurrentFunction<dim> 
  {
  public:
    WavePropagation(Vector<double> &k_wave,
                    Vector<double> &p_wave);
    
    void vector_value_list (const std::vector<Point<dim> > &points,
                            std::vector<Vector<double> >   &values,
                            const types::material_id &mat_id) const;
                            

    void curl_value_list(const std::vector<Point<dim> > &points,
                         std::vector<Vector<double> > &value_list,
                         const types::material_id &mat_id) const;
                         
    bool zero_vector() const { return false;}
    bool zero_curl() const { return false;}
    
  private:  
    Vector<double> k_wave;
    Vector<double> p_wave;
  };
//_________________________________polynomialTest________________________________________________________________  
  // Function for testing only.
  // Solution is:
  // A = (x^2, y^2, z^2).
  // curlA = (0, 0, 0).
  template<int dim>
  class polynomialTest : public EddyCurrentFunction<dim>
  {
  public:
    void vector_value_list (const std::vector<Point<dim> > &points,
                            std::vector<Vector<double> >   &values,
                            const types::material_id &mat_id) const;

    void rhs_value_list (const std::vector<Point<dim> > &points,
                         std::vector<Vector<double> > &value_list,
                         const types::material_id &mat_id) const;
    
    bool zero_vector() const { return false;}
    bool zero_rhs() const { return false;}
  private:  
    Vector<double> k_wave;
    Vector<double> p_wave;
  };
//______________________curlUniformField___________________________________________________________________________-
  // Very simple uniform background
  // This returns zero in the vector_value part
  // and the uniform field specified in the curl_value part.
  template<int dim>
  class curlUniformField : public EddyCurrentFunction<dim>
  {
  public:
    curlUniformField(const std::vector< Vector<double> > &uniform_field);

    void vector_value_list (const std::vector<Point<dim> > &points,
                            std::vector<Vector<double> >   &value_list) const;

    void curl_value_list (const std::vector<Point<dim> > &points,
                          std::vector<Vector<double> >   &value_list) const;
                          
    bool zero_rhs() const { return false;}

  private:
    const std::vector< Vector<double> > uniform_field;
  };
  
}
#endif

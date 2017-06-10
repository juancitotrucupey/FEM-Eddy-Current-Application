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

#include <backgroundfield.h>

using namespace dealii;

namespace backgroundField
{
  // CONDUCTINGSPHERE
  template<int dim>
  conductingSphere<dim>::conductingSphere(double sphere_radius,
                                          const std::vector<Vector<double>> &uniform_field)
  :
  sphere_radius(sphere_radius),
  uniform_field(uniform_field)
  {
    constant_B_magnitude = sqrt(uniform_field[0].norm_sqr() + uniform_field[1].norm_sqr());
    
    // calculate the required constants:
    // TODO
    double sigma = EquationData::param_sigma_conducting;
    double mu_c = EquationData::param_mur_conducting*EquationData::constant_mu0;
    double mu_n = EquationData::param_mur_background*EquationData::constant_mu0;
    double omega = EquationData::param_omega;
    constant_p = sigma*mu_c*omega;
    std::complex<double> temp(0, constant_p);
    std::complex<double> v = sqrt(temp)*sphere_radius;
    
    // Calculate some useful terms:
    temp = sqrt(2.0/(numbers::PI*v));
    std::complex<double> besselhalf_plus = temp*sinh(v);
    std::complex<double> besselhalf_minus = temp*cosh(v);
    
    // If conducting & non-conducting region have same mu, then it simplifies:
    
    if (EquationData::param_mur_background == EquationData::param_mur_conducting)
    {
      constant_C = 3.0*pow(sqrt(sphere_radius),3)/(v*besselhalf_plus);
      constant_D = pow(sphere_radius,3)*( 3.0*v*besselhalf_minus - (3.0+v*v)*besselhalf_plus )
                   / ( v*v*besselhalf_plus );
    }
    else
    {    
       constant_C = 3.0*mu_c*v*pow(sqrt(sphere_radius),3) / ( (mu_c-mu_n)*v*besselhalf_minus + (mu_n*(1.0+v*v)-mu_c)*besselhalf_plus );
       constant_D = pow(sphere_radius,3)*( (2.0*mu_c + mu_n)*v*besselhalf_minus - (mu_n*(1.0+v*v) + 2.0*mu_c)*besselhalf_plus )
                    / ( (mu_c - mu_n)*v*besselhalf_minus + (mu_n*(1.0+v*v)-mu_c)*besselhalf_plus );
    }

  }
  
  template <int dim>
  void conductingSphere<dim>::vector_value_list (const std::vector<Point<dim> > &points,
                                             std::vector<Vector<double> > &value_list,
                                             const types::material_id &mat_id) const
  {
    Assert(value_list.size() == points.size(), ExcDimensionMismatch(value_list.size(), points.size())); //Assert funcion de Deal.ii que permite abortar la ejecución si value_list.size() == points.size() no se cumple, ExcDimensionMismatch es un objeto de Deal.ii y se crea cuando Assert aborta la ejecución
    
    for (unsigned int k=0; k<points.size(); ++k)
    {
      const Point<dim> &p = points[k];
      double r, theta, phi;
      std::complex<double> factor;
      // Convert (x,y,z) to (r,theta,phi), spherical polars.
      r = sqrt(p.square());
      theta = acos(p(2)/r);
//       theta = atan2(p(1),hypot(p(0),p(2)));
      phi = atan2(p(1),p(0)); // TODO
      
      if (r < 1e-15)
      {
        factor=0.0;
        value_list[k]=0;
      }        
      else if (r < sphere_radius)
      {
        std::complex<double> temp(0, constant_p);
        std::complex<double> v=sqrt(temp)*r;
        
        std::complex<double> bessel3halfs_plus = sqrt(2.0/(numbers::PI*v))*(cosh(v) - (1.0/v)*sinh(v));
        factor = 0.5*constant_B_magnitude*constant_C*bessel3halfs_plus*sin(theta)/sqrt(r);
      }
      else
      {
        factor = 0.5*constant_B_magnitude*(r + constant_D/(r*r))*sin(theta);
      }
      // Convert back to cartesian 
      // & split real/imaginary parts:
      value_list[k](0) = -factor.real()*sin(phi);
      value_list[k](1) = factor.real()*cos(phi);
      value_list[k](2) = 0.0;
      value_list[k](3) = -factor.imag()*sin(phi);
      value_list[k](4) = factor.imag()*cos(phi);
      value_list[k](5) = 0.0;
    }
  }
  
  template <int dim>
  void conductingSphere<dim>::curl_value_list (const std::vector<Point<dim> > &points,
                                             std::vector<Vector<double> > &value_list,
                                             const types::material_id &mat_id) const
  {
    Assert(value_list.size() == points.size(), ExcDimensionMismatch(value_list.size(), points.size()));
    
    for (unsigned int k=0; k<points.size(); ++k)
    {
      
      const Point<dim> &p = points[k];
      double r, theta, phi;
      // Convert (x,y,z) to (r,theta,phi), spherical polars.
      r = sqrt(p.square());
      theta = acos(p(2)/r);
//       theta = atan2(p(1),hypot(p(0),p(2)));
      phi = atan2(p(1),p(0)); // TODO
      
      
      
      if (r < 1e-15)
      {
        value_list[k]=0;
      }        
      else if (r - sphere_radius <= 0)
      {
        std::complex<double> temp(0, constant_p);
        const std::complex<double> sqrt_ip = sqrt(temp);
        const std::complex<double> v = sqrt_ip*r;
        
        // Calculate the required bessel function values:
        temp = sqrt(2.0/(numbers::PI*v));
        const std::complex<double> bessel3halfs_plus = temp*(cosh(v) - (1.0/v)*sinh(v));
        const std::complex<double> bessel5halfs_plus = temp*( (1.0 + 3.0/(v*v))*sinh(v) - (3.0/v)*cosh(v));
        const std::complex<double> besselhalf_plus = temp*sinh(v);
        
        std::complex<double> factor_r
        = (1.0/sqrt(r*r*r))*constant_B_magnitude*constant_C*bessel3halfs_plus*cos(theta);
        
        std::complex<double> factor_theta
        = -(constant_B_magnitude*constant_C*sin(theta))
        *( bessel3halfs_plus/sqrt(r) + sqrt(r)*sqrt_ip*( besselhalf_plus + bessel5halfs_plus ) )/(4.0*r);
          
          // Convert to cartesian:
        std::complex<double> factor_x = factor_r*sin(theta)*cos(phi) + factor_theta*cos(theta)*cos(phi);
        std::complex<double> factor_y = factor_r*sin(theta)*sin(phi) + factor_theta*cos(theta)*sin(phi);
        std::complex<double> factor_z = factor_r*cos(theta) - factor_theta*sin(theta);
        
        value_list[k](0) = factor_x.real();
        value_list[k](1) = factor_y.real();
        value_list[k](2) = factor_z.real();
        value_list[k](3) = factor_x.imag();
        value_list[k](4) = factor_y.imag();
        value_list[k](5) = factor_z.imag();
      }
      else // r > sphere_radius
      {
        std::complex<double> factor_r = constant_B_magnitude*(1.0 + constant_D/(r*r*r))*cos(theta);
        std::complex<double> factor_theta = -constant_B_magnitude*(1.0 - constant_D/(2.0*r*r*r))*sin(theta);
        // Convert to cartesian:
        std::complex<double> factor_x = factor_r*sin(theta)*cos(phi) + factor_theta*cos(theta)*cos(phi);
        std::complex<double> factor_y = factor_r*sin(theta)*sin(phi) + factor_theta*cos(theta)*sin(phi);
        std::complex<double> factor_z = factor_r*cos(theta) - factor_theta*sin(theta);
        
        value_list[k](0) = factor_x.real();
        value_list[k](1) = factor_y.real();
        value_list[k](2) = factor_z.real();
        value_list[k](3) = factor_x.imag();
        value_list[k](4) = factor_y.imag();
        value_list[k](5) = factor_z.imag();
      }
    }
  }
  template<int dim>
  void conductingSphere<dim>::perturbed_field_value_list (const std::vector< Point<dim> > &points,
                                                          std::vector<Vector<double> > &value_list,
                                                          const types::material_id &mat_id) const
  {
    // Returns the value of the perturbed field:
    // H_{p} = H - H_{0}
    //
    // In general, this is only valid outside of the object, so return 0
    // for any position within the object.
    //
       
    Assert(value_list.size() == points.size(), ExcDimensionMismatch(value_list.size(), points.size()));

    FullMatrix<double> rhat(dim);
    FullMatrix<double> D2G(dim);
    FullMatrix<double> eye(IdentityMatrix(3));
    std::vector<Vector<double>> result(2, Vector<double> (dim));
    
    Vector<double> scaled_vector(dim);
    
    for (unsigned int k=0; k<points.size(); ++k)
    {
      const Point<dim> &p = points[k];

      double r = sqrt(p.square());
      if (r - sphere_radius < 1e-10)
      {
        for (unsigned int i=0; i<dim; ++i)
        {
          value_list[k](i) = 0.0;
          value_list[k](i+dim) = 0.0;
        }
      }
      else
      {        
        double factor = 1.0/(4.0*numbers::PI*r*r*r);
        // Construct D2G
        for (unsigned int i=0; i<dim; ++i)
        {
          scaled_vector(i)=p(i)/r;
        }
        rhat.outer_product(scaled_vector,scaled_vector);
        D2G=0;
        D2G.add(3.0,rhat,-1.0,eye);
        
        D2G.vmult(result[0], uniform_field[0]);
        D2G.vmult(result[1], uniform_field[1]); // D2G is real valued so no extra terms.
        
        result[0] *= factor;
        result[1] *= factor;
        // Now multiply by 2*pi*polarization_tensor
        // NOTE: the polarization tensor is diagonal: M = constant_D*identityMatrix.
        //       and it is complex valued.
        for (unsigned int i=0; i<dim; ++i)
        {
          value_list[k](i) = 2.0*numbers::PI*(constant_D.real()*result[0](i) - constant_D.imag()*result[1](i));
          value_list[k](i+dim) = 2.0*numbers::PI*(constant_D.imag()*result[0](i) + constant_D.real()*result[1](i));
        }
      }
    }
  }
  template<int dim>
  void conductingSphere<dim>::check_spherical_coordinates(const std::vector< Point<dim> > &points,
                                                          std::vector<Vector<double> > &value_list) const
  {
    Assert(value_list.size() == points.size(), ExcDimensionMismatch(value_list.size(), points.size()));
    
    for (unsigned int k=0; k<points.size(); ++k)
    {
      const Point<dim> &p = points[k];
      double r, theta, phi;
      // Convert (x,y,z) to (r,theta,phi), spherical polars.
      r = sqrt(p.square());
//       theta = atan2(p(1),hypot(p(0),p(2)));
      theta = acos(p(2)/r);
      phi = atan2(p(1),p(0));
      value_list[k](0)=r;
      value_list[k](1)=theta;
      value_list[k](2)=phi;
    }
  }

  template class conductingSphere<3>; 
  // END CONDUCTINGSPHERE
  
  // conductingObject_polarization_tensor
  template<int dim>
  conductingObject_polarization_tensor<dim>
  ::conductingObject_polarization_tensor(const std::vector<Vector<double> > &uniform_field,
                                         const std::vector<FullMatrix<double> > &polarizationTensor)
  :
  uniform_field(uniform_field),
  polarizationTensor(polarizationTensor)  
  {
  }
  
  
  template <int dim>
  void conductingObject_polarization_tensor<dim>
  ::curl_value_list (const std::vector<Point<dim> > &points,
                     std::vector<Vector<double> > &value_list,
                     const types::material_id &mat_id) const
  {
    Assert(value_list.size() == points.size(), ExcDimensionMismatch(value_list.size(), points.size()));
    
    // No analytical solution, just set uniform far field conditions. 
       
    std::vector<Vector<double> > perturbed_value_list(value_list.size(), Vector<double> (dim+dim));
    perturbed_field_value_list(points, perturbed_value_list, mat_id);
    
    for (unsigned int k=0; k<points.size(); ++k)
    {
      for (unsigned int d=0; d<dim; ++d)
      {
        // Seems to only work when we treat the BC as n x curl E = 
                
        value_list[k](d) = perturbed_value_list[k](d) + uniform_field[0](d);
        //EquationData::constant_mu0*
        value_list[k](d+dim) = perturbed_value_list[k](d+dim) + uniform_field[1](d);
      }
    }
  }
  template<int dim>
  void conductingObject_polarization_tensor<dim>::perturbed_field_value_list (const std::vector< Point<dim> > &points,
                                                        std::vector<Vector<double> > &value_list,
                                                        const types::material_id &mat_id) const
  {
    // Returns the value of the perturbed field:
    // H_{p} = H - H_{0}
    //
        
    Assert(value_list.size() == points.size(), ExcDimensionMismatch(value_list.size(), points.size()));

    FullMatrix<double> rhat(dim);
    FullMatrix<double> D2G(dim);
    FullMatrix<double> eye(IdentityMatrix(3));
    std::vector<Vector<double>> D2Gresult(2, Vector<double> (dim));
    std::vector<Vector<double>> PTresult(2, Vector<double> (dim));
    
    Vector<double> scaled_vector(dim);
    
    for (unsigned int k=0; k<points.size(); ++k)
    {
      const Point<dim> &p = points[k];
      
      double r = sqrt(p.square());
      // Avoid singular point
      if (r < 1e-10)
      {
        value_list[k]=0.0;
      }
      else
      {
        // perturbed field = D2G*M*H0 = V2*V1...
        // where:
        // V1 = M*H0, (Note [0] -> real part, [1] -> imag part).
        // V1[0] = M[0]*H0[0] - M[1]*H0[1];
        // V1[1] = M[0]*H-[1] + M[1]*H0[0];
        //
        // V2 = D2G*V1 (D2G is real valued).
        // V2[0] = D2G*V1[0];
        // V2[1] = D2G*V1[1];
        double factor = 1.0/(4.0*numbers::PI*r*r*r);
        // Construct D2G
        for (unsigned int i=0; i<dim; ++i)
        {
          scaled_vector(i)=p(i)/r;
        }
        rhat.outer_product(scaled_vector,scaled_vector);
        
        polarizationTensor[0].vmult(PTresult[0], uniform_field[0]);
        //polarizationTensor[0].vmult(PTresult[1], uniform_field[1]);
        
        polarizationTensor[1].vmult(PTresult[1], uniform_field[0]);
        //polarizationTensor[1].vmult_add(PTresult[0], uniform_field[1]);
        
        D2G=0;
        D2G.add(3.0,rhat,-1.0,eye);
        
        D2G.vmult(D2Gresult[0], PTresult[0]);
        D2G.vmult(D2Gresult[1], PTresult[1]); // D2G is real valued so no extra terms.
        
        D2Gresult[0] *= factor;
        D2Gresult[1] *= factor;
        
        for (unsigned int i=0; i<dim; ++i)
        {
          value_list[k](i) = D2Gresult[0](i);
          value_list[k](i+dim) = D2Gresult[1](i);
        }
      }
    }
  }

  template class conductingObject_polarization_tensor<3>; 
  // END conductingObject_polarization_tensor
  
  // WAVE PROPAGATION
  template<int dim>
  WavePropagation<dim>::WavePropagation(Vector<double> &k_wave,
                                        Vector<double> &p_wave)
  :
  k_wave(k_wave),
  p_wave(p_wave)
  {
    // Note: need k orthogonal to p
    //       with |k| < 1, |p| = 1.
    // Examples:
    // horizontal wave.
    //   k_wave(0) = 0.0;
    //   k_wave(1) = 0.1;
    //   k_wave(2) = 0.0;
    //   p_wave(0) = 1.0;
    //   p_wave(1) = 0.0;
    //   p_wave(2) = 0.0;
    // diagonal wave.
    //   k_wave(0) = -0.1;
    //   k_wave(1) = -0.1;
    //   k_wave(2) = 0.2;
    //   p_wave(0) = 1./sqrt(3.);
    //   p_wave(1) = 1./sqrt(3.);
    //   p_wave(2) = 1./sqrt(3.);
  
    // Make sure input is sane:
    Assert (k_wave.size() == 3, ExcDimensionMismatch (k_wave.size(), 3));
    Assert (k_wave.size() == k_wave.size(), ExcDimensionMismatch (k_wave.size(), k_wave.size()));
    const double delta = 1e-10;
    Assert (abs(p_wave.norm_sqr()-1.0) < delta, ExcIndexRange(p_wave.norm_sqr(), 1, 1));
  } 

  template <int dim>
  void WavePropagation<dim>::vector_value_list (const std::vector<Point<dim> > &points,
                                                std::vector<Vector<double> > &value_list,
                                                const types::material_id &mat_id) const
  {
    Assert(value_list.size() == points.size(), ExcDimensionMismatch(value_list.size(), points.size()));

    double exponent;
    for (unsigned int i=0; i<points.size(); ++i)
    {
      const Point<dim> &p = points[i];
      exponent = k_wave(0)*p(0)+k_wave(1)*p(1)+k_wave(2)*p(2);
      for (unsigned int d=0; d<dim; ++d)
      {
        // Real:
        value_list[i](d) = ( p_wave(d) )*std::cos(exponent);
        // Imaginary:
        value_list[i](d+dim) = ( p_wave(d) )*std::sin(exponent);
      }
    }
  }
  template <int dim>
  void WavePropagation<dim>::curl_value_list(const std::vector<Point<dim> > &points,
                                             std::vector<Vector<double> > &value_list,
                                             const types::material_id &mat_id) const
  {
    Assert(value_list.size() == points.size(), ExcDimensionMismatch(value_list.size(), points.size()));
    
    double exponent;
    for (unsigned int i=0; i<points.size(); ++i)
    {
      const Point<dim> &p = points[i];
      exponent = k_wave(0)*p(0)+k_wave(1)*p(1)+k_wave(2)*p(2);
      // Real:
      value_list[i](0) = -( k_wave(1)*p_wave(2) - k_wave(2)*p_wave(1) )*std::sin(exponent);
      value_list[i](1) = -( k_wave(2)*p_wave(0) - k_wave(0)*p_wave(2) )*std::sin(exponent);
      value_list[i](2) = -( k_wave(0)*p_wave(1) - k_wave(1)*p_wave(0) )*std::sin(exponent);
      // Imaginary:
      value_list[i](3) =  ( k_wave(1)*p_wave(2) - k_wave(2)*p_wave(1) )*std::cos(exponent);
      value_list[i](4) =  ( k_wave(2)*p_wave(0) - k_wave(0)*p_wave(2) )*std::cos(exponent);
      value_list[i](5) =  ( k_wave(0)*p_wave(1) - k_wave(1)*p_wave(0) )*std::cos(exponent);
    }
  }
  template class WavePropagation<3>; 
  // END WAVE PROPAGATION
  
  // POLYNOMIALTEST:
  // No constructor needed.
  template <int dim>
  void polynomialTest<dim>::vector_value_list (const std::vector<Point<dim> > &points,
                                              std::vector<Vector<double> > &value_list,
                                              const types::material_id &mat_id) const
  {
    Assert(value_list.size() == points.size(), ExcDimensionMismatch(value_list.size(), points.size()));

    for (unsigned int i=0; i<points.size(); ++i)
      {
        const Point<dim> &p = points[i];

        /* quadratic: */
        value_list[i](0) = p(0)*p(0);
        value_list[i](1) = p(1)*p(1);
        value_list[i](2) = p(2)*p(2);
        value_list[i](3) = p(0)*p(0);
        value_list[i](4) = p(1)*p(1);
        value_list[i](5) = p(2)*p(2);
      }
  }
  template <int dim>
  void polynomialTest<dim>::rhs_value_list (const std::vector<Point<dim> > &points,
                                            std::vector<Vector<double> > &value_list,
                                            const types::material_id &mat_id) const
  {
    Assert(value_list.size() == points.size(), ExcDimensionMismatch(value_list.size(), points.size()));
    for (unsigned int i=0; i<points.size(); ++i)
      {
        const Point<dim> &p = points[i];

        /* quadratic: */
        value_list[i](0) = p(0)*p(0);
        value_list[i](1) = p(1)*p(1);
        value_list[i](2) = p(2)*p(2);
        value_list[i](3) = p(0)*p(0);
        value_list[i](4) = p(1)*p(1);
        value_list[i](5) = p(2)*p(2);
      }
  }
  // END POLYNOMIALTEST
  template class polynomialTest<3>;
  
  // curlUniformField
  template<int dim>
  curlUniformField<dim>::curlUniformField(const std::vector<Vector<double> > &uniform_field)
  :
  uniform_field(uniform_field)
  {
  }

  template <int dim>
  void curlUniformField<dim>::vector_value_list (const std::vector<Point<dim> > &points,
                                                 std::vector<Vector<double> > &value_list) const
  {
    Assert(value_list.size() == points.size(), ExcDimensionMismatch(value_list.size(), points.size()));

    // Not required, just set to zero.
    for (unsigned int k=0; k<points.size(); ++k)
    {
      value_list[k]=0.0;
    }
  }

  template <int dim>
  void curlUniformField<dim>::curl_value_list (const std::vector<Point<dim> > &points,
                                               std::vector<Vector<double> > &value_list) const
  {
    Assert(value_list.size() == points.size(), ExcDimensionMismatch(value_list.size(), points.size()));

    // Simply set the field to be the uniform field given in the constructor.
    for (unsigned int k=0; k<points.size(); ++k)
    {
      for (unsigned int d=0; d<dim; ++d)
      {
        value_list[k](d) = uniform_field[0](d);
        value_list[k](d+dim) = uniform_field[1](d);
      }
    }
  }
  template class curlUniformField<3>;
  // END curlUniformField
  
}

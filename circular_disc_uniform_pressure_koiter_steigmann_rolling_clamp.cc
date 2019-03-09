//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//    Version 1.0; svn revision $LastChangedRevision: 1097 $
//LIC//
//LIC// $LastChangedDate: 2015-12-17 11:53:17 +0000 (Thu, 17 Dec 2015) $
//LIC// 
//LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================
#include <fenv.h> 
//Generic routines
#include "generic.h" 

// The equations
#include "C1_large_displacement_plate_models.h"

// The mesh
#include "meshes/triangle_mesh.h"
#include <signal.h>
#include <limits>
#include "slow_fourier_transform.h"

using namespace std;
using namespace oomph;
using MathematicalConstants::Pi;

namespace TestSoln
{
//Shape of the domain
double A = 1.0;
double B = 1.0;
bool save_disk_space = false;
// The coupling of the stretching energy
double eta = 0;
double p_mag = 0.; 
double p_prev= 0.; 
double nu = 0.5;
double h = (0.34/80.);
double h_prev= h; 
const unsigned Number_of_fourier_terms = 10;
const unsigned First_fourier_wavenumber = 11;

// Max times to half the pressure step
unsigned Max_number_pstep_halves = 2;
unsigned Max_number_hstep_halves = 2;

/*                     PARAMETRIC BOUNDARY DEFINITIONS                        */
// Here we create the geom objects for the Parametric Boundary Definition 
CurvilineCircleTop parametric_curve_top;
CurvilineCircleBottom parametric_curve_bottom;

// Bool for doing the eigen problem
bool do_eigen=true;

// Number of eigenvalues
unsigned neig = 30;
bool use_direct_solver = false;

// Max depth for eigenvalue searching
int max_depth = 16;

// bool to turn of do eigen after we find the onset
bool only_first_singular_value = false;

// Do we need to dump at every step
bool dump_at_every_step = false;

// Get s from x
double get_s_0(const Vector<double>& x)
{
// The arc length (parametric parameter) for the upper semi circular arc
 return atan2(-x[0],x[1]);
}

// Get s from x
double get_s_1(const Vector<double>& x)
{
// The arc length (parametric parameter) for the lower semi circular arc
return atan2(x[0],-x[1]);
}

// Assigns the value of pressure depending on the position (x,y)
void get_pressure(const Vector<double>& xi,const Vector<double>& ui,
 const DenseMatrix<double>& dui_dxj,  const Vector<double>& ni, 
 Vector<double>& pressure)
{
 // We are distributing the pressure over the deformed surface
 // so the area element will be dA = \sqrt(det(g)) da with da the area element on
 // the undeformed sheet. So a 'following' pressure will be p n_i dA
 // or \sqrt(g) n_i da = r,x * r,y p da  with  * the cross product
 Vector<double> non_unit_ni(3);
 non_unit_ni[0] = (dui_dxj(1,0))*dui_dxj(2,1) - dui_dxj(2,0)*(1.0+dui_dxj(1,1));
 non_unit_ni[1] =-(1.0+dui_dxj(0,0))*dui_dxj(2,1) + dui_dxj(2,0)*(dui_dxj(0,1));
 non_unit_ni[2] = (1.0+dui_dxj(0,0))*(1.0+dui_dxj(1,1)) - dui_dxj(0,1)*(dui_dxj(1,0));
  
 for(unsigned i=0; i<3;++i)
  {
   // N.B Non dimensional  p = p* L / E h where p* is dimensional
   pressure[i] = /*h**/p_mag*non_unit_ni[i]; 
  }
}

// Assigns the value of pressure depending on the position (x,y)
inline void get_d_pressure_dn(const Vector<double>& xi,const Vector<double>& ui,
 const DenseMatrix<double>& dui_dxj, const Vector<double>& ni, 
 DenseMatrix<double>& d_pressure_dn)
{/* EMPTY */}

// Assigns the value of pressure depending on the position (x,y)
void get_d_pressure_dr(const Vector<double>& xi,const Vector<double>& ui,
  const DenseMatrix<double>& dui_dxj,
 const Vector<double>& ni, DenseMatrix<double>& d_pressure_dn)
{/* EMPTY */}

// Assigns the value of pressure depending on the position (x,y)
inline void get_d_pressure_d_grad_u(const Vector<double>& xi,const Vector<double>& ui,
 const DenseMatrix<double>& dui_dxj,
 const Vector<double>& ni, RankThreeTensor<double>& d_pressure_du_grad)
{
  // This way doesn't need an intermediate variable
  d_pressure_du_grad(0,1,0) = p_mag*dui_dxj(2,1);
  d_pressure_du_grad(0,2,1) = p_mag*dui_dxj(1,0);
  d_pressure_du_grad(0,2,0) =-p_mag*(1.0+dui_dxj(1,1));
  d_pressure_du_grad(0,1,1) =-p_mag*dui_dxj(2,0);

  d_pressure_du_grad(1,0,0) =-p_mag*dui_dxj(2,1);
  d_pressure_du_grad(1,2,1) =-p_mag*(1.0+dui_dxj(0,0));
  d_pressure_du_grad(1,2,0) = p_mag*dui_dxj(0,1);
  d_pressure_du_grad(1,0,1) = p_mag*dui_dxj(2,0);

  d_pressure_du_grad(2,0,0) = p_mag*(1.0+dui_dxj(1,1));
  d_pressure_du_grad(2,1,1) = p_mag*(1.0+dui_dxj(0,0));
  d_pressure_du_grad(2,0,1) =-p_mag*dui_dxj(1,0);
  d_pressure_du_grad(2,1,0) =-p_mag*dui_dxj(0,1);
}

// The normal and tangential directions.
void get_normal_and_tangent(const Vector<double>& x, Vector<double>& n, 
 Vector<double>& t, DenseMatrix<double>& Dn, DenseMatrix<double>& Dt)
{
 // Fill in the normal and derivatives of the normal
 n[0] = x[0]/sqrt(x[0]*x[0]+x[1]*x[1]);
 n[1] = x[1]/sqrt(x[0]*x[0]+x[1]*x[1]);

 // The derivatives of the x and y components
 Dn(0,0) = x[1]*x[1] * pow(x[0]*x[0]+x[1]*x[1],-1.5);
 Dn(1,0) =-x[1]*x[0] * pow(x[0]*x[0]+x[1]*x[1],-1.5);
 Dn(0,1) =-x[0]*x[1] * pow(x[0]*x[0]+x[1]*x[1],-1.5);
 Dn(1,1) = x[0]*x[0] * pow(x[0]*x[0]+x[1]*x[1],-1.5);

  // Fill in the tangent and derivatives of the tangent
 t[0] =-x[1]/(x[0]*x[0]+x[1]*x[1]);
 t[1] = x[0]/(x[0]*x[0]+x[1]*x[1]);

 // The derivatives of the x and y components
 Dt(0,0)=(+2*x[0]*x[1])/pow(x[0]*x[0] + x[1]*x[1],2);
 Dt(0,1)=(-x[0]*x[0] + x[1]*x[1])/pow(x[0]*x[0] + x[1]*x[1],2);
 Dt(1,0)=(-x[0]*x[0] + x[1]*x[1])/pow(x[0]*x[0] + x[1]*x[1],2);
 Dt(1,1)=(-2*x[0]*x[1])/pow(x[0]*x[0] + x[1]*x[1],2);
}

//Exact solution for constant pressure, circular domain and resting boundary conditions
void get_exact_w(const Vector<double>& xi, Vector<double>& w)
{}

void error_metric(const Vector<double>& x, const 
  Vector<double>& u, const Vector<double>& u_exact, double& error, double& norm)
{
 using namespace SlowFourierTransform;
 // We use the theta derivative of the out of plane deflection
 error = pow((-x[1]*u[13] + x[0]*u[14]),2)/(x[0]*x[0]+x[1]*x[1]);
 norm  = pow(x[0]*u[13] +x[1]*u[14],2)/(x[0]*x[0]+x[1]*x[1]);
}

void fourier_transform_metric(const Vector<double>& x, const 
  Vector<double>& u, const Vector<double>& u_exact, Vector<double>& error,
  Vector<double> & norm)
{
 using namespace SlowFourierTransform;
 // We use the theta derivative of the out of plane deflection
 double r2=x[0]*x[0]+x[1]*x[1],r=sqrt(x[0]*x[0]+x[1]*x[1]);
 for (unsigned i=0, k; i<Number_of_fourier_terms;++i)
 {
  k=First_fourier_wavenumber+i;
  // Compute the foutier transform for wavenumber k
  error[i] = u[12]*Slow_fourier_sine_function(x[0],x[1],k);
  // Compute the norm for r^k, this way we only compute sqrt once 
  norm[i]  = u[12]*pow(r2,k/2)*(k%2 == 0 ? 1.0 : r);
 }
}
  

//Exact solution for constant pressure, circular domain and resting boundary conditions
void get_exact_w_radial(const Vector<double>& xi, Vector<double>& w)
{}

// Use the JFM nondim.
bool Use_jfm_nondimensionalisation = false;

double Singular_tolerance = 1e-14;
double Singular_prefactor = .1;

double Arnoldi_multiple = 2; 

// Use a fd jacobian
bool use_fd_jacobian = true;

// Dump the triangulateio object
bool Dump_triangulateio = false;
}

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////


//==start_of_problem_class============================================
/// Class definition
//====================================================================
template<class ELEMENT>
class UnstructuredFvKProblem : public virtual Problem
{

public:

/// Constructor
UnstructuredFvKProblem(double element_area = 0.1);

/// Destructor
~UnstructuredFvKProblem()
{
 delete (Surface_mesh_pt);
 delete (Bulk_mesh_pt);
 // Clean up mesh parameters
 delete Outer_boundary_pt;
 delete Outer_boundary_ellipse_pt;
 delete Outer_curvilinear_boundary_pt[0];
 delete Outer_curvilinear_boundary_pt[1];
 delete Inner_open_boundaries_pt[0];
 // Check it's not 0
 if(Boundary2_pt!=0)
  delete Boundary2_pt;
};

void actions_after_read_unstructured_meshes()
 {
 // Curved Edge upgrade
 upgrade_edge_elements_to_curve(0,Bulk_mesh_pt);
 upgrade_edge_elements_to_curve(1,Bulk_mesh_pt);
  
 // Rotate degrees of freedom
 rotate_edge_degrees_of_freedom(Bulk_mesh_pt);
  complete_problem_setup();
  apply_boundary_conditions();
 } 
/// Update after solve (empty)
void actions_after_newton_solve()
{
}

/// Update the problem specs before solve: Re-apply boundary conditions
/// Empty as the boundary conditions stay fixed
void actions_before_newton_solve()
{
//apply_boundary_conditions();
 Is_wrinkled = false;
}

/// \short Plot error when compared against a given exact solution.
///  Also returns the norm  of the error and that of the exact solution
virtual void compute_error(std::ostream &outfile,
                           FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                           Vector<double>& error, Vector<double>& norm)
 {
 //Initialise the norm and error
 for(unsigned i=0; i<norm.size();++i)
  {norm[i] =0.0;} 
 for(unsigned i=0; i<error.size();++i)
  {error[i] =0.0;} 

 //Per-element norm and error
 Vector<double> el_error,el_norm;

 //Loop over the elements
 unsigned n_element = Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
 {
  // Try to cast to FiniteElement
  KoiterSteigmannPlateEquations<2,2>*
el_pt=dynamic_cast<KoiterSteigmannPlateEquations<2,2>*>(Bulk_mesh_pt->element_pt(e));
  if (el_pt==0)
   {
    throw OomphLibError(
     "Can't execute compute_error(...) with multiple errors for non \
KoiterSteigmannPlateElements",
     OOMPH_CURRENT_FUNCTION,
     OOMPH_EXCEPTION_LOCATION);
   }

  // Reset elemental errors and norms
  for(unsigned i=0; i<norm.size();++i)
   {norm[i] =0.0;} 
  for(unsigned i=0; i<error.size();++i)
   {error[i] =0.0;} 
  el_pt->compute_error(outfile,exact_soln_pt,el_error,el_norm);
  //Add each element error to the global errors
  for(unsigned i=0; i<norm.size();++i)
   {norm[i] += el_error[i];} 
  for(unsigned i=0; i<error.size();++i)
   {error[i] += el_norm[i]; } 
  }
 }

/// Doc the solution
void doc_solution(const std::string& comment="");

void pin_first_point()
 {
  // Get first point
  Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(0,0);
    {
     oomph_info <<"Using first node. Pinning inplane displacements and setting to zero."<<std::endl;
     nod_pt -> pin(0);
     nod_pt -> pin(6);
     nod_pt -> set_value(0,0.0);
     nod_pt -> set_value(6,0.0);
    }
 }
void pin_centre_point()
 {
  // Loop over internal boundary
  unsigned num_nod=Bulk_mesh_pt->nboundary_node(Inner_boundary2);
  for (unsigned inod=0;inod<num_nod;inod++)
  {
   // Get nod
   Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(Inner_boundary2,inod);
   Vector<double> x(2,0.0);
   x[0]=nod_pt->x(0);
   x[1]=nod_pt->x(1);
   // If its in the centre
   const double tol = 1e-20;
   if(x[0]*x[0]+x[1]*x[1]<tol)
    {
     oomph_info <<"Found centre: "<<x<<std::endl;
     oomph_info <<"Pinning inplane displacements and setting to zero."<<std::endl;
     nod_pt -> pin(0);
     nod_pt -> pin(6);
     nod_pt -> set_value(0,0.0);
     nod_pt -> set_value(6,0.0);
    }
  }
 }
void output_centre_point()
 {
  // Loop over internal boundary
  unsigned num_nod=Bulk_mesh_pt->nboundary_node(2);
  for (unsigned inod=0;inod<num_nod;inod++)
  {
   // Get nod
   Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(2,inod);
   Vector<double> x(2,0.0);
   x[0]=nod_pt->x(0);
   x[1]=nod_pt->x(1);
   // If its in the centre
   if(x[0]*x[0]+x[1]*x[1]<1e-12)
    {
     oomph_info <<"Found centre: "<<x<<std::endl;
     oomph_info <<"Centre deflection: "<<nod_pt->raw_value(12)<<std::endl;
    }
  }
 }
/// \short Overloaded version of the problem's access function to
/// the mesh. Recasts the pointer to the base Mesh object to
/// the actual mesh type.
TriangleMesh<ELEMENT>* mesh_pt()
{
return dynamic_cast<TriangleMesh<ELEMENT>*> (Problem::mesh_pt()); 
}

/// Doc info object for labeling output
DocInfo Doc_info;

// Return status
bool is_wrinkled() const {return  Is_wrinkled;};

// Recursion to find the singular value
void find_singular_p_value(const unsigned nneg_initial,
 const std::pair<unsigned,unsigned>& nneg_bracket,
 const std::pair<double,double>& p_bracket, const double& thresh,
  bool& found_singular, int& depth)
 {
  // If we reach the bottom of the recursion or have found it return
  if(depth <= 0 || found_singular)
   { 
    oomph_info<<"#################################################"<<std::endl;
    return; 
   } 
  // Oomph info
  oomph_info<<"#################################################"<<std::endl;
  oomph_info <<"Starting at layer " << depth << " to find singular value."<<std::endl;
   // New pressure halfway between
  TestSoln::p_mag = (p_bracket.second+p_bracket.first)/2.0;
  oomph_info <<"New Pressure is "<<TestSoln::p_mag<<std::endl;
  // Solve the problem
  newton_solve();
  // Doc the solution
  doc_solution();
  oomph_info << std::endl;
  oomph_info << "---------------------------------------------" << std::endl;
  oomph_info << "Pressure (" << TestSoln::p_mag << ")" << std::endl;
  oomph_info << "Current h  (" << TestSoln::h << ")" << std::endl;
  oomph_info << "Poisson ratio (" << TestSoln::nu << ")" << std::endl;
  oomph_info << "Solution number (" <<Doc_info.number()-1 << ")" << std::endl;
  oomph_info << "---------------------------------------------" << std::endl;
  oomph_info << std::endl;
  // Solve the eigenproblem
  unsigned nneg;
  solve_eigen_problem(nneg, found_singular, thresh); 
  // Decrease the depth count
  --depth;

  // Where to next
  std::pair<double,double>  p_bracket_next;
  std::pair<unsigned,unsigned> nneg_bracket_next;
  // If the bracket between test and first is contains at least one crossing
  // OR there is more than one crossing between then always go smaller
  if(nneg != nneg_bracket.first || (nneg > nneg_bracket.first + 1
     || nneg + 1 < nneg_bracket.first) )
   {
    p_bracket_next.first = p_bracket.first;
    nneg_bracket_next.first = nneg_bracket.first;
    p_bracket_next.second = TestSoln::p_mag;
    nneg_bracket_next.second = nneg;
   } 
  else if (nneg == nneg_bracket.first && nneg != nneg_bracket.second)
   {
    p_bracket_next.first = TestSoln::p_mag;
    nneg_bracket_next.first = nneg;
    p_bracket_next.second = p_bracket.second;
    nneg_bracket_next.second = nneg_bracket.second;
   }
  // This shouldnt ever happen
  else
   {
    // Escape
    oomph_info << "More than one eigenvalue has changed sign,exiting."<<std::endl;
    oomph_info << nneg_bracket.first <<","<< nneg_bracket.second<<","<<p_bracket.first<<","<<p_bracket.second<<std::endl;
    depth = 0;
   }
 // Go deeper
 find_singular_p_value(nneg_initial, nneg_bracket_next, p_bracket_next, thresh, found_singular, depth);
 }

void loop_until_target_thickness(const double& initial_h, const double& target_h, const double& initial_dh)
  {
   oomph_info<<"#################################################"<<std::endl;
   oomph_info<<"Beginning Loop until Target Thickness"<<std::endl;
   // Check the input
   if( initial_dh == 0)
    {
     oomph_info<<"Initial dh is zero. This is not a valid dh."<<std::endl;
     oomph_info<<"Exiting."<<std::endl;
     return; 
    }

   // Check the input
   if( (target_h - initial_h) / initial_dh < 0 )
    {
     oomph_info<<"Initial thickness already past the target thickness: "<< target_h<<std::endl;
     oomph_info<<"Exiting."<<std::endl;
     return; 
    }
   else if( (target_h - initial_h) / initial_dh == 0 )
    {
     oomph_info<<"Initial pressure same as the target thickness: "<< target_h<<std::endl;
     oomph_info<<"Exiting."<<std::endl;
     return; 
    }
   // Set prev as the current mag
   oomph_info<<"Setting h_prev to be h."<<std::endl;
   TestSoln::h_prev = TestSoln::h;
 
   // Current pressure
   double current_h = initial_h, current_dh =initial_dh ,
          smallest_dh = initial_dh / pow(2.0,TestSoln::Max_number_hstep_halves);
   bool verbose = true;
   unsigned iter = 0 , maxiter = abs(round((target_h - initial_h)/smallest_dh));

   // Record of negative eigenvalues
   Vector<double> number_of_negative_eigenvalues;
   Vector<double> pressures;
  
   // Information about the loop
   if( verbose )
    {
     oomph_info<<"Aiming for target h: "<< target_h<<std::endl;
     oomph_info<<"Starting at initial h: "<< initial_h<<std::endl;
     oomph_info<<"Initial h step: "<< initial_dh<<std::endl;
     oomph_info<<"Minimum pressure step: "<< smallest_dh<<std::endl;
     oomph_info<<"Max iterations: "<< maxiter<<std::endl;
    }

   // In case we need a restart
   store_current_dof_values();
   // Loop until we reach target pressure
   do 
    {
     // Try a newton solve
     try 
      {
      // Output
      if( verbose )
       {
        oomph_info<<"Trying h: "<< TestSoln::h<<std::endl;
       }
       // Call the solve algorithm 
      newton_solve();
      }
     catch (NewtonSolverError)
      {
       // Go back a (half) step and continue
       if( iter !=0 )
        {
        oomph_info<<"h value "<< TestSoln::h << " did not converge. "
                  <<"Halving step size and trying again."<<std::endl;
        oomph_info<<"TestSoln::p: "<< TestSoln::p_mag<<std::endl;
         // Restore the the last converged step
         restore_dof_values();
         // Now store these!
         store_current_dof_values();
         // Now halve the step size and continue
         current_dh /= 2.0; 
         current_h -= current_dh;
         TestSoln::h = current_h;
         continue;
        }
       // If we failed on the first iteration
       else
       {
        oomph_info<<"Initial h did not converge. Please select a smaller"
                  <<" initial h."<<std::endl;
        oomph_info<<"Exiting."<<std::endl;
        return; 
       } 
      }
   // Output
//   if( verbose )
    {
     oomph_info<<"Current h: "<< current_h<<std::endl;
     oomph_info<<"TestSoln::h: "<< TestSoln::h<<std::endl;
     oomph_info<<"TestSoln::p: "<< TestSoln::p_mag<<std::endl;
     oomph_info<<"Current h_step: "<< current_dh<<std::endl;
     oomph_info<<"Current iteration: "<< iter<<"\n"<<std::endl;
    }

    // If we are dumping at every step 
    if(TestSoln::dump_at_every_step)
    {
    oomph_info<<"Trying dump"<<endl;
     // Dump the data 
     char filename[100];
     std::ofstream filestream;
     filestream.precision(15);
     sprintf(filename,"%s/kspe_circle_data%i-%f.dump",
             Doc_info.directory().c_str(),
             Doc_info.number(),
             TestSoln::p_mag
            );
     filestream.open(filename);
     dump(filestream);
     filestream.close();
    }
    
   // Document everything!
   TestSoln::h_prev = TestSoln::h;
   doc_solution();
   // Now say EVERYTHING
//   if(verbose)
   {
    oomph_info << std::endl;
    oomph_info << "---------------------------------------------" << std::endl;
    oomph_info << "Pressure (" << TestSoln::p_mag << ")" << std::endl;
    oomph_info << "Current h  (" << TestSoln::h << ")" << std::endl;
    oomph_info << "Poisson ratio (" << TestSoln::nu << ")" << std::endl;
    oomph_info << "Solution number (" <<Doc_info.number()-1 << ")" << std::endl;
    oomph_info << "---------------------------------------------" << std::endl;
    oomph_info << std::endl;
   }
   //In case we need a restart
   store_current_dof_values();
   // Increment 
   current_h += current_dh;
   TestSoln::h = current_h;
   ++iter;
   }
   while((initial_dh>0 ? current_h < target_h : current_h > target_h ) && current_dh/smallest_dh >= 1 && iter<maxiter);
   // HERE tidy up
 
   // Declare victory
   if((initial_dh>0 ? current_h > target_h : current_h < target_h ))
    {
     oomph_info<<"Completed loop in "<< iter<< " steps." << std::endl;
     oomph_info <<"COMPLETED SUCCESSFULLY."<<std::endl;
    }
   else if (current_dh/ smallest_dh <= 1)
    {
     oomph_info<<"Step size became too small. dh  : "<<current_dh<< std::endl;
    }
 }

// In this we keep the 'fvk nondimensional pressure constant' whilst modifying
// the thickness - which effectively keeps us steady on for FvK on the bif.
// diagram but will involve small changes for KS - should be favourable for
// remaining on the same branch and h stepping
void loop_until_target_thickness_keep_fvk_pressure_constant
  (const double& initial_h, const double& target_h, const double& initial_dh)
  {
   oomph_info<<"#################################################"<<std::endl;
   oomph_info<<"Beginning Loop until Target Thickness"<<std::endl;
   // Check the input
   if( initial_dh == 0)
    {
     oomph_info<<"Initial dh is zero. This is not a valid dh."<<std::endl;
     oomph_info<<"Exiting."<<std::endl;
     return; 
    }

   // Check the input
   if( (target_h - initial_h) / initial_dh < 0 )
    {
     oomph_info<<"Initial thickness already past the target thickness: "<< target_h<<std::endl;
     oomph_info<<"Exiting."<<std::endl;
     return; 
    }
   else if( (target_h - initial_h) / initial_dh == 0 )
    {
     oomph_info<<"Initial pressure same as the target thickness: "<< target_h<<std::endl;
     oomph_info<<"Exiting."<<std::endl;
     return; 
    }
   // Set prev as the current mag
   oomph_info<<"Setting h_prev to be h."<<std::endl;
   TestSoln::h_prev = TestSoln::h;
   TestSoln::p_prev = TestSoln::p_mag;
 
   // Current pressure
   double current_h = initial_h, current_dh =initial_dh ,
          smallest_dh = initial_dh / pow(2.0,TestSoln::Max_number_hstep_halves);
          
   bool verbose = true;
   unsigned iter = 0 , maxiter = abs(round((target_h - initial_h)/smallest_dh));

   // Record of negative eigenvalues
   Vector<double> number_of_negative_eigenvalues;
   Vector<double> pressures;
  
   // Information about the loop
   if( verbose )
    {
     oomph_info<<"Aiming for target h: "<< target_h<<std::endl;
     oomph_info<<"Starting at initial h: "<< initial_h<<std::endl;
     oomph_info<<"Initial h step: "<< initial_dh<<std::endl;
     oomph_info<<"Minimum pressure step: "<< smallest_dh<<std::endl;
     oomph_info<<"Max iterations: "<< maxiter<<std::endl;
    }

   // In case we need a restart
   store_current_dof_values();
   // Loop until we reach target pressure
   do 
    {
     // Try a newton solve
     try 
      {
      // Output
      if( verbose )
       {
        oomph_info<<"Trying h: "<< TestSoln::h<<std::endl;
       }
       // Call the solve algorithm 
      newton_solve();
      }
     catch (NewtonSolverError)
      {
       // Go back a (half) step and continue
       if( iter !=0 )
        {
        oomph_info<<"h value "<< TestSoln::h << " did not converge. "
                  <<"Halving step size and trying again."<<std::endl;
        oomph_info<<"TestSoln::p: "<< TestSoln::p_mag<<std::endl;
         // Restore the the last converged step
         restore_dof_values();
         // Now store these!
         store_current_dof_values();
         // Now halve the step size and continue
         current_dh /= 2.0; 
         current_h -= current_dh;
         TestSoln::p_mag =  TestSoln::p_mag * std::pow(current_h /TestSoln::h,3);
         TestSoln::h = current_h;
         oomph_info<<"Rewinding and halving step size"<< TestSoln::h <<std::endl;
         oomph_info<<"TestSoln::h "<< TestSoln::h <<std::endl;
         oomph_info<<"TestSoln::p: "<< TestSoln::p_mag<<std::endl;
         continue;
        }
       // If we failed on the first iteration
       else
       {
        oomph_info<<"Initial h did not converge. Please select a smaller"
                  <<" initial h."<<std::endl;
        oomph_info<<"Exiting."<<std::endl;
        return; 
       } 
      }
   // Output
//   if( verbose )
    {
     oomph_info<<"Current h: "<< current_h<<std::endl;
     oomph_info<<"TestSoln::h: "<< TestSoln::h<<std::endl;
     oomph_info<<"TestSoln::p: "<< TestSoln::p_mag<<std::endl;
     oomph_info<<"Current h_step: "<< current_dh<<std::endl;
     oomph_info<<"Current iteration: "<< iter<<"\n"<<std::endl;
    }

    // If we are dumping at every step 
    if(TestSoln::dump_at_every_step)
    {
    oomph_info<<"Trying dump"<<endl;
     // Dump the data 
     char filename[100];
     std::ofstream filestream;
     filestream.precision(15);
     sprintf(filename,"%s/kspe_circle_data%i-%f.dump",
             Doc_info.directory().c_str(),
             Doc_info.number(),
             TestSoln::p_mag
            );
     filestream.open(filename);
     dump(filestream);
     filestream.close();
    }
    
   // Document everything!
   TestSoln::p_prev = TestSoln::p_mag;
   doc_solution();
   // Now say EVERYTHING
//   if(verbose)
   {
    oomph_info << std::endl;
    oomph_info << "---------------------------------------------" << std::endl;
    oomph_info << "Pressure (" << TestSoln::p_mag << ")" << std::endl;
    oomph_info << "Current h  (" << TestSoln::h << ")" << std::endl;
    oomph_info << "Poisson ratio (" << TestSoln::nu << ")" << std::endl;
    oomph_info << "Solution number (" <<Doc_info.number()-1 << ")" << std::endl;
    oomph_info << "---------------------------------------------" << std::endl;
    oomph_info << std::endl;
   }
   //In case we need a restart
   store_current_dof_values();
   // Increment 
   current_h += current_dh;
   // Update p
   TestSoln::p_mag =  TestSoln::p_mag * std::pow(current_h / TestSoln::h,3);
   // Update h
   TestSoln::h = current_h;
   ++iter;
   }
   while((initial_dh>0 ? current_h < target_h : current_h > target_h ) && current_dh/smallest_dh >= 1 && iter<maxiter);
   // HERE tidy up
 
   // Declare victory
   if((initial_dh>0 ? current_h > target_h : current_h < target_h ))
    {
     oomph_info<<"Completed loop in "<< iter<< " steps." << std::endl;
     oomph_info <<"COMPLETED SUCCESSFULLY."<<std::endl;
    }
   else if (current_dh/ smallest_dh <= 1)
    {
     oomph_info<<"Step size became too small. dh  : "<<current_dh<< std::endl;
    }
 }
// Simple loop to target pressure
void loop_until_target_pressure(const double& initial_p, const double& target_p, const double& initial_dp)
  {
   oomph_info<<"#################################################"<<std::endl;
   oomph_info<<"Beginning Loop until Target Pressure"<<std::endl;
   // Check the input
   if( initial_dp == 0)
    {
     oomph_info<<"Initial dp is zero. This is not a valid dp."<<std::endl;
     oomph_info<<"Exiting."<<std::endl;
     return; 
    }

   // Check the input
   if( (target_p - initial_p) / initial_dp < 0 )
    {
     oomph_info<<"Initial pressure already past the target pressure: "<< target_p<<std::endl;
     oomph_info<<"Exiting."<<std::endl;
     return; 
    }
   else if( (target_p - initial_p) / initial_dp == 0 )
    {
     oomph_info<<"Initial pressure same as the target pressure: "<< target_p<<std::endl;
     oomph_info<<"Exiting."<<std::endl;
     return; 
    }
   // Set prev as the current mag
   oomph_info<<"Setting p_prev to be p_mag."<<std::endl;
   TestSoln::p_prev = TestSoln::p_mag;
 
   // Current pressure
   double current_p = initial_p, current_dp =initial_dp ,
          smallest_dp = initial_dp / pow(2.0,TestSoln::Max_number_pstep_halves);
   bool last_solve_success =false;
   bool verbose = true;
   unsigned iter = 0 , maxiter = abs(round((target_p - initial_p)/smallest_dp));

   // Record of negative eigenvalues
   Vector<double> number_of_negative_eigenvalues;
   Vector<double> pressures;
  
   // Information about the loop
   if( verbose )
    {
     oomph_info<<"Aiming for target pressure: "<< target_p<<std::endl;
     oomph_info<<"Starting at initial pressure: "<< initial_p<<std::endl;
     oomph_info<<"Initial pressure step: "<< initial_dp<<std::endl;
     oomph_info<<"Minimum pressure step: "<< smallest_dp<<std::endl;
     oomph_info<<"Max iterations: "<< maxiter<<std::endl;
    }

   // In case we need a restart
   store_current_dof_values();
   // Loop until we reach target pressure 
   do
    {
     // Try a newton solve
     try 
      {
      // Output
      if( verbose )
       {
        oomph_info<<"Trying pressure: "<< current_p<<std::endl;
       }
       // Call the solve algorithm 
      newton_solve();
      }
     catch (NewtonSolverError)
      {
       // Go back a (half) step and continue
       if( iter !=0 )
        {
        oomph_info<<"Pressure value "<<current_p << " did not converge."
                  <<"Halving step size and trying again."<<std::endl;
         // Restore the the last converged step
         restore_dof_values();
         // Now store these!
         store_current_dof_values();
         // Now halve the step size and continue
         current_dp /= 2.0; 
         // Go to halfway between
         current_p -= current_dp;
         // Revert 
         TestSoln::p_mag = current_p;
         continue;
        }
       // If we failed on the first iteration
       else
       {
        oomph_info<<"Initial pressure did not converge. Please select a smaller"
                  <<" initial pressure."<<std::endl;
        oomph_info<<"Exiting."<<std::endl;
        return; 
       } 
      }
    // Check to see if last solve was successful in a few steps
    if(current_dp < initial_dp && last_solve_success)
     {
      // Check if cnverging in few steps 
      if(Nnewton_iter_taken < 4)
       {
        // Double the step size
        current_dp *= 2;
       }
     }
   // Output
   if( verbose )
    {
     oomph_info<<"Current pressure: "<< current_p<<std::endl;
     oomph_info<<"Current pressure step: "<< current_dp<<std::endl;
     oomph_info<<"Current iteration: "<< iter<<"\n"<<std::endl;
    }

    // If we are dumping at every step 
    if(TestSoln::dump_at_every_step)
    {
    oomph_info<<"Trying dump"<<endl;
     // Dump the data 
     char filename[100];
     std::ofstream filestream;
     filestream.precision(15);
     sprintf(filename,"%s/kspe_circle_data%i-%f.dump",
             Doc_info.directory().c_str(),
             Doc_info.number(),
             TestSoln::p_mag
            );
     filestream.open(filename);
     dump(filestream);
     filestream.close();
    }
    
   // Document everything!
   TestSoln::p_prev = TestSoln::p_mag;
   doc_solution();
   // Now say EVERYTHING
//   if(verbose)
   {
    oomph_info << std::endl;
    oomph_info << "---------------------------------------------" << std::endl;
    oomph_info << "Pressure (" << TestSoln::p_mag << ")" << std::endl;
    oomph_info << "Current h  (" << TestSoln::h << ")" << std::endl;
    oomph_info << "Poisson ratio (" << TestSoln::nu << ")" << std::endl;
    oomph_info << "Solution number (" <<Doc_info.number()-1 << ")" << std::endl;
    oomph_info << "---------------------------------------------" << std::endl;
    oomph_info << std::endl;
   }
   //In case we need a restart
   store_current_dof_values();
   
   // If we are doing the eigen problem
   if(TestSoln::do_eigen)
    {
     // Solve the eigen problem
     unsigned nnegeig = 0;
     bool found_singular = false;
     solve_eigen_problem(nnegeig,found_singular,TestSoln::Singular_tolerance);

     // Perform push backs
     number_of_negative_eigenvalues.push_back(nnegeig);
     pressures.push_back(TestSoln::p_mag);

     // Output current information
     if(verbose)
     {
      oomph_info<<"Current number of negative Eigenvalues: " <<nnegeig<<std::endl;
      oomph_info<<"Current previous Eigenvalues: " 
                <<number_of_negative_eigenvalues<<std::endl;
     }
     // Get length
     unsigned length =pressures.size();
     oomph_info<<"Length is "<<length<<std::endl;
     // Now check to see if we need to search
     if(length>1)
      {
      if((number_of_negative_eigenvalues[length-2]) != (number_of_negative_eigenvalues[length-1]))
       {
        // Find the singular value and add it on
        int depth = TestSoln::max_depth;
        typedef std::pair<unsigned,unsigned> Negative_eig_bracket;
        typedef std::pair<double,double> Pressure_bracket;
        Negative_eig_bracket nneg_bracket(number_of_negative_eigenvalues[length-2],
          number_of_negative_eigenvalues[length-1]);
        Pressure_bracket p_bracket(pressures[length-2],pressures[length-1]);
        // Now find the singular p
        find_singular_p_value(nneg_bracket.first,nneg_bracket,p_bracket,
          TestSoln::Singular_tolerance,found_singular, depth);
        // Report our findings
        if (found_singular)
          {oomph_info<<"Found singular value at pressure = "<<TestSoln::p_mag<< 
                    std::endl;}
        else 
          {oomph_info<<"Max depth reached before singular value found"<<std::endl;}

       /// If we only want the first!
       if( TestSoln::only_first_singular_value)
        { 
        oomph_info<<"First singular value investigated, turning do_eigen off."
                 <<std::endl;
         TestSoln::do_eigen =  false;
        }
       } 
      }
    } 

    // Increment 
    current_p += current_dp;
    TestSoln::p_mag = current_p;
    ++iter;
    last_solve_success = true;
    }
   while((initial_dp>0 ? current_p < target_p : current_p > target_p ) && current_dp/smallest_dp >= 1 && iter<maxiter);
   // HERE tidy up
 
   // Declare victory
   if((initial_dp>0 ? current_p > target_p : current_p < target_p ))
    {
     oomph_info<<"Completed loop in "<< iter<< " steps." << std::endl;
     oomph_info <<"COMPLETED SUCCESSFULLY."<<std::endl;
    }
   else if (current_dp/ smallest_dp <= 1)
    {
     oomph_info<<"Step size became too small. dp  : "<<current_dp<< std::endl;
    }
   else 
    {
     oomph_info<<"Max iterations exceeded"<< std::endl;
    }
    oomph_info<<"#################################################\n"<<std::endl;
  }

// Solve the eigen problem - always finds the eigenvalue in the bracket that is
// closest to zero
void solve_eigen_problem(unsigned& number_of_negative_eigenvalues, 
  bool& found_singular, const double& thresh)
 {
  oomph_info<<"#################################################\n"<<std::endl;
  oomph_info << "Doing eigenproblem" << std::endl;
  // How many eigenvalues do we try to get
  // Set default to positive infinity
  complex<double> default_value;
  default_value.real()=std::numeric_limits<double>::infinity();
  default_value.imag()=0;
  DoubleVector broken_symmetry_vector;

  Vector<complex<double> > eigenvalues(TestSoln::neig,default_value);
  Vector<DoubleVector> eigenvectors;
  solve_eigenproblem(TestSoln::neig,eigenvalues,eigenvectors);

  // Determine how many found
  unsigned nfound=eigenvalues.size();
  oomph_info << "Found "<<nfound <<" eigenvalues.\n";

  //Now show the nonzero eigenvalues
  if(!TestSoln::use_direct_solver)
   {oomph_info << "Eigenvalues are: " << std::endl;}
  double min=DBL_MAX;
  unsigned i_min=0;
  for (unsigned ieig=0;ieig<nfound;ieig++)
   {
    // Check to see if its the minimum
    if((eigenvalues[ieig]).real()<std::numeric_limits<double>::max())
     {
     if(!TestSoln::use_direct_solver)
      {oomph_info << ieig << " " << eigenvalues[ieig] << " "; }
     if (std::abs(eigenvalues[ieig])<min)
      {
       i_min=ieig;
       min=std::abs(eigenvalues[ieig]);
      }

     // If negative, increment 
     if ((eigenvalues[ieig]).real()<0.0)
      {++(number_of_negative_eigenvalues);}
     // Output complex eigenvalues
     if (std::abs((eigenvalues[ieig]).imag())>1e-14)
      {oomph_info<<"Complex ev:"<<(eigenvalues[ieig]).real()<<(eigenvalues[ieig]).imag()<<"\n";}
     }
   }
  oomph_info << std::endl << "Minimum eigenvalue: " << i_min << " : " << min 
             << std::endl;
  
  oomph_info << std::endl << "Number of negative eigenvalues: "
             << number_of_negative_eigenvalues << " : " << min << std::endl;
  // Check its not postive infinity
  if((eigenvalues[i_min]).real()<std::numeric_limits<double>::max())
   {
   if (std::abs(eigenvalues[i_min])<thresh)
    {
    // Set the bool to true
    found_singular=true;
    broken_symmetry_vector=eigenvectors[i_min];
    // ...and assign real part to dofs
    add_eigenvector_to_dofs(TestSoln::Singular_prefactor,eigenvectors[i_min]);
    oomph_info << "Adding for singular eigenvector for eigenvalue " 
               << eigenvalues[i_min]
               << " to solution " << Doc_info.number() << std::endl;
    doc_solution();
    oomph_info << std::endl;
    oomph_info << "---------------------------------------------" << std::endl;
    oomph_info << "Pressure (" << TestSoln::p_mag << ")" << std::endl;
    oomph_info << "Current h  (" << TestSoln::h << ")" << std::endl;
    oomph_info << "Poisson ratio (" << TestSoln::nu << ")" << std::endl;
    oomph_info << "Solution number (" <<Doc_info.number()-1 << ")" << std::endl;
    oomph_info << "---------------------------------------------" << std::endl;
    oomph_info << std::endl;

    newton_solve();
    oomph_info << "Newton solve for singular eigenvector for eigenvalue " 
               << eigenvalues[i_min]
               << " to solution " << Doc_info.number() << std::endl;
    doc_solution();

    oomph_info << std::endl;
    oomph_info << "---------------------------------------------" << std::endl;
    oomph_info << "Pressure (" << TestSoln::p_mag << ")" << std::endl;
    oomph_info << "Current h  (" << TestSoln::h << ")" << std::endl;
    oomph_info << "Poisson ratio (" << TestSoln::nu << ")" << std::endl;
    oomph_info << "Solution number (" <<Doc_info.number()-1 << ")" << std::endl;
    oomph_info << "---------------------------------------------" << std::endl;
    oomph_info << std::endl;
    }
  }
  oomph_info<<"#################################################\n"<<std::endl;
 }

void print_jacobian();

private:
/// Helper function to apply boundary conditions
void apply_boundary_conditions();

/// \short Helper function to (re-)set boundary condition
/// and complete the build of  all elements
void complete_problem_setup();

/// Trace file to document norm of solution
ofstream Trace_file;

// Keep track of boundary ids
enum
{
 Outer_boundary0 = 0,
 Outer_boundary1 = 1,
 Inner_boundary0 = 2,
 Inner_boundary1 = 3,
 Inner_boundary2 = 4
};

/// Element area data
double Element_area;

/// Bool investigate the eigen cross?
bool Investigate_eigen_cross;

/// Bool use the centre point to pin it?
bool Use_centre_point_as_pin;

/// Bool is wrinkled
bool Is_wrinkled;

void upgrade_edge_elements_to_curve(const unsigned &b, Mesh* const & 
 bulk_mesh_pt);

void rotate_edge_degrees_of_freedom(Mesh* const &bulk_mesh_pt);

/// Pointer to "bulk" mesh
TriangleMesh<ELEMENT>* Bulk_mesh_pt;

/// Pointer to "surface" mesh
Mesh* Surface_mesh_pt;

/// Pointers to mesh objects
// Outer Ellipse objects
Ellipse* Outer_boundary_ellipse_pt;
TriangleMeshClosedCurve* Outer_boundary_pt;
Vector<TriangleMeshCurveSection*> Outer_curvilinear_boundary_pt;
// Inner Ellipse objects
Ellipse* Inner_boundary_ellipse_pt;
Vector<TriangleMeshClosedCurve*> Inner_boundary_pt;
Vector<TriangleMeshCurveSection*> Inner_curvilinear_boundary_pt;
// Inner Line objects
Vector<TriangleMeshOpenCurve *> Inner_open_boundaries_pt;
TriangleMeshPolyLine* Boundary2_pt;
Vector<TriangleMeshCurveSection *> Internal_curve_section1_pt;
}; // end_of_problem_class


template<class ELEMENT>
UnstructuredFvKProblem<ELEMENT>::UnstructuredFvKProblem(double element_area)
:
Element_area(element_area),
Investigate_eigen_cross(false),
Use_centre_point_as_pin(true),
Is_wrinkled(false)
{
// Parameters for internal boundary
double inner_radius = 4./5. ; // 0.7 ; // 4./5.;
double outer_element_area = element_area;
double inner_element_area = element_area * 10. ;//3.; // 10
// Preinitialise
Vector<double> zeta(1);
Vector<double> posn(2);

//Outer boundary
//--------------

double A = 1.0;
double B = 1.0;
Outer_boundary_ellipse_pt = new Ellipse(A, B);
// Initialize
Outer_boundary_pt = 0;
// Boundary specified in two parts
Outer_curvilinear_boundary_pt.resize(2);

//First bit
double zeta_start = 0.0;
double zeta_end = MathematicalConstants::Pi;
unsigned nsegment = (int)(MathematicalConstants::Pi/sqrt(outer_element_area));
Outer_curvilinear_boundary_pt[0] =
new TriangleMeshCurviLine(Outer_boundary_ellipse_pt, zeta_start,
zeta_end, nsegment, Outer_boundary0);

//Second bit
zeta_start = MathematicalConstants::Pi;
zeta_end = 2.0*MathematicalConstants::Pi;
nsegment = (int)(MathematicalConstants::Pi/sqrt(outer_element_area));
Outer_curvilinear_boundary_pt[1] =
new TriangleMeshCurviLine(Outer_boundary_ellipse_pt, zeta_start,
zeta_end, nsegment, Outer_boundary1);

Outer_boundary_pt =
new TriangleMeshClosedCurve(Outer_curvilinear_boundary_pt);


// Inner boundary
//--------------

A = inner_radius;
B = inner_radius;
Inner_boundary_ellipse_pt = new Ellipse(A, B);
// Boundary specified in two parts
Inner_boundary_pt.resize(1);
Inner_curvilinear_boundary_pt.resize(2);

//First bit
zeta_start = 0.0;
zeta_end = MathematicalConstants::Pi;
nsegment = (int)(inner_radius*MathematicalConstants::Pi/sqrt(outer_element_area));
Inner_curvilinear_boundary_pt[0] =
new TriangleMeshCurviLine(Inner_boundary_ellipse_pt, zeta_start,
zeta_end, nsegment, Inner_boundary0);

//Second bit
zeta_start = MathematicalConstants::Pi;
zeta_end = 2.0*MathematicalConstants::Pi;
nsegment = (int)(inner_radius*MathematicalConstants::Pi/sqrt(outer_element_area));
Inner_curvilinear_boundary_pt[1] =
new TriangleMeshCurviLine(Inner_boundary_ellipse_pt, zeta_start,
zeta_end, nsegment, Inner_boundary1);

 //Combine to internal curvilinear boundary
Inner_boundary_pt[0] =
new TriangleMeshClosedCurve(Inner_curvilinear_boundary_pt);

// Internal open boundaries
// Total number of open curves in the domain
unsigned n_open_curves = 1;
// We want internal open curves
Inner_open_boundaries_pt.resize(n_open_curves);

//Create the mesh
//---------------
//Create mesh parameters object
TriangleMeshParameters mesh_parameters(Outer_boundary_pt);
Boundary2_pt = 0;

// Internal bit - this means we can have a boundary which is just the centre
if(Use_centre_point_as_pin)
{
// We start by creating the internal boundaries
// The boundary 2 is defined by its two vertices
// Open curve 1
 Vector<Vector<double> > vertices(2,Vector<double>(2,0.0));
 vertices[0][0] = 0.0;
 vertices[0][1] = 0.0;
 
 vertices[1][0] = 0.5;
 vertices[1][1] = 0.0;
 unsigned boundary_id = Inner_boundary2;
 Boundary2_pt =
   new TriangleMeshPolyLine(vertices, boundary_id);
// 
// Each inteTriangleMeshPolyLine *boundary2_ptrnal open curve is defined by a vector of
// TriangleMeshCurveSection,
// on this example we only need one curve section for each internal boundary
 Internal_curve_section1_pt.resize(1);
 Internal_curve_section1_pt[0] = Boundary2_pt;

// The open curve that define this boundary is composed of just one
// curve section
 Inner_open_boundaries_pt[0] =
    new TriangleMeshOpenCurve(Internal_curve_section1_pt);

// Specify the internal open boundaries
mesh_parameters.internal_open_curves_pt() = Inner_open_boundaries_pt;
}

// Add internal boundary
mesh_parameters.internal_closed_curve_pt() = Inner_boundary_pt;
// Add region
Vector<double> outer(2);
outer[0] =  (1.0+inner_radius)/2.;
mesh_parameters.add_region_coordinates(1,outer);
// Specify the element area
mesh_parameters.set_target_area_for_region(1,outer_element_area);
mesh_parameters.element_area() = inner_element_area;

// Build an assign bulk mesh
Bulk_mesh_pt=new TriangleMesh<ELEMENT>(mesh_parameters);

// Create "surface mesh" that will contain only the prescribed-traction
// elements. The constructor creates the mesh without adding any nodes
// elements etc.
Surface_mesh_pt =  new Mesh;

//Add two submeshes to problem
add_sub_mesh(Bulk_mesh_pt);
add_sub_mesh(Surface_mesh_pt);

// Combine submeshes into a single Mesh
build_global_mesh();

// Curved Edge upgrade
upgrade_edge_elements_to_curve(0,Bulk_mesh_pt);
upgrade_edge_elements_to_curve(1,Bulk_mesh_pt);
 
// Rotate degrees of freedom
rotate_edge_degrees_of_freedom(Bulk_mesh_pt);

// Store number of bulk elements
complete_problem_setup();

char filename[100];
sprintf(filename, "RESLT/trace.dat");
Trace_file.open(filename);

oomph_info << "Number of equations: "
        << assign_eqn_numbers() << '\n';

// Dump the mesh
if(TestSoln::Dump_triangulateio)
{
 std::ofstream os;
 std::string nom ("dump");
 TriangulateIO tio = Bulk_mesh_pt->triangulateio_representation();
 Bulk_mesh_pt->write_triangulateio(tio,nom);
}
}



//==start_of_complete======================================================
/// Set boundary condition exactly, and complete the build of 
/// all elements
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::complete_problem_setup()
{   
unsigned nbound = Outer_boundary1 + 1;
 // Upcast to current element
 ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(0));
// Set the boundary conditions for problem: All nodes are
// free by default -- just pin the ones that have Dirichlet conditions
// here. 
//Just loop over outer boundary since inner boundary doesn't have boundary
//conditions

// Pin centre node
if(Use_centre_point_as_pin)
 pin_centre_point();
else 
 pin_first_point();

// Now for the rest of b/c
for(unsigned ibound=0;ibound<nbound;ibound++)
{
unsigned num_nod=Bulk_mesh_pt->nboundary_node(ibound);
for (unsigned inod=0;inod<num_nod;inod++)
{
 // Get nod
 Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
 Vector<double> x(2,0.0),w(18,0.0);
 x[0]=nod_pt->x(0);
 x[1]=nod_pt->x(1);

 // Pin unknown values (everything except for the second normal derivative)
 for(unsigned i=2;i<3;++i)
  {
   for(unsigned l=0;l<6;++l)
    {
    // Don't pin the second normal
    if(l != 3)  
     {
     unsigned index=el_pt->u_index_koiter_model(l,i);
     nod_pt->pin(index);
     nod_pt->set_value(index,w[index]);
     }
    }
  }
 }
} // end loop over boundaries 


// Complete the build of all elements so they are fully functional
unsigned n_element = Bulk_mesh_pt->nelement();
for(unsigned e=0;e<n_element;e++)
{
// Upcast from GeneralisedElement to the present element
ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

//Set the pressure function pointers and the physical constants
el_pt->pressure_fct_pt() = &TestSoln::get_pressure;
el_pt->nu_pt() = &TestSoln::nu;
el_pt->d_pressure_dn_fct_pt() = &TestSoln::get_d_pressure_dn;
el_pt->d_pressure_dr_fct_pt() = &TestSoln::get_d_pressure_dr;
el_pt->d_pressure_d_grad_u_fct_pt() = &TestSoln::get_d_pressure_d_grad_u;
el_pt->error_metric_fct_pt() = &TestSoln::error_metric;
el_pt->multiple_error_metric_fct_pt() = &TestSoln::fourier_transform_metric;
el_pt->thickness_pt() = &TestSoln::h;
// Use a finite difference jacobian
if(TestSoln::use_fd_jacobian)
 el_pt->enable_finite_difference_jacobian();
}

// Re-apply Dirichlet boundary conditions (projection ignores
// boundary conditions!)
apply_boundary_conditions();
}

//==start_of_apply_bc=====================================================
/// Helper function to apply boundary conditions
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::apply_boundary_conditions()
{

// Loop over all boundary nodes
//Just loop over outer boundary conditions
unsigned nbound = Outer_boundary1 + 1;

// Pin centre node
if(Use_centre_point_as_pin)
 pin_centre_point();
else 
 pin_first_point();

// Upcast to first element
ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(0));

// Pin Rotations
// Loop only over the first outer boundary (to avoid visiting the relevant node 
// twice) 
for(unsigned ibound=0;ibound<Outer_boundary1;ibound++)
{
 // Local block
 unsigned num_nod=Bulk_mesh_pt->nboundary_node(ibound);
 for(unsigned inod=0;inod<num_nod;++inod)
  {
  // First Node pt
  Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
  double x = nod_pt->x(0), y = nod_pt->x(1);
  // Pin Rotations by pinning Uy at Node at (1,0) 
  if(nod_pt->is_on_boundary(Outer_boundary0) && nod_pt->is_on_boundary(Outer_boundary1)
     && x > 0 )
    {
     oomph_info<<"Found node: ("<<x<<","<<y<<")\n Pinning to zero."<<std::endl;
     nod_pt->pin(6);
     nod_pt->set_value(6,0.0);
    }
  }
}

// Loop Boundaries
for(unsigned ibound=0;ibound<nbound;ibound++)
{
// Local block
unsigned num_nod=Bulk_mesh_pt->nboundary_node(ibound);
for (unsigned inod=0;inod<num_nod;inod++)
{
 // Get node
 Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
 
 // Extract nodal coordinates from node:
 Vector<double> x(2),w(18);
 x[0]=nod_pt->x(0);
 x[1]=nod_pt->x(1);

 // Pin unknown values (everything except for the second normal derivative)
 for(unsigned i=2;i<3;++i)
  {
   for(unsigned l=0;l<6;++l)
    {
    // Don't pin the second normal
    if(l != 3)  
     {
     unsigned index=el_pt->u_index_koiter_model(l,i);
     nod_pt->pin(index);
     nod_pt->set_value(index,w[index]);
     }
    }
  }
}
} 

} // end set bc

template <class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::
upgrade_edge_elements_to_curve(const unsigned &b, Mesh* const &bulk_mesh_pt) 
{
 // How many bulk elements adjacent to boundary b
 unsigned n_element = bulk_mesh_pt-> nboundary_element(b);
 // These depend on the boundary we are on
 CurvilineGeomObject* parametric_curve_pt; 
 // Define the functions for each part of the boundary
 switch (b)
  {
   // Upper boundary
   case 0:
    parametric_curve_pt = &TestSoln::parametric_curve_top;
   break;

   // Lower boundary
   case 1:
    parametric_curve_pt = &TestSoln::parametric_curve_bottom;
   break;

   default:
    throw OomphLibError(
     "I have encountered a boundary number that I wasn't expecting. This is very\
 peculiar.",
     "UnstructuredFvKProblem::upgrade_edge_elements_to_curve(...)",
     OOMPH_EXCEPTION_LOCATION);
   break;
  }
 
 // Loop over the bulk elements adjacent to boundary b
 for(unsigned e=0; e<n_element; e++)
  {
   // Get pointer to bulk element adjacent to b
   ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(
    bulk_mesh_pt->boundary_element_pt(b,e));
   
   // Loop over nodes
   const unsigned nnode=3;
   unsigned index_of_interior_node=3;

   // The edge that is curved
   MyC1CurvedElements::Edge edge;

   // Vertices positions
   Vector<Vector<double> > xn(3,Vector<double>(2,0.0));
 
   // Loop nodes
   for(unsigned n=0;n<nnode;++n)
    {
     // If it is on boundary
     Node* nod_pt = bulk_el_pt->node_pt(n);
     if(nod_pt->is_on_boundary(Outer_boundary0) || nod_pt->is_on_boundary(Outer_boundary1))
      {
       xn[n][0]=nod_pt->x(0);
       xn[n][1]=nod_pt->x(1);
      }
     // The edge is denoted by the index of the  opposite (interior) node
     else {index_of_interior_node = n;}
    }
   // Initialise s_ubar s_obar
   double s_ubar, s_obar;

   // s at the next (cyclic) node after interior
   s_ubar = parametric_curve_pt->get_zeta(xn[(index_of_interior_node+1) % 3]);
   // s at the previous (cyclic) node before interior
   s_obar = parametric_curve_pt->get_zeta(xn[(index_of_interior_node+2) % 3]);

   // Assign edge case
   switch(index_of_interior_node)
    {
     case 0: edge= MyC1CurvedElements::zero; 
      break;
     case 1: edge= MyC1CurvedElements::one; 
      break;
     case 2: edge= MyC1CurvedElements::two; 
      break;
     // Should break it here HERE
     default: edge= MyC1CurvedElements::none; 
      throw OomphLibError(
       "The edge number has been set to a value greater than two: either we have\
 quadrilateral elements or more likely the index_of_interior_node was never set\
 and remains at its default value.",
       "UnstructuredFvKProblem::upgrade_edge_elements_to_curve(...)",
       OOMPH_EXCEPTION_LOCATION);
      break;
     }
   if (s_ubar>s_obar)
    {
     oomph_info <<"Apparent clockwise direction of parametric coordinate."
                <<"This will probably result in an inverted element."
                <<"s_start="<<s_ubar<<"; s_end ="<<s_obar<<std::endl;
     throw OomphLibError(
       "The Edge coordinate appears to be decreasing from s_start to s_end. \
Either the parametric boundary is defined to be clockwise (a no-no) or \
the mesh has returned an inverted element (less likely)",
       "UnstructuredFvKProblem::upgrade_edge_elements_to_curve(...)",
       OOMPH_EXCEPTION_LOCATION);
    }

   // Upgrade it
   bulk_el_pt->upgrade_to_curved_element(edge,s_ubar,s_obar,
    parametric_curve_pt);     
  }
}// end upgrade elements

template <class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::
rotate_edge_degrees_of_freedom( Mesh* const &bulk_mesh_pt)
{
 // How many bulk elements
 unsigned n_element = bulk_mesh_pt-> nelement();
 
 // Loop over the bulk elements adjacent to boundary b
 for(unsigned e=0; e<n_element; e++)
  {
   // Get pointer to bulk element adjacent to b
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
 
//   el_pt->pin_all_in_plane_dofs();

   // Loop nodes  HERE
   unsigned nnode = 3;//nnode();
   unsigned nbnode=0 ;
   // Count the number of boundary nodes
   for (unsigned n=0; n<nnode;++n)
     {
      // Check if it is on the curved boundaries
      bool on_boundary_0=el_pt->node_pt(n)->is_on_boundary(Outer_boundary0);
      bool on_boundary_1=el_pt->node_pt(n)->is_on_boundary(Outer_boundary1);
      if(on_boundary_0 || on_boundary_1 )
       {nbnode+=1;}
     }

   // Now if we have nodes on boundary 
   if(nbnode>0)
    {
     // Set up vector
     Vector<unsigned> bnode (nbnode,0);
     unsigned inode(0);

     // Fill in the bnode Vector
     for (unsigned n=0; n<nnode;++n)
      {
       // Check it isn't on an internal boundary
      // Check if it is on the curved boundaries
      bool on_boundary_0=el_pt->node_pt(n)->is_on_boundary(Outer_boundary0);
      bool on_boundary_1=el_pt->node_pt(n)->is_on_boundary(Outer_boundary1);
      if(on_boundary_0 || on_boundary_1 )
       {
       // If it is on the boundary
       if(el_pt->node_pt(n)->is_on_boundary())
        {
         // Set up the Vector
         bnode[inode]=n;
         ++inode;
        }
       }
      }
    // Output that we have found element HERE
//    std::cout<<"Element "<<e<<" has "<<bnode<< " nodes on the boundary.\n";

    el_pt->set_up_rotated_dofs(nbnode,bnode,&TestSoln::get_normal_and_tangent);
   // Now rotate the nodes
   }
 }
}// end create traction elements

//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::doc_solution(const 
                                                    std::string& comment)
{ 
ofstream some_file;
char filename[100];

// Number of plot points
unsigned npts = 2;

if(Use_centre_point_as_pin)
 output_centre_point();

sprintf(filename,"RESLT/coarse_soln%i-%f.dat",Doc_info.number(),Element_area);
some_file.open(filename);
Bulk_mesh_pt->output(some_file,npts); 
some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \"" 
       << comment << "\"\n";
some_file.close();

//Document the extras if we are not worried about disk space
if(!TestSoln::save_disk_space)
{
// Number of plot points
npts = 6;

sprintf(filename,"RESLT/soln%i-%f.dat",Doc_info.number(),Element_area);
some_file.open(filename);
Bulk_mesh_pt->output(some_file,npts); 
some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \"" 
       << comment << "\"\n";
some_file.close();

// Output boundaries
//------------------
sprintf(filename,"RESLT/boundaries%i-%f.dat",Doc_info.number(),Element_area);
some_file.open(filename);
Bulk_mesh_pt->output_boundaries(some_file);
some_file.close();

// Output regions
unsigned n_region = Bulk_mesh_pt->nregion();
if (n_region > 1)
{
for (unsigned r = 0; r < n_region; r++)
{
 //Attempt to output elements in different regions
 sprintf(filename,"RESLT/region%i%i-%f.dat",r,Doc_info.number(),
   Element_area);
 some_file.open(filename);
 unsigned nel = Bulk_mesh_pt->nregion_element(r);
 for (unsigned e = 0; e < nel; e++)
  {
   Bulk_mesh_pt->region_element_pt(r,e)->output(some_file,npts);
  }
 some_file.close();
}
}
}
// // Doc error and return of the square of the L2 error
// //---------------------------------------------------
// //double error,norm,dummy_error,zero_norm;
 double dummy_error,zero_norm;
 if(!TestSoln::save_disk_space)
  sprintf(filename,"RESLT/error%i-%f.dat",Doc_info.number(),Element_area);
 else 
  sprintf(filename,"RESLT/error.dat");

 some_file.open(filename);
 
 Bulk_mesh_pt->compute_error(some_file,TestSoln::get_exact_w,
                        dummy_error,zero_norm);
 some_file.close();

 // Now do the multiple error function
 if(!TestSoln::save_disk_space)
  sprintf(filename,"RESLT/errors%i-%f.dat",Doc_info.number(),Element_area);
 else 
  sprintf(filename,"RESLT/error.dat");

 some_file.open(filename);

 Vector<double> dummy_errors(TestSoln::Number_of_fourier_terms), 
  zero_norms(TestSoln::Number_of_fourier_terms);
 Bulk_mesh_pt->compute_error(some_file,TestSoln::get_exact_w,
                        dummy_errors,zero_norms);
 some_file.close();
 
 // Tolerance for being on wrinkled branch
 const double tol = 1e-6;
 if(fabs(dummy_error) > tol)
  {Is_wrinkled = true; }

 // Doc L2 error and norm of solution
 oomph_info << "Error squared of computed solution: " << /*sqrt*/(dummy_error)<< std::endl;
 oomph_info << "Norm squared of computed solution: "  << /*sqrt*/(zero_norm)  << std::endl;
 oomph_info << "Fourier coefficients of computed solution: " << /*sqrt*/(dummy_errors)<< std::endl;
 oomph_info << "Fourier norms of computed solution: " << /*sqrt*/(zero_norms)<< std::endl;
 
 Trace_file << TestSoln::p_mag << " " << "\n ";

// Doc error and return of the square of the L2 error
//---------------------------------------------------
sprintf(filename,"RESLT/L2-norm%i-%f.dat",
        Doc_info.number(),
        Element_area);
some_file.open(filename);

some_file<<"### L2 Norm\n";
some_file<<"##  Format: err^2 norm^2 log(err/norm) \n";
// Print error in prescribed format
some_file<< dummy_error <<" "<< zero_norm <<" ";

some_file.close();

// Increment the doc_info number
Doc_info.number()++;

} // end of doc

//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::print_jacobian()
{
 // Print the Jacobian
 describe_dofs();
 oomph_info << "----------PRINTING JACOBIAN-----------" << std::endl;
 DoubleVector residuals;
 CRDoubleMatrix jac;
 get_jacobian(residuals,jac);
 jac.sort_entries();
 
 // Open filestream, high precision
 std::ofstream filestream;
 filestream.precision(15);
 // Filename
 char jacobian_filename[100];
 sprintf(jacobian_filename, "%s/%s.csv",Doc_info.directory().c_str(),
  (TestSoln::use_fd_jacobian? "jacobian_fd":"jacobian_exact"));
 filestream.open(jacobian_filename);
 filestream.precision(15);
 filestream << " # name: jac" << "\n"
            << "  # type: sparse matrix"<< "\n"
            << "  # nnz:" << jac.nnz() <<"\n"
            << "  # rows:" << jac.nrow() <<"\n"
            << "  # columns:" << jac.ncol() <<std::endl;

 // Output
  for (unsigned long i=0;i<jac.nrow();i++)
   {
    for (long j=jac.row_start()[i];j<jac.row_start()[i+1];j++)
     {
      // Output transpose so columns are grouped
      filestream << jac.column_index()[j]+1 << " " << i+1 << " " << jac.value()[j]
              << std::endl;
     }
   }
}

// Extend the namespace to write a checkpoint
namespace TestSoln{
// Problem_pt
UnstructuredFvKProblem<LargeDisplacementPlateC1CurvedBellElement<2,2,3,
KoiterSteigmannPlateEquations> >* problem_pt=0;

static void write_checkpoint()
 {
  // Dump the data 
  char filename[100];
  std::ofstream filestream;
  filestream.precision(15);
  sprintf(filename,"checkpoint.dump" );
  filestream.open(filename);
  problem_pt->restore_dof_values();
  problem_pt->dump(filestream);
  filestream.close();

  sprintf(filename,"checkpoint_p" );
  filestream.open(filename);
  filestream << TestSoln::p_prev;
  filestream.close();

  sprintf(filename,"checkpoint_nu" );
  filestream.open(filename);
  filestream << TestSoln::nu;
  filestream.close();

  sprintf(filename,"checkpoint_h" );
  filestream.open(filename);
  filestream << TestSoln::h;
  filestream.close();

  sprintf(filename,"checkpoint_nsoln" );
  filestream.open(filename);
  filestream <<problem_pt->Doc_info.number();
  filestream.close();
 }

 static void signal_handler(int signum)
  {
   // Write the checkpoint and exit
   write_checkpoint();
   exit(0);
  }

}

//=======start_of_main========================================
///Driver code for demo of inline triangle mesh generation
//============================================================
int main(int argc, char **argv)
{
 // Catch soft kill signal from Condor and checkpoint if we get evicted 
 signal(SIGTERM, TestSoln::signal_handler);

 feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);// | FE_UNDERFLOW) ;
 // Underflow error in ARPACK | FE_UNDERFLOW);

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Validation?
 CommandLineArgs::specify_command_line_flag("--validation");
 
 // Use FD Jacobian?
 CommandLineArgs::specify_command_line_flag("--use_fd_jacobian");

 // Dump Triangulateio?
 CommandLineArgs::specify_command_line_flag("--dump_triangulateio");

 // Physical Parameters
 CommandLineArgs::specify_command_line_flag("--p", &TestSoln::p_mag);

 CommandLineArgs::specify_command_line_flag("--nu", &TestSoln::nu);

 CommandLineArgs::specify_command_line_flag("--h", &TestSoln::h);

 // Stepping parameters
 unsigned n_step=2;
 CommandLineArgs::specify_command_line_flag("--n_step", &n_step);

 unsigned n_h_step = 0;
 CommandLineArgs::specify_command_line_flag("--n_h_step", &n_h_step);
 CommandLineArgs::specify_command_line_flag("--Max_stepsize_halves"
  , &TestSoln::Max_number_pstep_halves);
 CommandLineArgs::specify_command_line_flag("--keep_vk_pressure_constant");

 // P_Start
 double p_start=1e-8;
 // E = 1.36 e 6
 // h/L = 4.25e-3
 // So p*= 800 => 0.2 * 1.34* 4.25* 1e3 = 800 
 CommandLineArgs::specify_command_line_flag("--p_start", &p_start);

 // P_end
 double p_end=2e-5;//1.01e-3;
 CommandLineArgs::specify_command_line_flag("--p_end", &p_end);

 double h_end = 0.015;
 CommandLineArgs::specify_command_line_flag("--h_end", &h_end);
 
 
 // Element Area
 double element_area=0.1;
 CommandLineArgs::specify_command_line_flag("--element_area", &element_area);


 // Operational Flags
 CommandLineArgs::specify_command_line_flag("--save_disk_space");

 CommandLineArgs::specify_command_line_flag("--only_do_first_singular");

 CommandLineArgs::specify_command_line_flag("--use_jfm_nondim");

 CommandLineArgs::specify_command_line_flag("--use_LAPACK");

 CommandLineArgs::specify_command_line_flag("--no_eigen");

 unsigned max_depth;
 CommandLineArgs::specify_command_line_flag("--max_singular_steps", &max_depth);

 // Flag for dumping every iteration
 CommandLineArgs::specify_command_line_flag("--dump_at_every_step");

 CommandLineArgs::specify_command_line_flag("--singular_tolerance", 
  &TestSoln::Singular_tolerance);

 CommandLineArgs::specify_command_line_flag("--singular_prefactor", 
  &TestSoln::Singular_prefactor);

 double newton_solver_tolerance = 1e-8;
 CommandLineArgs::specify_command_line_flag("--singular_tolerance", 
  &newton_solver_tolerance);

 CommandLineArgs::specify_command_line_flag("--arnoldi_multiple", 
  &TestSoln::Arnoldi_multiple);
 
 // Parse command line
 CommandLineArgs::parse_and_assign(); 

 // Cast max depth to int
 if(CommandLineArgs::command_line_flag_has_been_set("--max_singular_steps"))
  {TestSoln::max_depth = int(max_depth);}

 // Save disk space flag
 TestSoln::save_disk_space =CommandLineArgs::command_line_flag_has_been_set("--save_disk_space");
 // Dump Triangulateio
 TestSoln::Dump_triangulateio =CommandLineArgs::command_line_flag_has_been_set("--dump_triangulateio");
 // Do h loop at const vk pressure 
 bool do_h_loop_const_vk_pressure = CommandLineArgs::command_line_flag_has_been_set("--keep_vk_pressure_constant");
 // Do fd jacobian
 TestSoln::use_fd_jacobian =CommandLineArgs::command_line_flag_has_been_set("--use_fd_jacobian");
 // Dump at every step flag
 TestSoln::dump_at_every_step=CommandLineArgs::command_line_flag_has_been_set("--dump_at_every_step");

 TestSoln::Use_jfm_nondimensionalisation=CommandLineArgs::command_line_flag_has_been_set("--use_jfm_nondim");
  
 // Only do first singular
 TestSoln::only_first_singular_value=CommandLineArgs::command_line_flag_has_been_set("--only_do_first_singular");

 // Set do eigen as false if no eigen set
 TestSoln::do_eigen=!(CommandLineArgs::command_line_flag_has_been_set("--no_eigen"));

 // Direct solver flag
 TestSoln::use_direct_solver =(CommandLineArgs::command_line_flag_has_been_set("--use_LAPACK"));

 // Check element area.
 if(element_area > 0.1)
  {
   oomph_info <<"Element area has been set to greater than 0.1, which causes \
problems for the curved edge elements. The prescribed area will be ignored."<<std::endl;
   element_area =0.1;
  }

 oomph_info<<"Setting Element area to :" <<element_area <<std::endl; 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Do a quick solve for initial guess
 TestSoln::p_mag=p_start;
 double first_p = p_start;

 // Problem instance
 UnstructuredFvKProblem<LargeDisplacementPlateC1CurvedBellElement<2,2,3,KoiterSteigmannPlateEquations> >problem(element_area);
 // Set pointer to the problem
 TestSoln::problem_pt = &problem;

 // Set up eigensolver parameters
 if(TestSoln::use_direct_solver) 
  {
   problem.eigen_solver_pt() = new LAPACK_QZ;
   TestSoln::neig = problem.ndof(); 
  }
 else
  {
   problem.eigen_solver_pt() = new ARPACK;
   // Downcast to ARPACK*
   ARPACK* arpack_pt=(dynamic_cast<ARPACK*>(problem.eigen_solver_pt()));
   arpack_pt->narnoldi() = TestSoln::Arnoldi_multiple * TestSoln::neig;
   // Look for LARGE eigenvalues in shift inverted system i.e smallest eigenvalues
   arpack_pt->get_eigenvalues_right_of_shift();
  }
 // Change some tolerances
 problem.max_residuals()=1e10;
 problem.max_newton_iterations()=30;
 problem.newton_solver_tolerance()=newton_solver_tolerance;

 // Check if a checkpoint file is present at start-up: use it if so 
 std::ifstream chk;
 chk.open("checkpoint.dump");
 if (chk.is_open()) {
  std::cerr<<"Restarting from checkpoint file."<<std::endl;
  problem.read(chk);    
  chk.close();

  std::cerr<<"Opening checkpoint_p."<<std::endl;
  chk.open("checkpoint_p");
  // Check the file opened
  if(!chk.is_open())
   {
    std::cerr<<"Failed to open checkpoint_p file."<<std::endl;
    std::cerr<<"Exiting."<<std::endl;
    exit(-1);
   }
  chk>>TestSoln::p_mag;
  oomph_info<<"The pressure read in is: "<<TestSoln::p_mag<<std::endl;
  chk.close();

  std::cerr<<"Opening checkpoint_nu."<<std::endl;
  chk.open("checkpoint_nu");
  // Check the file opened
  if(!chk.is_open())
   {
    std::cerr<<"Failed to open checkpoint_nu file."<<std::endl;
    std::cerr<<"Exiting."<<std::endl;
    exit(-1);
   }
  chk>>TestSoln::nu;
  oomph_info<<"The poisson ratio read in is: "<<TestSoln::nu<<std::endl;
  chk.close();

  std::cerr<<"Opening checkpoint h."<<std::endl;
  chk.open("checkpoint_h");
  // Check the file opened
  if(!chk.is_open())
   {
    std::cerr<<"Failed to open checkpoint_h file."<<std::endl;
    std::cerr<<"Exiting."<<std::endl;
    exit(-1);
   }
  oomph_info<<"The thickness read in is: "<<TestSoln::h<<std::endl;
  chk>>TestSoln::h;
  chk.close();

  std::cerr<<"Opening checkpoint doc info number."<<std::endl;
  chk.open("checkpoint_nsoln");
  // Check the file opened
  if(!chk.is_open())
   {
    std::cerr<<"Failed to open checkpoint_h file."<<std::endl;
    std::cerr<<"Exiting."<<std::endl;
    exit(-1);
   }
  chk>>problem.Doc_info.number();
  oomph_info<<"The info num read in is: "<<problem.Doc_info.number()<<std::endl;
  chk.close();

  // Now repick up the loop
  first_p =TestSoln::p_mag ;//i_start = n_step - unsigned((p_end-TestSoln::p_mag)/delta_p);
  oomph_info << "Restarting at pressure" << first_p <<std::endl;
  
  oomph_info << "Newton Solve to check"<<std::endl;
  problem.newton_solve();
 } 

//  problem.print_jacobian();

 // Loop to get to target pressure
 // problem.loop_until_target_pressure(first_p,p_end,(p_end-p_start)/(n_step));
 double h_start = TestSoln::h;
 if(n_step !=0)
  problem.loop_until_target_pressure(first_p,p_end,(p_end-p_start)/(n_step));

 if(n_h_step!=0 &&  !do_h_loop_const_vk_pressure)
  problem.loop_until_target_thickness(h_start,h_end,(h_end-h_start)/(n_h_step));

 else if(n_h_step!=0 && do_h_loop_const_vk_pressure)
  problem.loop_until_target_thickness_keep_fvk_pressure_constant(h_start,h_end
   ,(h_end-h_start)/(n_h_step));

 delete problem.eigen_solver_pt();
 remove("checkpoint.dump");
 remove("checkpoint_p");
 remove("checkpoint_h");
 remove("checkpoint_nu");
 remove("checkpoint_nsoln");

 oomph_info <<"Exiting Normally."<<std::endl;
} //End of main


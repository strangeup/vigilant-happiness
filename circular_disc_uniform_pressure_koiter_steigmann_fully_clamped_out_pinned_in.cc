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
double displ_outer = 0.0;

// Enum for displacement boundary condition type
enum Outer_stretch_type { 
 radial_stretch = 0, uniaxial_stretch = 1, random_mode_stretch = 3
 };

// Default mode
Outer_stretch_type mode = radial_stretch;


/*                     PARAMETRIC BOUNDARY DEFINITIONS                        */
// Here we create the geom objects for the Parametric Boundary Definition 
CurvilineCircleTop parametric_curve_top;
CurvilineCircleBottom parametric_curve_bottom;

// Use a Mathematica friendly alias for std::pow
template<typename T1, typename T2>
double Power(const T1& arg, const T2& exp)
{ return std::pow(arg,exp); }

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
{
}

// Assigns the value of pressure depending on the position (x,y)
void get_d_pressure_dr(const Vector<double>& xi,const Vector<double>& ui,
  const DenseMatrix<double>& dui_dxj,
 const Vector<double>& ni, DenseMatrix<double>& d_pressure_dn)
{
 
}

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

void get_my_prestretch_solution(const Vector<double>& xi, Vector<double>& w)
{
 const double x = xi[0], y = xi[1];
 
 w[0] = displ_outer*x; 
 w[1] = displ_outer; 
 w[6] = displ_outer*y; 
 w[8] = displ_outer; 
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

//Exact solution for constant pressure, circular domain and resting boundary conditions
void get_exact_w_radial(const Vector<double>& xi, Vector<double>& w)
{}

// Boolean for linear elasticity
bool use_linear_elasticity = false;

// Dimensionless Mooney rivlin constant C2 = C1 /10 and C1 +C2 = 1/6
double C2 = 1./66.; 

// Mooney Rivlin stress function 
void mooney_rivlin_stress(const Vector<double>& x, const Vector<double>& u,
  const DenseMatrix<double>& e, const DenseMatrix<double>& g, 
  DenseMatrix<double>& stress)
{
 // Constants
 const double c1 = 1/6. - TestSoln::C2;
 const double c2 = TestSoln::C2;

 // Matrix of cofactors of g   
 DenseMatrix<double> cof_g(2,2);
 cof_g(0,0) = g(1,1); 
 cof_g(1,1) = g(0,0); 
 cof_g(0,1) =-g(0,1); 
 cof_g(1,0) =-g(1,0);
 DenseMatrix<double> cof_e(2,2);
 cof_e(0,0) = e(1,1); 
 cof_e(1,1) = e(0,0); 
 cof_e(0,1) =-e(0,1); 
 cof_e(1,0) =-e(1,0);
  
 // HERE STOP USING g and e interchangebaly - they may have different scales
 const double det_g = g(0,0)*g(1,1)-g(0,1)*g(1,0); 
 const double det_e = e(0,0)*e(1,1)-e(0,1)*e(1,0); 
 // const double tr_g = g(0,0)+g(1,1); 
 const double tr_e = e(0,0)+e(1,1);
 // NB det(g) = 4 det(e) + 2 Tr(e) +1  
 // Determinant of g squared minus one
 const double det_g_2_m1 = (4*det_e+2*tr_e)*(2+4*det_e+2*tr_e); 

 // Now fill in the stress
 // Loop over indices
 DenseMatrix<double> i2(2,2,0.0);
 i2(0,0)=1.0;
 i2(1,1)=1.0;

 // Now Fill in the Stress
 for(unsigned alpha=0;alpha<2;++alpha)
  {
  for(unsigned beta=0;beta<2;++beta)
   {
   //  // 2nd Piola Kirchhoff (Membrane) Stress for Mooney Rivlin model
   //  stress(alpha,beta)=2*(c1+ c2 / det_g) * i2(alpha,beta) 
   //       + 2*((- c1 - tr_g *c2) / pow(det_g,2) + c2)*cof_g(alpha,beta);
    // For 2D:
    // Cof g = I2 + 2 Cof e
    // tr(g) = 2 + 2 tr(e)
    // Det(g) = 4 Det(e) + 2 Tr(e) +1 
    // 2nd Piola Kirchhoff (Membrane) Stress for Mooney Rivlin model
    // ( Modified so that all terms are order epsilon if possible *)
    stress(alpha,beta)=2*(c2*(-4*det_e-2*tr_e) / det_g) * i2(alpha,beta) 
         - 4*((c1 + 2*c2 +  2*tr_e *c2) / pow(det_g,2) - c2)*cof_e(alpha,beta)
         - 2*((c1 + 2*c2 +  2*tr_e *c2)*(-det_g_2_m1)/pow(det_g,2)+2*tr_e*c2)*i2(alpha,beta);
   }
  }
}

// Mooney Rivlin stiffness tensor 
void d_mooney_rivlin_stress(const Vector<double>& x, const Vector<double>& u,
  const DenseMatrix<double>& strain, const DenseMatrix<double>& g, 
  RankThreeTensor<double>& d_stress_du, RankFourTensor<double>& d_stress_dstrain)
{
 // HERE STOP USING g and e interchangebaly - they may have different scales
 // Constants
 const double c1 = 1/6. - TestSoln::C2;
 const double c2 = TestSoln::C2;

 // Matrix of cofactors of g   
 DenseMatrix<double> cof_g(2,2);
 cof_g(0,0) = g(1,1); 
 cof_g(1,1) = g(0,0); 
 cof_g(0,1) =-g(0,1); 
 cof_g(1,0) =-g(1,0);
  
 // Fill in det
 const double det_g = g(0,0)*g(1,1)-g(0,1)*g(1,0); 
 const double tr_g = g(0,0)+g(1,1); 

 // Now fill in the stress
 // Loop over indices
 DenseMatrix<double> i2(2,2,0.0);
 i2(0,0)=1.0;
 i2(1,1)=1.0;

 // Now Fill in the Stress
 for(unsigned alpha=0;alpha<2;++alpha)
  {
  for(unsigned beta=0;beta<2;++beta)
   {
    // 2nd Piola Kirchhoff (Membrane) Stress for Mooney Rivlin model
   // stress(alpha,beta)=2*(c1+ c2 / det_g) * i2(alpha,beta) 
   //      + 2*((- c1 - tr_g *c2) / pow(det_g,2) + c2)*cof_g(alpha,beta);
   for(unsigned gamma=0;gamma<2;++gamma)
    {
    for(unsigned delta=0;delta<2;++delta)
     {
     // 2nd Piola Kirchhoff (Membrane) Stress for Mooney Rivlin model
     // d (cof_g) dg
     d_stress_dstrain(alpha,beta,gamma,delta)=2*( 
       - 2*(c2 / pow(det_g,2)) * i2(alpha,beta) * cof_g(gamma,delta) 
       + 2*(c2 - (c1 + tr_g *c2) / pow(det_g,2))*
       // dcof(g)/dg = (cof(g) \otimes cof(g) - cof(g) . dg / dg . cof(g))/det(g)
         (cof_g(alpha,beta)*cof_g(gamma,delta) - cof_g(alpha,gamma)*cof_g(delta,beta))
          / det_g
       + 4*((c1 + tr_g *c2) / pow(det_g,3))*cof_g(alpha,beta)*cof_g(gamma,delta)
       - 2*(c2 / pow(det_g,2))*cof_g(alpha,beta)*i2(gamma,delta));
     }
    }
   }
  }
}

// Random displacement from mathematica
void random_fourier_mode_displacement(const double& theta, Vector<double>& ui)
 {
 // Create aliases
 double (*Sin)(double theta)(& std::sin);
 double (*Cos)(double theta)(& std::cos);
 /**************************************************************/
 /* Runtime is:Wed 21 Feb 2018 11:09:00*/
 /*
 h = (3 + Cos[2*(b + 3*theta)] + Cos[2*(d + 3*theta)] + Cos[2*(a +  4*theta)])/3
 a = 0.09420264555394686
 b = 0.9530505658231787
 d = 0.1828786107259499
 */

 // Randomly picked phases
 const double a =  0.09420264555394686;
 
 const double b =  0.9530505658231787;
 
 const double d =  0.1828786107259499;
 
 // Function with randomly picked fourier modes for radial displacement
 // The average of this function over 2Pi is 1
 const double random_pert =(3 + Cos(2*(b + 3*theta)) + Cos(2*(d + 3*theta)) + \
 Cos(2*(a + 4*theta)))/3. ;
 const double d_random_pert =(-2*(3*Sin(2*(b + 3*theta)) + 3*Sin(2*(d + 3*theta)) + \
 4*Sin(2*(a + 4*theta))))/3. ;
 const double d2_random_pert =(-4*(9*Cos(2*(b + 3*theta)) + 9*Cos(2*(d + 3*theta)) \
 + 16*Cos(2*(a + 4*theta))))/3. ;
/**************************************************************/

 // Now
 ui[0] = displ_outer*(random_pert*Cos(theta));
 ui[2] = displ_outer*(d_random_pert*Cos(theta) - random_pert*Sin(theta));
 ui[5] = displ_outer*(-random_pert*Cos(theta) - 2*d_random_pert*Sin(theta)
    + d2_random_pert*Cos(theta));

 ui[6] =displ_outer*(random_pert*Sin(theta));
 ui[8] =displ_outer*(d_random_pert*Sin(theta) + random_pert*Cos(theta));
 ui[11]=displ_outer*(-random_pert*Sin(theta) + 2*d_random_pert*Cos(theta) 
    + d2_random_pert*Sin(theta));
 }

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
apply_boundary_conditions();
}

/// Doc the solution
void doc_solution(const std::string& comment="");

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

void set_values_to_applied_prestretch()
{
// Complete the build of all elements so they are fully functional
unsigned n_element = Bulk_mesh_pt->nelement();
for(unsigned e=0;e<n_element;e++)
{
 // Upcast from GeneralisedElement to the present element
 ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
 unsigned nnode = el_pt->nnode();
 for(unsigned i=0;i<nnode;++i)
  {
   // Containers
   Vector<double> x(2),w(18,0.0);
   // Get node pointer
   Node* nod_pt = el_pt->node_pt(i);
   // Get x 
   x[0]=nod_pt->x(0);
   x[1]=nod_pt->x(1);
   // Get test
   TestSoln::get_my_prestretch_solution(x,w);
   // Set value
  for(unsigned l=0 ; l<6;++l)
   {
   for(unsigned j=0 ; j<3;++j)
    {
    unsigned index=el_pt->u_index_koiter_model(l,j);
    nod_pt->set_value(index,w[index]);
    }
   }
  }
  // Pin the internal dofs
  const unsigned n_internal_dofs = el_pt->number_of_internal_dofs();
  for(unsigned i=0;i<n_internal_dofs;++i)
   {
   for(unsigned j=0;j<3;++j)
    {
    // Containers
    Vector<double> x(2,0.0), s(2,0.0), w(18,0.0);
    // Get internal data
    Data* internal_data_pt = el_pt->internal_data_pt(1);
    el_pt->get_internal_dofs_location(i,s);
    // Get test
    el_pt->get_coordinate_x(s,x); 
    TestSoln::get_my_prestretch_solution(x,w);
    // HERE need a function that can give us this lookup
     internal_data_pt->set_value(i+j*n_internal_dofs,w[6*j+0]);
    }
   } 
  apply_boundary_conditions(); //Just in case 
}

} // end set bc
private:

// Mesh parameters
Ellipse* Outer_boundary_ellipse_pt;
TriangleMeshClosedCurve* Outer_boundary_pt;
Vector<TriangleMeshCurveSection*> Outer_curvilinear_boundary_pt;
Vector<TriangleMeshOpenCurve *> Inner_open_boundaries_pt;
TriangleMeshPolyLine* Boundary2_pt;

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
 Inner_boundary0 = 2
};

double Element_area;

void upgrade_edge_elements_to_curve(const unsigned &b, Mesh* const & 
 bulk_mesh_pt);

void rotate_edge_degrees_of_freedom(Mesh* const &bulk_mesh_pt);

/// \short Set pointer to prescribed-flux function for all elements
/// in the surface mesh
void set_prescribed_traction_pt();

/// Pointer to "bulk" mesh
TriangleMesh<ELEMENT>* Bulk_mesh_pt;

/// Pointer to "surface" mesh
Mesh* Surface_mesh_pt;

}; // end_of_problem_class


template<class ELEMENT>
UnstructuredFvKProblem<ELEMENT>::UnstructuredFvKProblem(double element_area)
:
Element_area(element_area)
{
Vector<double> zeta(1);
Vector<double> posn(2);

//Outer boundary
//--------------

double A = 1.0;
double B = 1.0;
// Initialize
Outer_boundary_ellipse_pt = new Ellipse(A, B);
Outer_boundary_pt = 0;
Outer_curvilinear_boundary_pt.resize(2);
// We want internal open curves
// Internal open boundaries
// Total number of open curves in the domain
unsigned n_open_curves = 1;
Inner_open_boundaries_pt.resize(n_open_curves);
Boundary2_pt =  0;


//First bit
double zeta_start = 0.0;
double zeta_end = MathematicalConstants::Pi;
unsigned nsegment = (int)(MathematicalConstants::Pi/sqrt(element_area));
Outer_curvilinear_boundary_pt[0] =
new TriangleMeshCurviLine(Outer_boundary_ellipse_pt, zeta_start,
zeta_end, nsegment, Outer_boundary0);

//Second bit
zeta_start = MathematicalConstants::Pi;
zeta_end = 2.0*MathematicalConstants::Pi;
nsegment = (int)(MathematicalConstants::Pi/sqrt(element_area));
Outer_curvilinear_boundary_pt[1] =
new TriangleMeshCurviLine(Outer_boundary_ellipse_pt, zeta_start,
zeta_end, nsegment, Outer_boundary1);

Outer_boundary_pt =
new TriangleMeshClosedCurve(Outer_curvilinear_boundary_pt);


// Internal bit - this means we can have a boundary which is just the centre
// We start by creating the internal boundaries
// The boundary 2 is defined by its two vertices
// Open curve 1
 Vector<Vector<double> > vertices(2,Vector<double>(2,0.0));
 vertices[0][0] = 0.0;
 vertices[0][1] = 0.0;
 
 vertices[1][0] = 0.5;
 vertices[1][1] = 0.0;
 unsigned boundary_id = Inner_boundary0;
 Boundary2_pt =  new TriangleMeshPolyLine(vertices, boundary_id);

// Each internal open curve is defined by a vector of
// TriangleMeshCurveSection,
// on this example we only need one curve section for each internal boundary
 Vector<TriangleMeshCurveSection *> internal_curve_section1_pt(1);
 internal_curve_section1_pt[0] = Boundary2_pt;

// The open curve that define this boundary is composed of just one
// curve section
 Inner_open_boundaries_pt[0] =
    new TriangleMeshOpenCurve(internal_curve_section1_pt);

//Create the mesh
//---------------
//Create mesh parameters object
TriangleMeshParameters mesh_parameters(Outer_boundary_pt);

// Specify the internal open boundaries
mesh_parameters.internal_open_curves_pt() = Inner_open_boundaries_pt;
// Specify the element area
mesh_parameters.element_area() = element_area;

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

 // if Axisymmetric pin
 if(TestSoln::mode == TestSoln::radial_stretch)
   {
   // Set the radial displacement
   w[0]= TestSoln::displ_outer*x[0];
   w[2]=-TestSoln::displ_outer*x[1];
   w[5]=-TestSoln::displ_outer*x[0];
   w[6]= TestSoln::displ_outer*x[1];
   w[8]= TestSoln::displ_outer*x[0];
   w[11]=-TestSoln::displ_outer*x[1];
   }
 else  if(TestSoln::mode == TestSoln::uniaxial_stretch)
   {
   // Set the radial displacement
   w[0]= TestSoln::displ_outer*x[0];
   w[2]=-TestSoln::displ_outer*x[1];
   w[5]=-TestSoln::displ_outer*x[0];
   }
 else //  if(TestSoln::mode == TestSoln::uniaxial_stretch)
  {
   // The angle
   const double theta = atan2(x[1],x[0]);
   TestSoln::random_fourier_mode_displacement(theta,w); 
  }
 // Pin unknown values (everything except for the second normal derivative)
 for(unsigned i=2;i<3;++i)
  {
   for(unsigned l=0;l<6;++l)
    {
    // Don't pin the second normal (clamped)
    if(l != 3)  
     {
     unsigned index=el_pt->u_index_koiter_model(l,i);
     nod_pt->pin(index);
     nod_pt->set_value(index,w[index]);
     }
    }
  }
 // Pin unknown values (everything except for the second normal derivative)
 for(unsigned i=0;i<2;++i)
  {
   for(unsigned l=0;l<6;++l)
    {
    // Pin the value, tangent deriv. and second tangent deriv. (resting)
    if(l == 0 || l==2 || l == 5)
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
if(!TestSoln::use_linear_elasticity)
{
 el_pt->stress_fct_pt() = &TestSoln::mooney_rivlin_stress;
 el_pt->d_stress_fct_pt() = &TestSoln::d_mooney_rivlin_stress;
}
el_pt->nu_pt() = &TestSoln::nu;
el_pt->d_pressure_dn_fct_pt() = &TestSoln::get_d_pressure_dn;
el_pt->d_pressure_dr_fct_pt() = &TestSoln::get_d_pressure_dr;
el_pt->d_pressure_d_grad_u_fct_pt() = &TestSoln::get_d_pressure_d_grad_u;
el_pt->thickness_pt() = &TestSoln::h;

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
// Upcast to first element
ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(0));

for(unsigned ibound=0;ibound<nbound;ibound++)
{
unsigned num_nod=Bulk_mesh_pt->nboundary_node(ibound);
for (unsigned inod=0;inod<num_nod;inod++)
{
 // Get node
 Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
 
 // Extract nodal coordinates from node:
 Vector<double> x(2),w(18,0.0);
 x[0]=nod_pt->x(0);
 x[1]=nod_pt->x(1);

 // if Axisymmetric pin
 if(TestSoln::mode == TestSoln::radial_stretch)
   {
   // Set the radial displacement
   w[0]= TestSoln::displ_outer*x[0];
   w[2]=-TestSoln::displ_outer*x[1];
   w[5]=-TestSoln::displ_outer*x[0];
   w[6]= TestSoln::displ_outer*x[1];
   w[8]= TestSoln::displ_outer*x[0];
   w[11]=-TestSoln::displ_outer*x[1];
   }
 else  if(TestSoln::mode == TestSoln::uniaxial_stretch)
   {
   // Set the radial displacement
   w[0]= TestSoln::displ_outer*x[0];
   w[2]=-TestSoln::displ_outer*x[1];
   w[5]=-TestSoln::displ_outer*x[0];
   }
 else //  if(TestSoln::mode == TestSoln::uniaxial_stretch)
  {
   // The angle
   const double theta = atan2(x[1],x[0]);
   TestSoln::random_fourier_mode_displacement(theta,w); 
  }

 // Pin unknown values (everything except for the second normal derivative)
 for(unsigned i=2;i<3;++i)
  {
   for(unsigned l=0;l<6;++l)
    {
    // Don't pin the second normal (clamped)
    if(l != 3)  
     {
     unsigned index=el_pt->u_index_koiter_model(l,i);
     nod_pt->pin(index);
     nod_pt->set_value(index,w[index]);
     }
    }
  }
 // Pin unknown values (everything except for the second normal derivative)
 for(unsigned i=0;i<2;++i)
  {
   for(unsigned l=0;l<6;++l)
    {
    // Pin the value, tangent deriv. and second tangent deriv. (resting)
    if(l == 0 || l==2 || l == 5)
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
     if(nod_pt->is_on_boundary(0) || nod_pt->is_on_boundary(1))
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
 
   // Loop nodes  HERE
   unsigned nnode = 3;//nnode();
   unsigned nbnode=0 ;
   // Count the number of boundary nodes
   for (unsigned n=0; n<nnode;++n)
     {
      // Check it isn't on an internal boundary
      bool on_boundary_2=el_pt->node_pt(n)->is_on_boundary(2);
      bool on_boundary_3=el_pt->node_pt(n)->is_on_boundary(3);
      if(!(on_boundary_2 || on_boundary_3))
       {nbnode+=unsigned(el_pt->node_pt(n)->is_on_boundary());}
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
       bool on_boundary_2=el_pt->node_pt(n)->is_on_boundary(2);
       bool on_boundary_3=el_pt->node_pt(n)->is_on_boundary(3);
       if(!(on_boundary_2 || on_boundary_3))
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
}// end create rotate elements

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

// // Doc error and return of the square of the L2 error
// //---------------------------------------------------
// //double error,norm,dummy_error,zero_norm;
  double dummy_error,zero_norm;
 sprintf(filename,"RESLT/error%i-%f.dat",Doc_info.number(),Element_area);
 some_file.open(filename);
 
 Bulk_mesh_pt->compute_error(some_file,TestSoln::get_exact_w,
                        dummy_error,zero_norm);
 some_file.close();
 
 // Doc L2 error and norm of solution
 oomph_info << "Error of computed solution: " << /*sqrt*/(dummy_error)<< std::endl;
 oomph_info << "Norm of computed solution: "  << /*sqrt*/(zero_norm)  << std::endl;
 
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

// Only divide by norm if its nonzero
some_file.close();
}
// Increment the doc_info number
Doc_info.number()++;
} // end of doc


// Namespace extension
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

  sprintf(filename,"checkpoint_C2" );
  filestream.open(filename);
  filestream <<TestSoln::C2;
  filestream.close();

  sprintf(filename,"checkpoint_uouter" );
  filestream.open(filename);
  filestream <<TestSoln::displ_outer;
  filestream.close();

  sprintf(filename,"checkpoint_is_linear" );
  filestream.open(filename);
  filestream <<TestSoln::use_linear_elasticity;
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
 // Print the git version using the Makefile MACRO
 printf("Version: %s\n", GITVERSION);
 printf("Binary build date: %s @ %s\n", __DATE__, __TIME__);
 // Catch soft kill signal from Condor and checkpoint if we get evicted 
 signal(SIGTERM, TestSoln::signal_handler);

 feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
// Dodgy turn off all errors fedisableexcept(FE_ALL_EXCEPT); // somewhere early in main(), before the exception is happening 

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Define possible command line arguments and parse the ones that
 // were actually specified

 // Validation?
 CommandLineArgs::specify_command_line_flag("--validation");

 // Directory for solution
 string output_dir="RSLT";
 CommandLineArgs::specify_command_line_flag("--dir", &output_dir);

 // Poisson Ratio
 //CommandLineArgs::specify_command_line_flag("--nu", &TestSoln::nu);
 
 // Applied Pressure
 CommandLineArgs::specify_command_line_flag("--p", &TestSoln::p_mag);

 double displ_target = 0;
 TestSoln::displ_outer = 0;
 CommandLineArgs::specify_command_line_flag("--u_outer",&displ_target);

 CommandLineArgs::specify_command_line_flag("--do_uniaxial_outer_pin");

 CommandLineArgs::specify_command_line_flag("--do_modal_outer_pin");

 CommandLineArgs::specify_command_line_flag("--nu", &TestSoln::nu);

 CommandLineArgs::specify_command_line_flag("--h", &TestSoln::h);

 CommandLineArgs::specify_command_line_flag("--C2", &TestSoln::C2);

 // n_step
 unsigned n_step=100;
 CommandLineArgs::specify_command_line_flag("--n_step", &n_step);
 
 // Element Area
 double element_area=0.1;
 CommandLineArgs::specify_command_line_flag("--element_area", &element_area);
 // P_end
 double p_start=0;
 // E = 1.36 e 6
 // h/L = 4.25e-3
 // So p*= 800 => 0.2 * 1.34* 4.25* 1e3 = 800 
 CommandLineArgs::specify_command_line_flag("--p_start", &p_start);

 // Number if prestetch steps
 unsigned n_prestretch_step=0;
 CommandLineArgs::specify_command_line_flag("--n_prestretch_steps",&n_prestretch_step);

 // P_end
 double p_end=0.2;
 CommandLineArgs::specify_command_line_flag("--p_end", &p_end);

 // Flag for saving disk space
 CommandLineArgs::specify_command_line_flag("--save_disk_space");

 // Flag for dumping every iteration
 CommandLineArgs::specify_command_line_flag("--dump_at_every_step");

 CommandLineArgs::specify_command_line_flag("--dry_run");
 
 CommandLineArgs::specify_command_line_flag("--use_linear_elasticity");

 // Flag for using analytic initial guess
 CommandLineArgs::specify_command_line_flag("--use_analytic_radial_prestretch_as_initial_guess");

 // Parse command line
 CommandLineArgs::parse_and_assign(); 

 // Save disk space flag
 TestSoln::save_disk_space =CommandLineArgs::command_line_flag_has_been_set("--save_disk_space");
 TestSoln::use_linear_elasticity=CommandLineArgs::command_line_flag_has_been_set("--use_linear_elasticity");
 bool use_analytic_radial_applied_prestretch =
   CommandLineArgs::command_line_flag_has_been_set("--use_analytic_radial_prestretch_as_initial_guess");
 bool dry_run =  CommandLineArgs::command_line_flag_has_been_set("--dry_run");
 // Check the mode selection is sensible
 if(CommandLineArgs::command_line_flag_has_been_set("--do_uniaxial_outer_pin")
  &&  CommandLineArgs::command_line_flag_has_been_set("--do_modal_outer_pin"))
  {
   oomph_info <<"Contradictory outer pin options requested."<< std::endl;
   oomph_info <<"Please select only a single option for outer pin type."<< std::endl;
   oomph_info <<"Note that the outer pin defaults to constant radial displacement."<< std::endl;
   oomph_info <<"Exiting."<< std::endl;
   exit(0);
  }
  else if (CommandLineArgs::command_line_flag_has_been_set("--do_uniaxial_outer_pin"))
  {
   TestSoln::mode = TestSoln::uniaxial_stretch;
  } 
  else if(CommandLineArgs::command_line_flag_has_been_set("--do_modal_outer_pin"))
  {
   TestSoln::mode = TestSoln::random_mode_stretch;
  }
  else
  {
   oomph_info <<"Defaulting to constant radial displacement, with magnitude: "
              <<displ_target<< std::endl;
   TestSoln::mode = TestSoln::radial_stretch;
  } 

 // Dump at every step flag
 bool dump_at_every_step=CommandLineArgs::command_line_flag_has_been_set("--dump_at_every_step");
  
 // Check element area.
 if(element_area > 0.1)
  {
   oomph_info <<"Element area has been set to greater than 0.1, which causes \
problems for the curved edge elements. The prescribed area will be ignored."<<std::endl;
   element_area =0.1;
  }
 // Check element area.
 if(TestSoln::C2 >= 1./6. || TestSoln::C2 <0)
  {
   oomph_info <<"C2 is meant to be a physical parameter bewteen 0 and 1/6.! \
Fix this and rerun. Exiting Script."<<std::endl;
   exit(-1);
  }

 oomph_info<<"Setting Element area to :" <<element_area <<std::endl; 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Do a quick solve for initial guess
 TestSoln::p_mag=p_start;
 unsigned i_start = 0;

 // Problem instance
 // UnstructuredFvKProblem<KoiterSteigmannC1CurvedBellElement<2,2,5> >problem(element_area);
 UnstructuredFvKProblem<LargeDisplacementPlateC1CurvedBellElement<2,2,3,KoiterSteigmannPlateEquations> >problem(element_area);
 // Set pointer to the problem
 // Set pointer to the problem
 TestSoln::problem_pt = &problem;
 // Change some tolerances
 problem.max_residuals()=1e10;
 problem.max_newton_iterations()=30;
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
  chk.close();
  oomph_info<<"The pressure read in is: "<<TestSoln::p_mag<<std::endl;


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
  chk.close();
  oomph_info<<"The poisson ratio read in is: "<<TestSoln::nu<<std::endl;

  std::cerr<<"Opening checkpoint h."<<std::endl;
  chk.open("checkpoint_h");
  // Check the file opened
  if(!chk.is_open())
   {
    std::cerr<<"Failed to open checkpoint_h file."<<std::endl;
    std::cerr<<"Exiting."<<std::endl;
    exit(-1);
   }
  chk>>TestSoln::h;
  chk.close();
  oomph_info<<"The thickness read in is: "<<TestSoln::h<<std::endl;

  std::cerr<<"Opening checkpoint C2."<<std::endl;
  chk.open("checkpoint_C2");
  // Check the file opened
  if(!chk.is_open())
   {
    std::cerr<<"Failed to open checkpoint_h file."<<std::endl;
    std::cerr<<"Exiting."<<std::endl;
    exit(-1);
   }
  chk>>TestSoln::C2;
  chk.close();
  oomph_info<<"The C2 parameter read in is: "<<TestSoln::C2<<std::endl;

  std::cerr<<"Opening checkpoint is_linear."<<std::endl;
  chk.open("checkpoint_is_linear");
  // Check the file opened
  if(!chk.is_open())
   {
    std::cerr<<"Failed to open checkpoint_h file."<<std::endl;
    std::cerr<<"Exiting."<<std::endl;
    exit(-1);
   }
  chk>>TestSoln::use_linear_elasticity;
  chk.close();
  oomph_info<<"The linear elasticty bool read in is: "<<TestSoln::use_linear_elasticity<<std::endl;

  std::cerr<<"Opening checkpoint u_outer."<<std::endl;
  chk.open("checkpoint_uouter");
  // Check the file opened
  if(!chk.is_open())
   {
    std::cerr<<"Failed to open checkpoint_uouter file."<<std::endl;
    std::cerr<<"Exiting."<<std::endl;
    exit(-1);
   }
  chk>>TestSoln::displ_outer;
  chk.close();
  oomph_info<<"The u_outer read in is: "<<TestSoln::displ_outer<<std::endl;

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
  chk.close();
  oomph_info<<"The info num read in is: "<<problem.Doc_info.number()<<std::endl;

  // Now repick up the loop
  double delta_p = (p_end-p_start)/(n_step);
  i_start = n_step - unsigned((p_end-TestSoln::p_mag)/delta_p);
  oomph_info<<"Restarting at iteration "<<i_start<<std::endl;
 } 


 // DISPLACEMENT LOOP
 // Looping to get to the 
 if (n_prestretch_step > 0)
  {
  const double u_step = (displ_target) / (n_prestretch_step);
  // Loop until we get there
  const double tol = 1e-12;
  oomph_info<<"Start displ_outer=" << TestSoln::displ_outer<<"\n";
  oomph_info<<"Target displ_outer=" << displ_target<<"\n";
  oomph_info<<"Step size=" << u_step<<"\n";
  oomph_info<<"Tolerance=" << tol <<"\n";
  unsigned iter = 1;
  // Now loop
  while(fabs(TestSoln::displ_outer - displ_target)>tol)
   {
   oomph_info<<"Solving for displ_outer=" << TestSoln::displ_outer<<"\n";
   oomph_info<<"Iteration " << iter<<std::endl;
   //Just in case we need a restart
   problem.store_current_dof_values();

   try {
   oomph_info<<"Doing Newton Solve"<<endl;
   if(!dry_run)
    problem.newton_solve();
   
   }
   catch(...)
    {
     oomph_info<<"Error caught in newton solve"<<endl;
     // Dump the data 
     char filename[100];
     std::ofstream filestream;
     filestream.precision(15);
     sprintf(filename,"%s/kspe_circle_data%i-%f.dump",
             problem.Doc_info.directory().c_str(),
             problem.Doc_info.number(),
             TestSoln::p_mag
            );
     filestream.open(filename);
     problem.dump(filestream);
     filestream.close();
     exit(0); 
   }
  
   // Document
   oomph_info<<"Trying doc"<<endl;
   if(!dry_run)
    problem.doc_solution();
   oomph_info << std::endl;
   oomph_info << "---------------------------------------------" << std::endl;
   oomph_info << "Prestretch (" << TestSoln::displ_outer << ")" << std::endl;
   oomph_info << "Current n_step  (" << n_step << ")" << std::endl;
   oomph_info << "Poisson ratio (" << TestSoln::nu << ")" << std::endl;
   oomph_info << "Solution number (" <<problem.Doc_info.number()-1 << ")" << std::endl;
   oomph_info << "---------------------------------------------" << std::endl;
   oomph_info << std::endl;
   //Increment
   TestSoln::displ_outer += u_step; 
   ++iter; 
   }
  }
 // IF we are not doing any steps
 else
  {
   oomph_info <<"Setting u_outer to u_target."<<std::endl;
   // Set the displ_outer to the target value
   TestSoln::displ_outer =  displ_target;
   // Use analytic radial prestretch as an initial guess
   if(use_analytic_radial_applied_prestretch)
    { problem.set_values_to_applied_prestretch();}
  }

 // PRESSURE LOOP
 // Set prev as the current mag
 TestSoln::p_prev = TestSoln::p_mag;
 oomph_info<<"Set p_prev to be p_mag"<<std::endl;
 // Loop Curvatures
 for(unsigned i=i_start;i<n_step+1;++i)
  {
  //Increment
  oomph_info<<"Starting iteration "<< i <<std::endl;
  oomph_info<<"Solving for pressure=" << TestSoln::p_mag<<"\n";
  //Just in case we need a restart
  problem.store_current_dof_values();

  try {
  oomph_info<<"Doing Newton Solve"<<endl;
  if(!dry_run)
   problem.newton_solve();
  
  }
  catch(...)
   {
    oomph_info<<"Error caught in newton solve"<<endl;
    // Dump the data 
    char filename[100];
    std::ofstream filestream;
    filestream.precision(15);
    sprintf(filename,"%s/kspe_circle_data%i-%f.dump",
            problem.Doc_info.directory().c_str(),
            problem.Doc_info.number(),
            TestSoln::p_mag
           );
    filestream.open(filename);
    problem.dump(filestream);
    filestream.close();
    exit(0); 
  }
 
  // If we are dumping at every step 
  if(dump_at_every_step && !dry_run)
  {
  oomph_info<<"Trying dump"<<endl;
   // Dump the data 
   char filename[100];
   std::ofstream filestream;
   filestream.precision(15);
   sprintf(filename,"%s/kspe_circle_data%i-%f.dump",
           problem.Doc_info.directory().c_str(),
           problem.Doc_info.number(),
           TestSoln::p_mag
          );
   filestream.open(filename);
   problem.dump(filestream);
   filestream.close();
  }
  
  // Document
  oomph_info<<"Trying doc"<<endl;
  if(!dry_run)
   problem.doc_solution();
  oomph_info << std::endl;
  oomph_info << "---------------------------------------------" << std::endl;
  oomph_info << "Pressure (" << TestSoln::p_mag << ")" << std::endl;
  oomph_info << "Prestretch (" << TestSoln::displ_outer << ")" << std::endl;
  oomph_info << "Current n_step  (" << n_step << ")" << std::endl;
  oomph_info << "Poisson ratio (" << TestSoln::nu << ")" << std::endl;
  oomph_info << "Solution number (" <<problem.Doc_info.number()-1 << ")" << std::endl;
  oomph_info << "---------------------------------------------" << std::endl;
  oomph_info << std::endl;
  // Backup previous pressure  
  TestSoln::p_prev = TestSoln::p_mag;
  TestSoln::p_mag+=(p_end-p_start)/(n_step);
   
}
remove("checkpoint.dump");
remove("checkpoint_p");
remove("checkpoint_h");
remove("checkpoint_nu");
remove("checkpoint_C2");
remove("checkpoint_is_linear");
remove("checkpoint_nsoln");
oomph_info <<"COMPLETED SUCCESSFULLY."<<std::endl;
} //End of main


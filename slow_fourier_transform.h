#ifndef SLOW_FOURIER_TRANSFORM
#define SLOW_FOURIER_TRANSFORM

#include <boost/math/special_functions/binomial.hpp>
namespace oomph {
namespace SlowFourierTransform {
// The smart way
inline double Slow_fourier_sine_function(const double& x, const double& y, const unsigned& n)
{
 // The idea is that if we are integrating over the whole domain we can Fourier
 // transform and then average over r using r^k Cos[k\theta] (which is easily
 // representable) as our Fourier operator
 //
 // So \int_\Omega \{ f(r,\theta) r^k Cos[k\theta]\} r dr d\theta
 // which = \int_\Omega \{ g(x,y) p(x,y)\} dx dy
 double sine_function = 0.0;

 // Loop over odd k 
 for(int l=0,k=n-1; k>=0; ++l,k=n - (2*l+1)) 
  {
   // +ve for odd l -ve for even l
   sine_function += std::pow(x,k)*std::pow(y,n-k)*
    boost::math::binomial_coefficient<double>(n,k)*(l%2==0 ?1:-1); 
  }
 return sine_function; 
}

// The smart way
inline double Slow_fourier_cosine_function(const double& x, const double& y, const unsigned& n)
{
 // The idea is that if we are integrating over the whole domain we can Fourier
 // transform and then average over r using r^k Cos[k\theta] (which is easily
 // representable) as our Fourier operator
 //
 // So \int_\Omega \{ f(r,\theta) r^k Cos[k\theta]\} r dr d\theta
 // which = \int_\Omega \{ g(x,y) p(x,y)\} dx dy
 double cosine_function =0;

 // Loop over even k 
 for(int l=0,k=n; k>=0; ++l,k=n - (2*l)) 
  {
   // +ve for odd l -ve for even l
   cosine_function += std::pow(x,k)*std::pow(y,n-k)*
    boost::math::binomial_coefficient<double>(n,k)*(l%2==0 ?1:-1); 
  }
 return cosine_function; 
} 
// The Cartesian definitions of r^k Cos[k \theta]
// For k is 1 to 30
inline double R_pow_2k_cos_ktheta(const double& x, const double& y, const unsigned& k)
{
 double (*Power)(double base, int exponent);
 double (*Sqrt)(double theta);
 Power = & std::pow;
 Sqrt = & std::sqrt;

 // Switch statement - hopefully will be optimised out by compiler
 switch(k)
  { 
 case 1 : 
  return x*Sqrt(Power(x,2) + Power(y,2));
 break;
 case 2 : 
  return Power(x,4) - Power(y,4);
 break;
 case 3 : 
  return Power(Power(x,2) + Power(y,2),1.5)*(Power(x,3) - 3*x*Power(y,2));
 break;
 case 4 : 
  return Power(Power(x,2) + Power(y,2),2)*
     (Power(x,4) - 6*Power(x,2)*Power(y,2) + Power(y,4));
 break;
 case 5 : 
  return Power(Power(x,2) + Power(y,2),2.5)*
     (Power(x,5) - 10*Power(x,3)*Power(y,2) + 5*x*Power(y,4));
 break;
 case 6 : 
  return Power(Power(x,2) + Power(y,2),3)*
     (Power(x,6) - 15*Power(x,4)*Power(y,2) + 15*Power(x,2)*Power(y,4) - 
       Power(y,6));
 break;
 case 7 : 
  return Power(Power(x,2) + Power(y,2),3.5)*
     (Power(x,7) - 21*Power(x,5)*Power(y,2) + 35*Power(x,3)*Power(y,4) - 
       7*x*Power(y,6));
 break;
 case 8 : 
  return Power(Power(x,2) + Power(y,2),4)*
     (Power(x,8) - 28*Power(x,6)*Power(y,2) + 70*Power(x,4)*Power(y,4) - 
       28*Power(x,2)*Power(y,6) + Power(y,8));
 break;
 case 9 : 
  return Power(Power(x,2) + Power(y,2),4.5)*
     (Power(x,9) - 36*Power(x,7)*Power(y,2) + 126*Power(x,5)*Power(y,4) - 
       84*Power(x,3)*Power(y,6) + 9*x*Power(y,8));
 break;
 case 10 : 
  return Power(Power(x,2) + Power(y,2),5)*
     (Power(x,10) - 45*Power(x,8)*Power(y,2) + 210*Power(x,6)*Power(y,4) - 
       210*Power(x,4)*Power(y,6) + 45*Power(x,2)*Power(y,8) - Power(y,10));
 break;
 case 11 : 
  return Power(Power(x,2) + Power(y,2),5.5)*
     (Power(x,11) - 55*Power(x,9)*Power(y,2) + 330*Power(x,7)*Power(y,4) - 
       462*Power(x,5)*Power(y,6) + 165*Power(x,3)*Power(y,8) - 11*x*Power(y,10));
 break;
 case 12 : 
  return Power(Power(x,2) + Power(y,2),6)*
     (Power(x,12) - 66*Power(x,10)*Power(y,2) + 495*Power(x,8)*Power(y,4) - 
       924*Power(x,6)*Power(y,6) + 495*Power(x,4)*Power(y,8) - 
       66*Power(x,2)*Power(y,10) + Power(y,12));
 break;
 case 13 : 
  return Power(Power(x,2) + Power(y,2),6.5)*
     (Power(x,13) - 78*Power(x,11)*Power(y,2) + 715*Power(x,9)*Power(y,4) - 
       1716*Power(x,7)*Power(y,6) + 1287*Power(x,5)*Power(y,8) - 
       286*Power(x,3)*Power(y,10) + 13*x*Power(y,12));
 break;
 case 14 : 
  return Power(Power(x,2) + Power(y,2),7)*
     (Power(x,14) - 91*Power(x,12)*Power(y,2) + 1001*Power(x,10)*Power(y,4) - 
       3003*Power(x,8)*Power(y,6) + 3003*Power(x,6)*Power(y,8) - 
       1001*Power(x,4)*Power(y,10) + 91*Power(x,2)*Power(y,12) - Power(y,14));
 break;
 case 15 : 
  return Power(Power(x,2) + Power(y,2),7.5)*
     (Power(x,15) - 105*Power(x,13)*Power(y,2) + 1365*Power(x,11)*Power(y,4) - 
       5005*Power(x,9)*Power(y,6) + 6435*Power(x,7)*Power(y,8) - 
       3003*Power(x,5)*Power(y,10) + 455*Power(x,3)*Power(y,12) - 
       15*x*Power(y,14));
 break;
 case 16 :
  return Power(Power(x,2) + Power(y,2),8)*
     (Power(x,16) - 120*Power(x,14)*Power(y,2) + 1820*Power(x,12)*Power(y,4) - 
       8008*Power(x,10)*Power(y,6) + 12870*Power(x,8)*Power(y,8) - 
       8008*Power(x,6)*Power(y,10) + 1820*Power(x,4)*Power(y,12) - 
       120*Power(x,2)*Power(y,14) + Power(y,16));
 break;
 case 17 :
  return Power(Power(x,2) + Power(y,2),8.5)*
     (Power(x,17) - 136*Power(x,15)*Power(y,2) + 2380*Power(x,13)*Power(y,4) - 
       12376*Power(x,11)*Power(y,6) + 24310*Power(x,9)*Power(y,8) - 
       19448*Power(x,7)*Power(y,10) + 6188*Power(x,5)*Power(y,12) - 
       680*Power(x,3)*Power(y,14) + 17*x*Power(y,16));
 break;
 case 18 : 
  return Power(Power(x,2) + Power(y,2),9)*
     (Power(x,18) - 153*Power(x,16)*Power(y,2) + 3060*Power(x,14)*Power(y,4) - 
       18564*Power(x,12)*Power(y,6) + 43758*Power(x,10)*Power(y,8) - 
       43758*Power(x,8)*Power(y,10) + 18564*Power(x,6)*Power(y,12) - 
       3060*Power(x,4)*Power(y,14) + 153*Power(x,2)*Power(y,16) - Power(y,18));
 break;
 case 19 :
  return Power(Power(x,2) + Power(y,2),9.5)*
     (Power(x,19) - 171*Power(x,17)*Power(y,2) + 3876*Power(x,15)*Power(y,4) - 
       27132*Power(x,13)*Power(y,6) + 75582*Power(x,11)*Power(y,8) - 
       92378*Power(x,9)*Power(y,10) + 50388*Power(x,7)*Power(y,12) - 
       11628*Power(x,5)*Power(y,14) + 969*Power(x,3)*Power(y,16) - 
       19*x*Power(y,18));
 break;
 case 20 : 
  return Power(Power(x,2) + Power(y,2),10)*
     (Power(x,20) - 190*Power(x,18)*Power(y,2) + 4845*Power(x,16)*Power(y,4) - 
       38760*Power(x,14)*Power(y,6) + 125970*Power(x,12)*Power(y,8) - 
       184756*Power(x,10)*Power(y,10) + 125970*Power(x,8)*Power(y,12) - 
       38760*Power(x,6)*Power(y,14) + 4845*Power(x,4)*Power(y,16) - 
       190*Power(x,2)*Power(y,18) + Power(y,20));
 break;
 case 21 : 
  return Power(Power(x,2) + Power(y,2),10.5)*
     (Power(x,21) - 210*Power(x,19)*Power(y,2) + 5985*Power(x,17)*Power(y,4) - 
       54264*Power(x,15)*Power(y,6) + 203490*Power(x,13)*Power(y,8) - 
       352716*Power(x,11)*Power(y,10) + 293930*Power(x,9)*Power(y,12) - 
       116280*Power(x,7)*Power(y,14) + 20349*Power(x,5)*Power(y,16) - 
       1330*Power(x,3)*Power(y,18) + 21*x*Power(y,20));
 break;
 case 22 : 
  return Power(Power(x,2) + Power(y,2),11)*
     (Power(x,22) - 231*Power(x,20)*Power(y,2) + 7315*Power(x,18)*Power(y,4) - 
       74613*Power(x,16)*Power(y,6) + 319770*Power(x,14)*Power(y,8) - 
       646646*Power(x,12)*Power(y,10) + 646646*Power(x,10)*Power(y,12) - 
       319770*Power(x,8)*Power(y,14) + 74613*Power(x,6)*Power(y,16) - 
       7315*Power(x,4)*Power(y,18) + 231*Power(x,2)*Power(y,20) - Power(y,22));
 break;
 case 23 : 
  return Power(Power(x,2) + Power(y,2),11.5)*
     (Power(x,23) - 253*Power(x,21)*Power(y,2) + 8855*Power(x,19)*Power(y,4) - 
       100947*Power(x,17)*Power(y,6) + 490314*Power(x,15)*Power(y,8) - 
       1144066*Power(x,13)*Power(y,10) + 1352078*Power(x,11)*Power(y,12) - 
       817190*Power(x,9)*Power(y,14) + 245157*Power(x,7)*Power(y,16) - 
       33649*Power(x,5)*Power(y,18) + 1771*Power(x,3)*Power(y,20) - 
       23*x*Power(y,22));
 break;
 case 24 : 
  return Power(Power(x,2) + Power(y,2),12)*
     (Power(x,24) - 276*Power(x,22)*Power(y,2) + 10626*Power(x,20)*Power(y,4) - 
       134596*Power(x,18)*Power(y,6) + 735471*Power(x,16)*Power(y,8) - 
       1961256*Power(x,14)*Power(y,10) + 2704156*Power(x,12)*Power(y,12) - 
       1961256*Power(x,10)*Power(y,14) + 735471*Power(x,8)*Power(y,16) - 
       134596*Power(x,6)*Power(y,18) + 10626*Power(x,4)*Power(y,20) - 
       276*Power(x,2)*Power(y,22) + Power(y,24));
 break;
 case 25 : 
  return Power(Power(x,2) + Power(y,2),12.5)*
     (Power(x,25) - 300*Power(x,23)*Power(y,2) + 12650*Power(x,21)*Power(y,4) - 
       177100*Power(x,19)*Power(y,6) + 1081575*Power(x,17)*Power(y,8) - 
       3268760*Power(x,15)*Power(y,10) + 5200300*Power(x,13)*Power(y,12) - 
       4457400*Power(x,11)*Power(y,14) + 2042975*Power(x,9)*Power(y,16) - 
       480700*Power(x,7)*Power(y,18) + 53130*Power(x,5)*Power(y,20) - 
       2300*Power(x,3)*Power(y,22) + 25*x*Power(y,24));
 break;
 case 26 : 
  return Power(Power(x,2) + Power(y,2),13)*
     (Power(x,26) - 325*Power(x,24)*Power(y,2) + 14950*Power(x,22)*Power(y,4) - 
       230230*Power(x,20)*Power(y,6) + 1562275*Power(x,18)*Power(y,8) - 
       5311735*Power(x,16)*Power(y,10) + 9657700*Power(x,14)*Power(y,12) - 
       9657700*Power(x,12)*Power(y,14) + 5311735*Power(x,10)*Power(y,16) - 
       1562275*Power(x,8)*Power(y,18) + 230230*Power(x,6)*Power(y,20) - 
       14950*Power(x,4)*Power(y,22) + 325*Power(x,2)*Power(y,24) - Power(y,26));
 break;
 case 27 : 
  return Power(Power(x,2) + Power(y,2),13.5)*
     (Power(x,27) - 351*Power(x,25)*Power(y,2) + 17550*Power(x,23)*Power(y,4) - 
       296010*Power(x,21)*Power(y,6) + 2220075*Power(x,19)*Power(y,8) - 
       8436285*Power(x,17)*Power(y,10) + 17383860*Power(x,15)*Power(y,12) - 
       20058300*Power(x,13)*Power(y,14) + 13037895*Power(x,11)*Power(y,16) - 
       4686825*Power(x,9)*Power(y,18) + 888030*Power(x,7)*Power(y,20) - 
       80730*Power(x,5)*Power(y,22) + 2925*Power(x,3)*Power(y,24) - 
       27*x*Power(y,26));
 break;
 case 28 : 
  return Power(Power(x,2) + Power(y,2),14)*
     (Power(x,28) - 378*Power(x,26)*Power(y,2) + 20475*Power(x,24)*Power(y,4) - 
       376740*Power(x,22)*Power(y,6) + 3108105*Power(x,20)*Power(y,8) - 
       13123110*Power(x,18)*Power(y,10) + 30421755*Power(x,16)*Power(y,12) - 
       40116600*Power(x,14)*Power(y,14) + 30421755*Power(x,12)*Power(y,16) - 
       13123110*Power(x,10)*Power(y,18) + 3108105*Power(x,8)*Power(y,20) - 
       376740*Power(x,6)*Power(y,22) + 20475*Power(x,4)*Power(y,24) - 
       378*Power(x,2)*Power(y,26) + Power(y,28));
 break;
 case 29 : 
  return Power(Power(x,2) + Power(y,2),14.5)*
     (Power(x,29) - 406*Power(x,27)*Power(y,2) + 23751*Power(x,25)*Power(y,4) - 
       475020*Power(x,23)*Power(y,6) + 4292145*Power(x,21)*Power(y,8) - 
       20030010*Power(x,19)*Power(y,10) + 51895935*Power(x,17)*Power(y,12) - 
       77558760*Power(x,15)*Power(y,14) + 67863915*Power(x,13)*Power(y,16) - 
       34597290*Power(x,11)*Power(y,18) + 10015005*Power(x,9)*Power(y,20) - 
       1560780*Power(x,7)*Power(y,22) + 118755*Power(x,5)*Power(y,24) - 
       3654*Power(x,3)*Power(y,26) + 29*x*Power(y,28));
 break;
 case 30 : 
  return Power(Power(x,2) + Power(y,2),15)*
     (Power(x,30) - 435*Power(x,28)*Power(y,2) + 27405*Power(x,26)*Power(y,4) - 
       593775*Power(x,24)*Power(y,6) + 5852925*Power(x,22)*Power(y,8) - 
       30045015*Power(x,20)*Power(y,10) + 86493225*Power(x,18)*Power(y,12) - 
       145422675*Power(x,16)*Power(y,14) + 145422675*Power(x,14)*Power(y,16) - 
       86493225*Power(x,12)*Power(y,18) + 30045015*Power(x,10)*Power(y,20) - 
       5852925*Power(x,8)*Power(y,22) + 593775*Power(x,6)*Power(y,24) - 
       27405*Power(x,4)*Power(y,26) + 435*Power(x,2)*Power(y,28) - Power(y,30));
 break;
 default :
    throw OomphLibError( "Cannot access slow fourier function > 30\
it is not defined.", OOMPH_CURRENT_FUNCTION,  
     OOMPH_EXCEPTION_LOCATION);
 break;
  }  
}
}

}
#endif

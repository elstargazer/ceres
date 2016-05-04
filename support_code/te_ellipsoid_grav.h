//============================================================================
// Name        : GravityTE.cpp
// Author      : Anton Ermakov
// Version     :
// Copyright   : Your copyright notice
// Description : Computing gravitational acceleration of a triaxial ellipsoid
//============================================================================

#include <iostream>
#include <time.h>

#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/exceptions.h>

#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>

#include "../support_code/local_math.h"

using namespace std;
using namespace dealii;

template<int dim>
double GetK(Point<dim, double> pt, vector<double> ell)
// This function solves the equation:
// x^2/(a^2+k) + y^2/(b^2+k) + z^2/(c^2+k) == 1
{
	double x, y, z, a, b, c, k;
	if (dim == 3)
	{
	    x = pt(0);
	    y = pt(1);
	    z = pt(2);

	    a = ell[0];
	    b = ell[1];
	    c = ell[2];
	}
	else if (dim ==2)
	{
	    x = pt(0);
	    y = 100.0;
	    z = pt(1);

	    a = ell[0];
	    b = ell[0];
	    c = ell[1];
	}

	double eps = 1e-12;
	if ((x*x/(a*a)+y*y/(b*b)+z*z/(c*c)) > (1+eps))
	{
		double k0;
		double r2 = x*x+y*y+z*z;

		k = r2 - a*a;

		double tol = 1e-14;
		double d = tol+1;

		while (d>tol)
		{
			k0 = k;
			k = x*x/(1+a*a/k) + y*y/(1+b*b/k) + z*z/(1+c*c/k);
			d = fabs((k-k0)/k);
		}
	}
	else
		k=0;

	return k;
}

template<int dim>
void GetDK(Point<dim, double> pt, vector<double> ell,
		double* DkDx, double* DkDy, double* DkDz)
// This function find the derivatives of k w.r.t x, y and z,
// where k is the solution of x^2/(a^2+k) + y^2/(b^2+k) + z^2/(c^2+k) == 1
{
	double x, y, z, a, b, c, k;
//	template double GetK<int>(Point<dim, double> pt, vector<double> ell)); // explicit instantiation.

	if (dim == 3)
	{
	    x = pt(0);
	    y = pt(1);
	    z = pt(2);

	    a = ell[0];
	    b = ell[1];
	    c = ell[2];

		k = GetK(pt,ell);
	}
	else if (dim == 2)
	{
	    x = pt(0);
	    y = 100.0;
	    z = pt(1);

	    a = ell[0];
	    b = ell[0];
	    c = ell[1];

		k = GetK(pt,ell);
	}

	double bot = pow(x/(a*a+k),2.0) +
			     pow(y/(b*b+k),2.0) +
			     pow(z/(c*c+k),2.0);

	*DkDx = 2*x/(bot*(a*a+k));
	*DkDy = 2*y/(bot*(b*b+k));
	*DkDz = 2*z/(bot*(c*c+k));
}

template <unsigned int dim>
class TEAnalyticGravity
// Class for computing gravitational acceleration analytically for
// a triaxial ellipsoid
{
public:

	std::vector<double> a_vec;
	std::vector<double> b_vec;
	std::vector<double> c_vec;
	std::vector<double> rho_vec;

	unsigned int nlayers;

	void setup_vars(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> rho);
	void get_gravity_layer(const dealii::Point<dim,double> &p,  std::vector<double> &g, unsigned int nl);
	void get_gravity(const dealii::Point<dim,double> &p,  std::vector<double> &g);

};

template <unsigned int dim>
// Setting dimensions and density of the triaxial ellipsoid
void TEAnalyticGravity<dim>::setup_vars (std::vector<double> ai, std::vector<double> bi, std::vector<double> ci, std::vector<double> rhoi)
{
	a_vec    = ai;
	b_vec    = bi;
	c_vec    = ci;
	rho_vec = rhoi;

	nlayers = a_vec.size();
}

template <unsigned int dim>
// This function computes gravitational acceleration of the triaxial ellisoid
void TEAnalyticGravity<dim>::get_gravity_layer(const dealii::Point<dim> &p, std::vector<double> &g, unsigned int nl)
{
	    double a = a_vec[nl];
	    double b = b_vec[nl];
	    double c = c_vec[nl];

	    vector<double> ell;
	    ell.push_back(a);
	    ell.push_back(b);
	    ell.push_back(c);

	    // define constants
		double G = 6.67e-11;

		double ax, ay, az;
		double x, y, z;

		if (dim == 3)
		{
		    x = p(0);
		    y = p(1);
		    z = p(2);
		}
		else if (dim == 2)
		{
		    x = p(0);
		    y = 100.0;
		    z = p(1);
		}
		std::cout << "ell axes = " << ell[0] << " " << ell[1] << " " << ell[2] << std::endl;
		std::cout << "compute gravity at " << x << " " << y << " " << z << std::endl;

		// check wheather or not point p is inside or outside of the ellipsoid
		double w = x*x/a/a+y*y/b/b+z*z/c/c;

		// define frequently used expressions
		double amc2=a*a-c*c;
		double amb2=a*a-b*b;
		double bmc2=b*b-c*c;

		double eratio = (a*a-b*b)/(a*a-c*c);
		double abc = a*b*c;

		double eps = 1e-12;

		if (w < (1-eps))

		// if w < 1, point p is inside the ellipsoid
		{
			double psi=asin(sqrt((a*a-c*c)/(a*a)));

			// compute elliptic intgrals of the first and second kind
			double F =  boost::math::ellint_1(sqrt(eratio),psi);
			double E =  boost::math::ellint_2(sqrt(eratio),psi);

			// compute components of the gravity acceleration
			ax=(2*abc)/sqrt(amc2)*(-2*x/amb2)*F+
					2*abc/sqrt(amc2)*(2*x/amb2)*E;

			ay=(2*abc)/sqrt(amc2)*(+2*y/amb2)*F+
					2*abc/sqrt(amc2)*(-amc2*2*y/(amb2*bmc2))*E+
					2*abc/sqrt((a*a)*(b*b)*(c*c))*((c*c)*2*y/bmc2);

			az=2*a*b*c/sqrt(amc2)*(+2*z/bmc2)*E+
					2*abc/sqrt((a*a)*(b*b)*(c*c))*(-(b*b)*2*z/bmc2);

			// assemble gravity vector
			g.clear();
		    if (dim == 2)
		    {
			    g.push_back(ax*PI*G);
			    g.push_back(az*PI*G);
		    }
		    else if (dim == 3)
		    {
			    g.push_back(ax*PI*G);
				g.push_back(ay*PI*G);
				g.push_back(az*PI*G);
		    }
		}
		else
		// if w > 1, the point p is outside of the ellipsoid
		{
			// compute k
//			template double GetK<dim>(Point<dim, double> pt, vector<double> ell);
			double k=GetK(p,ell);

			double a2pk = a*a+k;
			double b2pk = b*b+k;
			double c2pk = c*c+k;

			double psi=asin(sqrt((a*a-c*c)/a2pk));

			double DkDx;
			double DkDy;
			double DkDz;

			// compute k derivatives kx, ky, kz
//			template void GetDK<dim>(Point<dim, double> pt, vector<double> ell, double*, double*, double*);
			GetDK(p,ell,&DkDx,&DkDy,&DkDz);

			// compute elliptic integrals of the first and second kind
			double F =  boost::math::ellint_1(sqrt(eratio),psi);
			double E =  boost::math::ellint_2(sqrt(eratio),psi);

			// next bunch of line is to compute all terms in the derivation of the
			// potential expression
			double A = (2*abc)/sqrt(amc2)*(1-x*x/amb2+y*y/amb2);

			double tempA = 4*abc/(amb2*sqrt(amc2));

			std::cout << "vars1 = " << abc << " " << amc2 << " " << amb2 << " " << tempA << std::endl;

			double DADx = -tempA*x;
			double DADy =  tempA*y;

			double DFDpsi = pow((1 - eratio * pow(sin(psi),2)),-0.5);
			double DpsiDk = -0.5*sqrt(amc2/c2pk)/a2pk;

			double DFDk = DFDpsi * DpsiDk;

			double DFDx = DFDk * DkDx;
			double DFDy = DFDk * DkDy;
			double DFDz = DFDk * DkDz;

			double B = 2*abc/sqrt(amc2)*(x*x/amb2-amc2*y*y/(amb2*bmc2)+z*z/bmc2);

			double DBDx =  4*abc*x/(amb2*sqrt(amc2));
			double DBDy = -4*abc*sqrt(amc2)*y/(amb2*sqrt(bmc2));
			double DBDz =  4*abc*z/(sqrt(amc2)*bmc2);

			double DEDpsi = sqrt(1-eratio * pow(sin(psi),2));

			double DEDk = DEDpsi * DpsiDk;

			double DEDx = DEDk * DkDx;
			double DEDy = DEDk * DkDy;
			double DEDz = DEDk * DkDz;

		    double C = 2*abc/sqrt(a2pk*b2pk*c2pk);
		    double D = (c2pk*y*y/bmc2-b2pk*z*z/bmc2);

		    double DCDk = -abc*(a2pk*b2pk+a2pk*c2pk+b2pk*c2pk)/
		        pow(a2pk*b2pk*c2pk,1.5);

		    double DCDx = DCDk * DkDx;
		    double DCDy = DCDk * DkDy;
		    double DCDz = DCDk * DkDz;

		    double DDDk =  (y*y-z*z)/bmc2;
		    double DDDy =  2*c2pk*y/bmc2;
		    double DDDz = -2*b2pk*z/bmc2;

		    double DDDxtot = DDDk * DkDx;
		    double DDDytot = DDDk * DkDy + DDDy;
		    double DDDztot = DDDk * DkDz + DDDz;

		    // sum up all term to get gravitational acceleration

			double ax= DADx * F + DFDx * A + DBDx * E + B * DEDx + DCDx * D + C * DDDxtot;
			double ay= DADy * F + DFDy * A + DBDy * E + B * DEDy + DCDy * D + C * DDDytot;
		    double az= 0        + DFDz * A + DBDz * E + B * DEDz + DCDz * D + C * DDDztot;

			std::cout << "vars = "<< DADx << " " << F << " " << DFDx << " " << A << " " << DBDx << std::endl;

		    // assemble gravitational acceleration vector
		    g.clear();

		    if (dim == 2)
		    {
			    g.push_back(ax*PI*G);
			    g.push_back(az*PI*G);
		    }
		    else if (dim == 3)
		    {
			    g.push_back(ax*PI*G);
				g.push_back(ay*PI*G);
				g.push_back(az*PI*G);
		    }
		}
}

template<unsigned int dim>
void TEAnalyticGravity<dim>::get_gravity(const dealii::Point<dim> &p, std::vector<double> &g)
{
	std::vector<double> g_temp(0);
	for(unsigned int i=0; i<nlayers; i++)
	{
		get_gravity_layer(p, g_temp, i);
		if (i==0)
			for(unsigned int j=0; j<dim; j++)
			    g[j] += g_temp[j] * rho_vec[i];
		else
			for(unsigned int j=0; j<dim; j++)
		        g[j] += g_temp[j] * (rho_vec[i] - rho_vec[i-1]);
	}
}




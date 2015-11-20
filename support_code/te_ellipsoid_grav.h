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

using namespace std;
using namespace dealii;

double GetK(Point<3, double> pt, double a, double b, double c)
// This function solves the equation:
// x^2/(a^2+k) + y^2/(b^2+k) + z^2/(c^2+k) == 1
{
	double k;
	double x = pt(0);
	double y = pt(1);
	double z = pt(2);

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

void GetDK(Point<3, double> pt, double a, double b, double c,
		double* DkDx, double* DkDy, double* DkDz)
// This function find the derivatives of k w.r.t x, y and z,
// where k is the solution of x^2/(a^2+k) + y^2/(b^2+k) + z^2/(c^2+k) == 1
{
	double k = GetK(pt,a,b,c);
	double x = pt(0);
	double y = pt(1);
	double z = pt(2);

	double bot = pow(x/(a*a+k),2.0) +
			     pow(y/(b*b+k),2.0) +
			     pow(z/(c*c+k),2.0);

	*DkDx = 2*x/(bot*(a*a+k));
	*DkDy = 2*y/(bot*(b*b+k));
	*DkDz = 2*z/(bot*(c*c+k));
}

template <unsigned int dim>
class AnalyticGravity
// Class for computing gravitational acceleration analytically for
// a triaxial ellipsoid
{
public:
	double a, b, c, rho;

	void setup_vars(double a, double b, double c, double rho);
	void get_gravity(const dealii::Point<dim> &p,  std::vector<double> &g);

};

template <unsigned int dim>
// Setting dimensions and density of the triaxial ellipsoid
void AnalyticGravity<dim>::setup_vars (double ai, double bi, double ci, double rhoi)
{
	a=ai;
	b=bi;
	c=ci;
	rho=rhoi;
}

template <unsigned int dim>
// This function computes gravitational acceleration of the triaxial ellisoid
void AnalyticGravity<dim>::get_gravity (const dealii::Point<dim> &p, std::vector<double> &g)
{
	    Assert (dim == 3, ExcNotImplemented());

	    // define constants
		double G = 6.67e-11;
		double pi = 2.0 * asin(1.0);

		double ax, ay, az;
		double x, y, z;

		x = p(0);
		y = p(1);
		z = p(2);

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
			g.push_back(ax*pi*rho*G);
			g.push_back(ay*pi*rho*G);
			g.push_back(az*pi*rho*G);
		}
		else
		// if w > 1, the point p is outside of the ellipsoid
		{
			// compute k
			double k=GetK(p,a,b,c);

			double a2pk = a*a+k;
			double b2pk = b*b+k;
			double c2pk = c*c+k;

			double psi=asin(sqrt((a*a-c*c)/a2pk));

			double DkDx;
			double DkDy;
			double DkDz;

			// compute k derivatives kx, ky, kz
			GetDK(p,a,b,c,&DkDx,&DkDy,&DkDz);

			// compute elliptic integrals of the first and second kind
			double F =  boost::math::ellint_1(sqrt(eratio),psi);
			double E =  boost::math::ellint_2(sqrt(eratio),psi);

			// next bunch of line is to compute all terms in the derivation of the
			// potential expression
			double A = (2*abc)/sqrt(amc2)*(1-x*x/amb2+y*y/amb2);

			double tempA = 4*abc/(amb2*sqrt(amc2));

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

		    // assemble gravitational acceleration vector
			g.push_back(ax*pi*rho*G);
			g.push_back(ay*pi*rho*G);
			g.push_back(az*pi*rho*G);
		}
}

int main()
{
    time_t time_start = time (NULL);
    time_t time_finish;

	double a, b, c;
	double x, y, z;
	double rho;

	a = 3.0;
	b = 2.0;
	c = 1.0;

	int npts = 40.0;
	double fac = 444.;

	for(int i=1;i<npts;i++)
	{
		for(int j=1;j<npts;j++)
		{
			for(int k=1;k<npts;k++)
			{
				x = i/npts*fac;
				y = j/npts*fac;
				z = k/npts*fac;

				rho = 1000.0;

				Point<3,double> p = Point<3,double>(x,y,z);
				std::vector<double> g;

				AnalyticGravity<3> * aGrav = new AnalyticGravity<3>;  //sets up pointer to gravity object
				aGrav->setup_vars(a,b,c,rho);
				aGrav->get_gravity(p,g);

//				printf("ax = %23.16E\n", g[0]);
//				printf("ay = %23.16E\n", g[1]);
//				printf("az = %23.16E\n", g[2]);
			}
		}
	}

//	printf("ax = %23.16E\n", g[0]);
//	printf("ay = %23.16E\n", g[1]);
//	printf("az = %23.16E\n", g[2]);
	time_finish = time (NULL);
	printf("Time elapsed = %ld [sec]\n",time_finish-time_start);

	return 0;
}



#ifndef EXAMPLES_H
#define EXAMPLES_H

// Library of example routines used to test the operation of the root-finder class
// R. Sheehan 8 - 3 - 2013

namespace examples{

	// test functions
	double f1(double x);
	double f2(double x);
	double f3(double x);

	static const double lambda = 1.55; // wavelength
	static const double W = 2.5; // waveguide width
	static const double nc = 3.38; // core refractive index
	static const double ns = 3.17; // substrate refractive index
	static const double ncl = 1.0; // cladding refractive index
	static const double p=(atan(1.0)); // pi/4
	static const double Two_PI=(8.0*p); // 2*pi
	static const double PI=(4.0*p); // pi = 3.1415926785......
	static const double k0 = ((2.0*PI)/lambda); // wavenumber in slab waveguide
	static const double kncsqr=template_funcs::DSQR(k0*nc);
	static const double knssqr=template_funcs::DSQR(k0*ns);
	static const double knclsqr=template_funcs::DSQR(k0*ncl);

	// real world example
	// find propagation constants of slab waveguide by searching for roots of dispersion equation
	double disp_eqn(double x); 

	void ex_1(); 
	void ex_2(); 
	void ex_3(); 
	void ex_4(); 

}

#endif
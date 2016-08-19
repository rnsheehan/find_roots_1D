#ifndef ATTACH_H
#include "Attach.h"
#endif

double examples::f1(double x)
{
	// This function should have a root on [-4, 4]

	return exp(x) + pow(2.0,-x) + 2.0*sin(x) - 6.0; 
}

double examples::f2(double x)
{
	// This function should have a root on [0, 4]

	return 2.0*x*cos(2.0*x) - template_funcs::DSQR(x - 2.0); 
}

double examples::f3(double x)
{
	// This function should have a root on [-2, 5]

	return exp(x) - 3.0*template_funcs::DSQR(x); 
}

double examples::disp_eqn(double x)
{
	// The dispersion equation for the slab waveguide 
	// propagation constants will be found in the range [(k0*ns),(k0*nc)]

	double h, p, q; 
	double t1, t2, t3, t4, xsqr; 

	xsqr = template_funcs::DSQR(x); 
	h = sqrt(kncsqr - xsqr); 
	p = sqrt(xsqr - knssqr); 
	q = sqrt(xsqr - knclsqr); 

	t1 = (p*q) / template_funcs::DSQR(h); 
	t2 = q/h;
	t3 = p/h; 
	t4 = W*h; 

	return ((1.0-t1)*sin(t4)-(t2+t3)*cos(t4)); 
}

void examples::ex_1()
{
	// Find the roots of function f1

	find_root rtsearch(f1,-4.0,4.0,20); 

	rtsearch.bracket_roots(); 
	rtsearch.bisection_search(1.0e-5); 
	rtsearch.newton_raphson_search(1.0e-5);
}

void examples::ex_2()
{
	// Find the roots of function f1

	find_root rtsearch2(f2,0,4.0,20); 

	rtsearch2.bracket_roots(); 
	rtsearch2.bisection_search(1.0e-5); 
	rtsearch2.newton_raphson_search(1.0e-5);

	/*rtsearch2.set_interval(0.0,5.0,30); 
	rtsearch2.bracket_roots(); 
	rtsearch2.bisection_search(1.0e-5); 
	rtsearch2.newton_raphson_search(1.0e-5);*/
}

void examples::ex_3()
{
	// Find the roots of function f3

	find_root rtsearch3(f3,-2.0,5.0,20); 

	rtsearch3.bracket_roots(); 
	rtsearch3.bisection_search(1.0e-5);
	rtsearch3.newton_raphson_search(1.0e-5);
}

void examples::ex_4()
{
	// Find the roots of function disp_eqn

	find_root rtsearch4(disp_eqn,(k0*ns),(k0*nc),50); 

	rtsearch4.bracket_roots(); 
	rtsearch4.bisection_search(1.0e-5);
	rtsearch4.newton_raphson_search(1.0e-5);
}

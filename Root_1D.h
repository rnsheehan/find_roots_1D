#ifndef ROOT_1D_H
#define ROOT_1D_H

// Declare a class that can be used to determine the locations of the roots of one-dimensional functions
// currently uses either bisection or Newton-Rhaphson methods
// R. Sheehan 8 - 3 - 2013

typedef double (*dfunc)(double x); // pointer to a function of type double, this was the easiest way I could think of for passing a general function

// Data type for an interval [xlower, xupper]
class interval{
public:
	// Constructor
	interval(); 
	interval(double xl, double xu); 

	// Methods

	void set_xl_xu(double xl, double xu); 

	bool has_bounds(){return interval_defined; }

	double get_x_lower(){return xlower; }
	double get_x_upper(){return xupper; }

private:
	bool interval_defined; 
	double xlower; 
	double xupper; 
};

// Declare the root finder class
// Search the interval [a, b] for roots of the function func(.)
// R. Sheehan 8 - 3 - 2013

class find_root{
public:
	// Constructor
	find_root(); 
	find_root(double (*func)(double), double a, double b, int ntest); 

	// Methods
	void set_interval(double a, double b, int ntest); // change the interval properties

	//void test_func();

	void bracket_roots(); // function for bracketing the roots
	void bisection_search(double toler); // bisection method root search
	void newton_raphson_search(double toler); // Newton-Raphson method root search

private:
	// user does not need access to this function, it is called by newton_raphson_search
	double derivative(double x, double dx); // compute the value of the derivative of f at the point x by Richardson extrapolation

private:
	int nsub; // number of subintervals to search inside [a, b]
	int nroots; // the number of roots found on the interval [a, b]

	interval search_space; // this will define the interval [a, b]
	std::vector<interval> sub_intervals; // an array of sub-intervals known to contain a root

	dfunc the_func; // the function whose roots are being sought
	
	std::vector<double> the_roots; // vector to store the roots
};

#endif
#ifndef ATTACH_H
#include "Attach.h"
#endif

// Defintions of the members of interval and find_root

// interval object
// Constructors
interval::interval()
{
	// Default constructor
	interval_defined = false; 
	xlower = xupper = 0.0; 
}

interval::interval(double xl, double xu)
{
	// Constructor
	// construct an interval over the range [xl, xu]
	set_xl_xu(xl,xu); 
}

//Methods
void interval::set_xl_xu(double xl, double xu)
{
	// set the values of xlower and xupper
	// ensure that xl < xu

	try{
	
		if(fabs(xl-xu) > 1.0e-12){

			// xl != xu => interval can be created

			xlower = std::min(xl, xu);

			xupper = std::max(xl, xu); 

			interval_defined = true; 
		}
		else{
			// xl == xu throw exception

			std::string reason = "Error: void interval::set_xl_xu(double xl, double xu)\n"; 
			reason += "Attempting to construct an interval in which xl == xu\n"; 
			reason += "xl = " + template_funcs::toString(xl,4) + "\n"; 
			reason += "xu = " + template_funcs::toString(xl,4) + "\n"; 

			throw std::invalid_argument(reason); 
		}
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}	
}

// root search object
// Constructors
find_root::find_root()
{
	// Default constructor
	nsub = nroots = 0;
}

find_root::find_root(double (*func)(double), double a, double b, int ntest)
{
	// Constructor for single variable functions
	// search over ntest subintervals of [a, b] for the roots of func(x)

	try{

		if(ntest > 0){

			// require ntest > 0, otherwise throw exception

			nsub = ntest; 

			nroots = 0; 

			search_space.set_xl_xu(a, b); 

			the_func = func; 
		}
		else{
			std::string reason = "Error: find_root::find_root(double (*func)(double), double a, double b, int ntest)\n"; 
			reason += "Invalid number of search sub-intervals entered\n"; 
			reason += "ntest = " + template_funcs::toString(ntest) + "\n"; 

			throw std::invalid_argument(reason); 
		}	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

// Methods

//void find_root::test_func()
//{
//	// method to ensure that the function under consideration has been passed correctly
//	// not required
//
//	cout<<"Function value  = "<<the_func(7.0)<<endl;
//}

void find_root::set_interval(double a, double b, int ntest)
{
	// change / update the interval properties

	try{

		if(ntest > 0){

			// require ntest > 0, otherwise throw exception

			nsub = ntest; 

			nroots = 0; 

			search_space.set_xl_xu(a, b); 
		}
		else{
			std::string reason = "Error: void find_root::set_interval(double a, double b, int ntest)\n"; 
			reason += "Invalid number of search sub-intervals entered\n"; 
			reason += "ntest = " + template_funcs::toString(ntest) + "\n"; 

			throw std::invalid_argument(reason); 
		}	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}

}

void find_root::bracket_roots()
{
	// Find the sub-intervals of [a, b] known to contain a root of the equation f(x) = 0
	// Requires that a search_space [a, b] be defined

	try{

		if( search_space.has_bounds() ){
			// search_space is properly bounded, search for sub-intervals containing roots can proceed
			// nroots is not known a-priori, each of nsub sub-intervals must be tested to see if contains a root
			// in general user chooses some nsub > nroot, correct value can be found by experimentation

			double fl, fu; 
			double xl, xu; 
			double dx = ((search_space.get_x_upper() - search_space.get_x_lower()) / (nsub - 1)); 

			xl = search_space.get_x_lower(); 
			xu = xl + dx; 
			for(int i=1; i<=nsub; i++){
		
				// Evaluate the function on the sub-interval
				if(i==1){
					fl = template_funcs::Signum( the_func(xl) ); 
				}
				else{
					fl = fu; // minimise the number of function evaluations
				}

				fu = template_funcs::Signum(the_func(xu));

				// Perform the bisection test
				if(fl*fu < 0.0){
					// the sub-interval contains a root so it is stored
					nroots++; 
					sub_intervals.push_back( interval(xl, xu) ); 
				}

				// update the endpoints of the sub-interval
				xl = xu; 
				xu += dx; 
			}

			std::cout<<"The function contains "<<nroots<<" roots on the interval [ "<<search_space.get_x_lower()<<" , "<<search_space.get_x_upper()<<" ]\n"; 
			if(nroots > 0){
				std::cout<<"The roots are located in \n";
				for(int i=0; i<nroots; i++){
					std::cout<<"[ "<<sub_intervals[i].get_x_lower()<<" , "<<sub_intervals[i].get_x_upper()<<" ]\n";
				}
				std::cout<<"\n";
			}
		}
		else{
			// search_space is not bounded properly, throw exception
			std::string reason = "Error: void find_root::bracket_roots()\n"; 
			reason += "search_space is not properly bounded\n"; 
			reason += "root search cannot proceed\n"; 
			
			throw std::logic_error(reason); 
		}		
	}
	catch(std::logic_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

void find_root::bisection_search(double toler)
{
	// Locate the roots of the function f over nroots sub-intervals known to contain roots
	// to within a value of toler using the bisection method algorithm

	try{
	
		bool c1 = nroots > 0 ? true : false; // ensure that the functions contains at least one root inside search_space
		bool c2 = toler > 1.0e-16 ? true : false; // ensure that the desired tolerance is a realistic value

		if(c1 && c2){
	
			int iter = 0, max_iter = 500; 
			bool cgt = false; 
			double r=0.0, xl, xu, dx, fl, fr; 

			for(int i = 0; i<nroots; i++){
				// Loop over each sub-interval and locate the root there in

				iter = 0; 
				cgt = false; 

				xl = sub_intervals[i].get_x_lower(); 
				xu = sub_intervals[i].get_x_upper(); 

				fl = template_funcs::Signum( the_func(xl) ); 

				while(iter < max_iter){
			
					dx = 0.5*(xu-xl); 

					r = xl + dx; 

					if(fabs(dx) < toler){
						std::cout<<"Bisection algorithm has converged to within tolerance in "<<iter<<" iterations\n"; 
						the_roots.push_back(r); 
						cgt = true; 
						break; 
					}
			
					// Apply the bisection test
					// update the endpoints to the interval known to
					// contain the root
					fr = template_funcs::Signum(the_func(r)); 

					if(fl*fr > 0.0){
						xl = r; 
						fl = fr;
					}
					else{
						xu = r; 
					}

					iter++; 
				}

				if(cgt){
					std::cout<<"The root of f is located at "<<r<<"\n"; 
					std::cout<<"f ( "<<r<<" ) = "<<the_func(r)<<"\n\n";
				}
				else{
					std::cout<<"Bisection method failed to converge to a root\n\n"; 
				}

			}

		}
		else{
			// tolerance value is too low or interval does not contain any roots
			std::string reason = "Error: void find_root::bisection_search(double toler)\n"; 
			
			if(!c1){
				reason += "search_space [" + template_funcs::toString(search_space.get_x_lower(),2) + " , " 
					+ template_funcs::toString(search_space.get_x_upper(), 2) + " ] contains no roots of f\n"; 
			}

			if(!c2){
				reason += "Desired tolerance value " + template_funcs::toString(toler, 5) + " is not practical\n"; 
			}

			throw std::invalid_argument(reason); 
		}
	}
	catch(std::invalid_argument &e){
		std::cerr<<e.what()<<"\n";  
	}
}

void find_root::newton_raphson_search(double toler)
{
	// Locate the roots of the function f over sub-intervals known to contain roots
	// to within a value of toler using the Newton-Rhaphson method

	try{

		bool c1 = nroots > 0 ? true : false; // ensure that the functions contains at least one root inside search_space
		bool c2 = toler > 1.0e-16 ? true : false; // ensure that the desired tolerance is a realistic value

		if(c1 && c2){
	
			int iter = 0, max_iter = 500; 
			bool cgt = false; 
			double rnew = 0.0, rold = 0.0, xl, xu, dx, fl, fr; 

			for(int i = 0; i<nroots; i++){
				// Loop over each sub-interval and locate the root there in

				iter = 0; 
				cgt = false; 

				xl = sub_intervals[i].get_x_lower();

				xu = sub_intervals[i].get_x_upper(); 

				rold = xl; 

				while(iter < max_iter){
			
					dx = the_func(rold) / derivative(rold,0.1); 

					rnew = rold - dx; // update the root approximation

					if(fabs(dx) < toler){
						std::cout<<"Newton-Raphson algorithm has converged to within tolerance in "<<iter<<" iterations\n"; 
						the_roots.push_back(rnew); 
						cgt = true; 
						break; 
					}
			
					// Ensure that N-R alg stays inside the original interval by using a bisection step
					if(rnew > xu || rnew < xl){
					
						//cout<<"Newton-Raphson algorithm has jumped out of the search interval\n"; 
						rold = xl+0.5*(xu-xl); 

						fl = template_funcs::Signum( the_func(xl) ); 
					
						fr = template_funcs::Signum( the_func(rold) ); 
					
						if(fl*fr > 0.0){
							xl = rold;
						}
						else{
							xu = rold; 
						}

					}
					else{
						rold = rnew; 
					}

					iter++; 
				}

				if(cgt){
					std::cout<<"The root of f is located at "<<rnew<<"\n"; 
					std::cout<<"f ( "<<rnew<<" ) = "<<the_func(rnew)<<"\n\n";
				}
				else{
					std::cout<<"Newton-Raphson method failed to converge to a root\n\n"; 
				}

			}

		}
		else{
			// tolerance value is too low or interval does not contain any roots
			std::string reason = "Error: void find_root::newton_raphson_search(double toler)\n"; 
			
			if(!c1){
				reason += "search_space [" + template_funcs::toString(search_space.get_x_lower(),2) + " , " 
					+ template_funcs::toString(search_space.get_x_upper(), 2) + " ] contains no roots of f\n"; 
			}

			if(!c2){
				reason += "Desired tolerance value " + template_funcs::toString(toler, 5) + " is not practical\n"; 
			}

			throw std::invalid_argument(reason); 
		}
	
	}
	catch(std::invalid_argument &e){
		std::cerr<<e.what()<<"\n";  
	}
}

double find_root::derivative(double x, double dx)
{
	// Use Richardson extrapolation to estimate the derivative of a
	// function at the point x
	// this implementation is on the dfridr algorithm given in
	// NRinC by Press et al. 
	// this implementation uses less function evaluations than the 
	// one provided in the lectures 

	if(dx == 0.0){
		std::cout<<"Derivative cannot be computed with zero step size\n";
		return 1.0; // avoid division by zero
	}
	else{
		const int ntab = 10; 
		const double con = 1.4, con2 = template_funcs::DSQR(con); 
		const double big = std::numeric_limits<double>::max(); 
		const double safe = 2.0; 

		int i, j; 
		double err, errt, fac, hh, ans; 

		double a[ntab][ntab]; // keep the function values in a table

		hh = dx; 

		a[0][0] = (the_func(x+hh)-the_func(x-hh)) / (2.0*hh); // first approximation to f'(x)

		err = big; 

		for(i=1; i<ntab; i++){

			hh /= con; 
			
			a[0][i] = (the_func(x+hh)-the_func(x-hh)) / (2.0*hh); // approximation to f'(x) with smaller step size
			
			fac = con2; 
			
			// extrapolate the derivative to higher orders without extra function evaluations
			for(j=1; j<=i; j++){

				a[j][i] = (a[j-1][i]*fac-a[j-1][i-1]) / (fac-1.0); 
				
				fac = con2*fac; 
				
				errt = std::max( fabs(a[j][i]-a[j-1][i]), fabs(a[j][i]-a[j-1][i-1]) ); 
				
				// compute the new error with the error from the previous step
				if(errt <= err){
					err = errt; 
					ans = a[j][i]; // update the derivative estimate
				}

			}

			// if error has increased significantly stop
			if(fabs(a[i][i]-a[i-1][i-1]) >= safe*err){
				break; 
			}

		}

		//cout<<"The value of the derivative at x = "<<x<<" is "<<ans<<"\n"; 

		return ans; 
	}
}

#ifndef CURVEBASE_HPP
#define CURVEBASE_HPP

class Curvebase {
	public:
        Curvebase(double, double, bool);
		
		double x(double); //Function of s, the arc length parameter
		double y(double);
		void set_tol(double);

	protected: //Exposed to derived classes but not to other classes
		double const pmin, pmax;
		bool const rev; // orientation of the curve
		
		virtual double xp(double) = 0; // pure virtual function, has to be implemented but will not be implemented in base class
		virtual double yp(double) = 0; // check if a<p<b
		virtual double dxp(double) = 0;
		virtual double dyp(double) = 0;
		
		double integrate(double); // arc length integral, has to be defined in base class, can be overwritten i guess?
	
	private:
		double ASI_routine(double, double, double, double, double, double, double, double);
		double integrand(double); //Function to go in the integral of the standard definition of integrate
		double simp(double,double,double,double,double);
		double newton(double, double);
		double tollerance;
};

#endif
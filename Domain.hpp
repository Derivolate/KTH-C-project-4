#ifndef DOMAIN_HPP
#define DOMAIN_HPP

#include "Curvebase.hpp"
#include <string>

class Domain{
	public:
		Domain(Curvebase&, Curvebase&, Curvebase&, Curvebase&, double =  1e-12);
		Domain(const Domain&);
		~Domain();

		Domain& operator =(Domain&); //Assignment operator
	
		void generate_grid(int, int, int); // pass m*n grid and interpolation option. 1-linear, 2-stretched vertical indices
		void print_grid(bool,bool,bool,std::string ="");
	private:
		Curvebase * sides[4];
        double * x_ = nullptr, * y_ = nullptr;
        int n_ = 0, m_ = 0;
		double tol;
		double phi1(double); //transition function for grid generation phix(0) = 1, phix(1) = 0
		double phi2(double); //transition function for grid generation phiy(1) = 1, phix(0) = 0
};

#endif
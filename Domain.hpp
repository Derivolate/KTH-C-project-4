#ifndef DOMAIN_HPP
#define DOMAIN_HPP

#include "Curvebase.hpp"
#include <string>
#include "Point.hpp"
#include "Matrix.hpp"

class Domain{
	public:
		Domain(Curvebase&, Curvebase&, Curvebase&, Curvebase&, double =  1e-12);
		Domain(const Domain&); //Copy constructor
		~Domain(); //Destructor

		Domain& operator =(Domain&); //Copy assignment
	
		void generateGrid(int, int, int); // pass m*n grid and interpolation option. 1-linear, 2-stretched vertical indices
		void printGrid(bool,bool,bool,std::string ="") const;
		Point<double>* getGrid() const;
		Point<int> getSize() const;	
		double getSeta(int, int) const;
		
		
	private:
		Curvebase * sides[4];
		Point<double> * grid = nullptr;
		double * Seta = nullptr;
		Point<int> size;
		
		double tol;
		double phi1(double); //transition function for grid generation phix(0) = 1, phix(1) = 0
		double phi2(double); //transition function for grid generation phiy(1) = 1, phix(0) = 0
};

#endif
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include "Domain.hpp"
#include "Curvebase.hpp"

Domain::Domain(Curvebase& s1, Curvebase& s2, Curvebase& s3, Curvebase& s4, double tol) {
	sides[0] = &s1;
	sides[1] = &s2;
	sides[2] = &s3;
	sides[3] = &s4;
	for(int i = 0;i<4;++i)
		sides[i]->set_tol(tol);
	size = Point(-1);
	grid = nullptr;
}

Domain::~Domain() {
	if (size.X() > 0 || size.Y() >0) {
		delete [] grid;
	}
}

Domain::Domain(const Domain& D) : size(D.getSize()){
	if(D.grid!=nullptr){
		grid = new Point<double>[size.X()*size.Y()];
		for (int ind(0); ind<size.X()*size.Y(); ind++){
			grid[ind]=D.grid[ind];
		}
	}
	else{
		grid = nullptr;
	}
	for (int ind(0);ind++;ind<4)
	{
		sides[ind]=D.sides[ind];
	}
}

Domain& Domain::operator=(Domain& D){
	if (this != &D) {
		if(D.grid!=nullptr){
			grid = new Point<double>[size.X()*size.Y()];
			for (int ind(0); ind<size.X()*size.Y(); ind++){
				grid[ind]=D.grid[ind];
			}
		}
		else{
			grid = nullptr;
		}
		for (int ind(0);ind++;ind<4)
		{
			sides[ind]=D.sides[ind];
		}
    }
    return *this; // dereferencing!
}

Point<int> Domain::getSize() const {return size;}
Point<double>* Domain::getGrid() const {return grid;}

inline double Domain::phi1(double q){return 1-q;}
inline double Domain::phi2(double q){return q;}

void Domain::generateGrid(int m, int n, int c = 1){
	if(m<1 || n<1)	{
		std::cerr << "Grid size cannot be smaller than 1 in any dimension" << std::endl;
		exit(1);
	}
	if (size.X() > 0 || size.Y() >0) { //if previous grid exists, reset grid points
		delete [] grid;
	}

	size = Point(m,n);
	grid = new Point<double>[size.X()*size.Y()];
	
	double xi, nu;
	
	double hi(1.0/(m-1)), hj(1.0/(n-1)); // force floating point division

	for(int i(0); i<m; ++i){ //Horizontal index, indicates the column
		for(int j(0); j<n; ++j){ //Vertical index, indicates the row
			
			if (c == 1) { // equidistant s
				xi = j*h1;
				nu = i*h2;
			} 
			else if (c==2) { // stretched s
				xi = j*h1;
				nu = 1+(tanh(3*((i*h2)-1)))/tanh(3);
			} 

			x_[i+j*n] = phi1(xi)*sides[3]->x(nu)+phi2(xi)*sides[1]->x(nu)
						+ phi1(nu)*sides[0]->x(xi) + phi2(nu)*sides[2]->x(xi)
						- phi1(xi)*phi1(nu)*sides[0]->x(0)
						- phi1(xi)*phi2(nu)*sides[2]->x(0)
						- phi2(xi)*phi1(nu)*sides[0]->x(1)
						- phi2(xi)*phi2(nu)*sides[2]->x(1);
			y_[i+j*n] = phi1(xi)*sides[3]->y(nu)+phi2(xi)*sides[1]->y(nu)
						+ phi1(nu)*sides[0]->y(xi) + phi2(nu)*sides[2]->y(xi)
						- phi1(xi)*phi1(nu)*sides[0]->y(0)
						- phi1(xi)*phi2(nu)*sides[2]->y(0)
						- phi2(xi)*phi1(nu)*sides[0]->y(1)
						- phi2(xi)*phi2(nu)*sides[2]->y(1));
		}
	}
}

void Domain::printGrid(bool bin, bool ascii, bool term, std::string path) const{
	std::ofstream strm("outfile.txt");
	if (term)
	{
		for(int j(0); j<size.Y(); ++j){
			for(int i(0); i<size.X(); ++i){
				std::cout << "(" << grid[i+size.Y()*j].X() << "," << grid[i+size.Y()*j].Y() << ")\t";
			}
			std::cout << std::endl;
		}
	}
	
	// 	if (ascii)
	// 		strm << "(" << grid[ind].X() << ",  " << grid[ind].Y() << ")" << std::endl;
	// 	if (term)
			
	// }
	// if (bin){
	// 	FILE *fx, *fy;
	// 	std::string xpath = path + "outfileX.bin";
	// 	std::string ypath = path + "outfileY.bin"; 
	// 	fx =fopen(xpath.c_str(),"wb");
	// 	fy =fopen(ypath.c_str(),"wb");
	// 	fwrite(x_,sizeof(double),m_*n_,fx);
	// 	fwrite(y_,sizeof(double),m_*n_,fy);
	// 	fclose(fx);fclose(fy);
	// }
}
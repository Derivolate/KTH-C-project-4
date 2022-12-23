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
	m_ = n_ = 0;
	x_ = y_ = nullptr;
}

Domain::~Domain() {
	if (m_ > 0 || n_ >0) {
		delete [] x_;
		delete [] y_;
	}
}

int Domain::get_m(){return m_;}
int Domain::get_n(){return n_;}

double* Domain::get_x(){return x_;}
double* Domain::get_y(){return y_;}

inline double Domain::phi1(double q){return 1-q;}
inline double Domain::phi2(double q){return q;}


void Domain::generateGrid(int m, int n, int c = 1){
	if(m<1 || n<1)	{
		std::cerr << "Grid size cannot be smaller than 1 in any dimension" << std::endl;
		exit(1);
	}
	if(m_>0 || n_>0) { //if previous grid exists, reset grid points
		delete [] x_;
		delete [] y_;
	}

	x_ = new double[m*n]; 
	y_ = new double[m*n];
	double xi, nu;
	m_ = m, n_= n;
	
	double h1(1.0/(n-1)), h2(1.0/(m-1)); // force floating point division

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
						- phi2(xi)*phi2(nu)*sides[2]->y(1);
		}
	}
}

void Domain::printGrid(bool bin, bool ascii, bool term, std::string path){
	std::ofstream strm("outfile.txt");
	for(int ind = 0; ind<n_*m_; ++ind){
		if (ascii)
			strm << "(" << x_[ind] << ",  " << y_[ind] << ")" << std::endl;
		if (term)
			std::cout << "(" << x_[ind] << ",  " << y_[ind] << ")" << std::endl;
	}
	if (bin){
		FILE *fx, *fy;
		std::string xpath = path + "outfileX.bin";
		std::string ypath = path + "outfileY.bin"; 
		fx =fopen(xpath.c_str(),"wb");
		fy =fopen(ypath.c_str(),"wb");
		fwrite(x_,sizeof(double),m_*n_,fx);
		fwrite(y_,sizeof(double),m_*n_,fy);
		fclose(fx);fclose(fy);
	}
}
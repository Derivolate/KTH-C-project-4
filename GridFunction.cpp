#include "Matrix.hpp"
#include "Domain.hpp"
#include "GridFunction.hpp"
GridFunction::GridFunction(Domain* _grid) : grid(_grid){}

GridFunction::~GridFunction(){
    //Delete u
    //Delete grid
}

void GridFunction::printFkt(){
    u.printMatrix();
}

// This functionality cannot be in the constructor because it calls a virtual function
void GridFunction::fillGrid(){
    int n(grid->get_n()), m(grid->get_m());
    double* x(grid->get_x());
    double* y(grid->get_y());
    double val;
    u = Matrix(n);
    for(int i(0);i<n;++i){
        for(int j(0);j<m;++j){
            val = u_function(x[i+m*j],y[i+m*j]); //TODO double check if these indices are correct
            u.setElem(val,i,j);
        }
    }
}
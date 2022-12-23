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

void GridFunction::Dx() {
    double* x(grid->get_x());
    double* y(grid->get_y());
    int m(u.getSize()), n(u.getSize());
    double h1(1.0/(n-1)), h2(1.0/(m-1));
    Matrix dxxi(n), dxeta(n), dyxi(n), dyeta(n), duxi(n), dueta(n);
    for (int i(0); i<m; i++){
            dxxi.setElem((x[i+2*m] - 4*x[i+m] + 3*x[i])/(3*h1),i,0);
            dyxi.setElem((y[i+2*m] - 4*y[i+m] + 3*y[i])/(3*h1),i,0);
            dxxi.setElem((x[i+2*m] - 4*x[i+m] + 3*x[i])/(3*h1),i,n);
            dyxi.setElem((y[i+2*m] - 4*y[i+m] + 3*y[i])/(3*h1),i,n);

    for (int i(0); i<m; i++){
        for (int j(0); j<n; j++){
            if (j == 0){ // Deal with boundary0 
                dxxi.setElem((x[i+2*m] - 4*x[i+m] + 3*x[i])/(3*h1),i,j);
                dyxi.setElem((y[i+2*m] - 4*y[i+m] + 3*y[i])/(3*h1),i,j);
                if (i == 0){
                    dxeta.setElem((x[i+2] - 4*x[i+1] + 3*x[i])/(3*h2),i,j);
                    dyeta.setElem((y[i+2] - 4*y[i+1] + 3*y[i])/(3*h2),i,j);
                }
                else if (i == m-1){
                    dxeta.setElem((x[i-2] - 4*x[i-1] + 3*x[i])/(3*h2),i,j);
                    dyeta.setElem((y[i-2] - 4*y[i-1] + 3*y[i])/(3*h2),i,j);
                }
                else{
                    dxeta.setElem((x[i-2] - 4*x[i-1] + 3*x[i])/(3*h2),i,j);
                    dyeta.setElem((y[i-2] - 4*y[i-1] + 3*y[i])/(3*h2),i,j);
                }
                if
                u.setElem(dx, i,j);
            }    
        }   
    }
}
#include "Matrix.hpp"
#include "Domain.hpp"
#include "GridFunction.hpp"
#include <cmath>

GridFunction::GridFunction(std::shared_ptr<Domain> _dom) : domain(_dom),u(_dom->getSize()){}

void GridFunction::printFkt(){
    u.printMatrix();
}

// This functionality cannot be in the constructor because it calls a virtual function
void GridFunction::fillGrid(){
    Point<int> size =  domain->getSize();
    Point<double>* coords = domain->getGrid();
    Point<double> coord;
    double val;
    for(int i(0);i<u.getSize().Y();++i){
        for(int j(0);j<u.getSize().X();++j){
            coord = coords[i+u.getSize().Y()*j];
            val = u_function(coord);
            u.setElem(val,i,j);
        }
    }
}

double GridFunction::u_function(Point<double> P)
{
    return std::sin(std::pow(P.X()/10,2))*std::cos(P.X()/10)+P.Y();
}
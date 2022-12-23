#include <iostream>
#include <string>

#include "Hline.hpp"
#include "Vline.hpp"
#include "BumpedCurve.hpp"
#include "Domain.hpp"
#include "GridFunction.hpp"
#include "Matrix.hpp"
#include <memory>

using namespace std;

int main(int argc, char **argv)
{
    Point size = Point(5,5);
    Matrix M = Matrix(size);
    Matrix N(M);
    
    M.setElem(2,0,0);
    N.setElem(3,2,0);

    // M.printMatrix();
    Matrix P(N);
    M.printMatrix();
    P.printMatrix();
    P+=M;
    P.printMatrix();
    
    // P*=3;
    // P.printMatrix();
    // M%=P;
    // M.printMatrix();
    // // // M.printMatrix();

    double ox(0),oy(0),width(2), height(2);
    
    // BumpedCurve border1 = BumpedCurve(0,1,true,ox,oy,width,-6,0,3,3,0.5);
    Hline border1 = Hline(0,1,true,ox,oy,width);
    Vline border2 = Vline(0,1,true,ox+width,oy,height);
    Hline border3 = Hline(0,1,true,ox,oy+height,height);
    Vline border4 = Vline(0,1,true,ox,oy,height);

    Domain domain = Domain(border1,border2,border3,border4);
    domain.generateGrid(20,15,1);
    Domain test1(domain);
    Domain test2 = domain;


    std::shared_ptr<Domain> domain_ptr= std::make_shared<Domain>(domain);

    GridFunction Fkt = GridFunction(domain_ptr);
    Fkt.fillGrid();
    domain.printGrid(true,true,true);
    Fkt.printFkt();

}
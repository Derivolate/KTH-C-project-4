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

int main(int argc, char **argv){
    double ox(-10),oy(0),width(15), height(3);
    
    BumpedCurve border1 = BumpedCurve(0,1,true,ox,oy,width,-6,0,3,3,0.5);
    Vline border2 = Vline(0,1,true,ox+width,oy,height);
    Hline border3 = Hline(0,1,true,ox,oy+height,height);
    Vline border4 = Vline(0,1,true,ox,oy,height);

    Domain domain = Domain(border1,border2,border3,border4);
    domain.generateGrid(20,20,2);
    
    std::shared_ptr<Domain> domain_ptr= std::make_shared<Domain>(domain);
    std::cout << "SETA " <<  domain.getSeta(0,2) << std::endl;
   


    GridFunction Fkt = GridFunction(domain_ptr);
    

    std::cout << "============ Filling Grid ============" << std::endl;
    Fkt.fillGrid();
    GridFunction dFdy = Fkt.Dy();
    GridFunction dFdx = Fkt.Dx();
    GridFunction LF = dFdy.Dy() + dFdx.Dx();
    
    std::cout << "============ Printing Grid ============" << std::endl;
    domain.printGrid(true,false,true);
    std::cout << "============ Printing Function ============" << std::endl;
    std::cout << "----------- Grid Function -----------" << std::endl;
    Fkt.printFkt("u_out.bin");
    std::cout << "----------- DX -----------" << std::endl;
    dFdx.printFkt("dx_out.bin");
    std::cout << "----------- DY -----------" << std::endl;
    dFdy.printFkt("dy_out.bin");
    std::cout << "----------- DDxy -----------" << std::endl;
    LF.printFkt("ddxy_out.bin");
}
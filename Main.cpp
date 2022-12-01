#include <iostream>
#include <string>

#include "Hline.hpp"
#include "Vline.hpp"
#include "BumpedCurve.hpp"
#include "Domain.hpp"

using namespace std;

int main(int argc, char **argv)
{
    double ox(-10),oy(0),width(15), height(3);
    
    BumpedCurve border1 = BumpedCurve(0,1,true,ox,oy,width,-6,0,3,3,0.5);
    Vline border2 = Vline(0,1,true,ox+width,oy,height);
    Hline border3 = Hline(0,1,true,ox,oy+height,height);
    Vline border4 = Vline(0,1,true,ox,oy,height);

    Domain domain = Domain(border1,border2,border3,border4);
    
    domain.generate_grid(20,50,1);

    domain.print_grid(true,true,true);
}
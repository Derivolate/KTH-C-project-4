#include <cmath>
#include "ProjFunction.hpp"
#include "Domain.hpp"
// using namespace std;

ProjFunction::ProjFunction(Domain* grid) : GridFunction(grid){
    fillGrid();
}

double ProjFunction::u_function(double x, double y)
{
    return std::sin(std::pow(x/10,2))*std::cos(x/10)+y;
}
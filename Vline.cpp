#include "Vline.hpp"

//A horizontal curve running from p = a to p = b, to the right (dir = true) or to left (dir = false) starting at Point begin running for length len
Vline::Vline(double a, double b, bool dir, double x, double y, double len) 
    : Curvebase(a, b, dir), bx(x), by(y), length(len), ex(x), ey(y+len) {
}

double Vline::xp(double p){
    return bx;
}
double Vline::yp(double p){
    return by+(p-pmin)*(ey-by)/(pmax-pmin);
}

double Vline::dxp(double p){
    return 0;
}

double Vline::dyp(double p){
    return (ey-by)/(pmax-pmin);
}
#ifndef VLINE_HPP
#define VLINE_HPP

#include "Curvebase.hpp"

class Vline : public Curvebase{
    public :
        Vline(double, double, bool, double, double, double);
    protected :
        double xp(double);
        double yp(double);
        double dxp(double);
        double dyp(double);
    private :
        double bx,by,length; //The x-coordinate and y-coordinate of the beginning of the line, and the length of the line
        double ex,ey; //The x-coordinate and y-coordinate of the end of the line
};

#endif
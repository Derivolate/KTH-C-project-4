#ifndef HLINE_HPP
#define HLINE_HPP

#include "Curvebase.hpp"

class Hline : public Curvebase {
    public :
        Hline(double, double, bool, double, double, double);
    
    protected :
        double xp(double) override;
        double yp(double) override;
        double dxp(double) override;
        double dyp(double) override;

    private :
        double bx,by,length; //The x-coordinate and y-coordinate of the beginning of the line, and the length of the line
        double ex,ey; //The x-coordinate and y-coordinate of the end of the line
};

#endif
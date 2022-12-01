#ifndef BCURVE_HPP
#define BCURVE_HPP

#include "Curvebase.hpp"

class BumpedCurve : public Curvebase {
    public :
        BumpedCurve(double,double,bool,double,double,double,double,double,double,double,double);
    
    protected :
        double xp(double) override;
        double yp(double) override;
        double dxp(double) override;
        double dyp(double) override;

    private :
        double bx,by,length; //The x-coordinate and y-coordinate of the beginning of the line, and the length of the line
        double ex,ey; //The x-coordinate and y-coordinate of the end of the line
        double s1, s2, x1, x2, h, xc;
};
#endif
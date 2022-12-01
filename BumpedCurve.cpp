#include "BumpedCurve.hpp"
#include <iostream>
#include <cmath>
//a, b: minimal and maximal value for p
//rev:  if the curve is reversed or not, currently unused
//x0, y0: coordinates of the start point of the curve
//len: length of the curve if the bump has a height of 0
//x1_,x2_: x positions of the beginning and end of the bump
//s1_,s2_: slope of the beginning and end of the bump
//h_: height of the bump
//xc: the x-coordinate of the crossing of the two tangent curves
BumpedCurve::BumpedCurve(double a, double b, bool rev, double x0, double y0, double len, double x1_, double x2_, double s1_, double s2_, double h_)
    :Curvebase(a, b, rev),bx(x0),by(y0),ex(x0+len),ey(y0),x1(x1_),x2(x2_),s1(s1_),s2(s2_),h(h_),xc(((s1_/s2_)*x1_+x2_)/(1+(s1_/s2_))){}

double  BumpedCurve::xp(double p)
{
    return bx+(p-pmin)*(ex-bx)/(pmax-pmin);
}

double  BumpedCurve::yp(double p)
{
    double yval;
    double xval(xp(p));
    if (xval<xc)
    {
        yval = -s1*(xval-x1);
    }
    else 
    {
        yval =  s2*(xval-x2);
    }
    return h/(1+exp(yval));
}

double  BumpedCurve::dxp(double p)
{
    return (ex-bx)/(pmax-pmin);
}

double  BumpedCurve::dyp(double p)
{
    double sval, expval;
    double xval(xp(p));
    if (xval<xc)
    {
        expval = exp(-s1*(xval-x1));
        sval = -s1;
    }
    else
    {
        expval = exp(s2*(xval-x2));
        sval = s2;
    }
    return -(h*expval*sval)/pow(1+expval,2);
}
#include "Curvebase.hpp"
#include <cmath>
#include <iostream>

Curvebase::Curvebase(double a, double b, bool dir) : pmin(a), pmax(b), rev(dir) {}

Curvebase::~Curvebase(){}

void Curvebase::set_tol(double tol)
{
    tollerance = tol;
}

double Curvebase::integrate(double p) {
    // Calculate values to kickstart the algorithm
    double a = pmin, c = p;
    double b = a+c/2;
    double fa = integrand(a), fb = integrand(b), fc = integrand(c);
    double I = simp(fa,fb,fc,a,c);
    // Run the algorithm
    return ASI_routine(a,b,c,tollerance,fa,fb,fc,I);
}



double Curvebase::ASI_routine(double a, double b, double  c, double tol, double fa, double fb, double fc, double I1){
    //Calculate the midpoints
    double ab((a+b)/2), bc((b+c)/2);
    //Evaluate the midpoints
    double fab(integrand(ab)), fbc(integrand(bc));
    //Compute the midpoint integrals
    double Iab(simp(fa,fab,fb,a,b)), Ibc(simp(fb, fbc, fc, b,c));
    double I2 = Iab + Ibc;
    
    //Compare the new result, and return it if it is accurate enough
    double err(fabs(I1-I2));
    if( err < 15*tol)
    { 
        return I2;
    }
    else
    {
        //Recursive function call if the result is not accurate enough
        return ASI_routine(a,ab,b,tol/2,fa,fab,fb,Iab) + ASI_routine(b,bc,c,tol/2,fb,fbc,fc,Ibc);
    }
}

//Inline function of the integrand
inline double Curvebase::integrand(double q){
    return sqrt(pow(dxp(q),2)+pow(dyp(q),2));
}
//Inline function for simpsons rule using pre-calculated function evaluations
inline double Curvebase::simp(double fa,double fb,double fc,double a,double c)
{
    return (fa+4*fb+fc)*(c-a)/6;
}

double Curvebase::newton(double s, double guess){
    double p(guess); //initial guess used to start Newtons method. 
    double dp(1.0);  //delta between the steps used to evaluate tolerance (can be improved by averaging over last x steps)
    int i(0); //iteration count to avoid non converging sequences
    int timeout(1e5); //timeout at which the newton method stops assuming no convergance
    double arclenth(integrate(pmax));
    while(fabs(dp) > tollerance){ //begin newton iterations
        dp = (integrate(p) - s*arclenth)/sqrt(dxp(p)*dxp(p)+dyp(p)*dyp(p));
        p -= dp;
        if (p<pmin){p = pmin;}
        else if(p>pmax){p = pmax;}
        ++i;
    }
    if (!(pmin<=p and p<=pmax)){
        std::cerr<<"Value of p outside of allowed range"<<std::endl;
        exit(1);
    }
    return p;
}

double Curvebase::x(double s){
    double guess((pmin+pmax)/2);
    double p0;
    if (!rev) p0 = newton(1.0-s, guess);
    else p0 = newton(s,guess);
    return xp(p0);
}

double Curvebase::y(double s){
    double guess((pmin+pmax)/2);
    double p0;
    if (!rev) p0 = newton(1.0-s, guess);
    else p0 = newton(s,guess);
    return yp(p0);
}
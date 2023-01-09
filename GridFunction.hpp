#ifndef GRIDFUNCTION_HPP
#define GRIDFUNCTION_HPP
#include <cmath>
#include "Matrix.hpp"
#include "Domain.hpp"
#include <memory>

class GridFunction {
    public:
        //Smart pointer can take care on its own of proper copying and deletion. Rule of 0 applies
        // GridFunction(const GridFunction&) = default;
        // ~GridFunction() = default;
        // GridFunction& operator=(const GridFunction&) = default;   //Assignment operator

        GridFunction(std::shared_ptr<Domain>);
        GridFunction(GridFunction const &);
        GridFunction& operator+=(const GridFunction&);  //Pointwise grid addition assignment operator
        // GridFunction& operator+(const GridFunction&) const;  //Pointwise grid addition operator
        GridFunction& operator*=(const GridFunction&);  //Pointwise grid multiplication assignment operator
        // GridFunction& operator*(const GridFunction&) const;  //Pointwise grid multiplication operator
        GridFunction operator+(const GridFunction&);
        
        GridFunction Dx();    //X derivatives
        GridFunction Dy();    //Y derivative
        GridFunction DDxy(GridFunction * Dx, GridFunction * Dy);  //Laplace function
        
        void printFkt(std::string) const;
        void fillGrid();

    private:
        std::shared_ptr<Domain> domain;
        Matrix u;
        Matrix phix, phiy;

        double u_function(Point<double>);

        // Derivative Tools
        template<int dx, int dy> double center_(int i, int j, Matrix const & fctn);
        template<int dx, int dy> double asym_(int i, int j, Matrix const & fctn);
        template<int dx, int dy> double face_(int i, int j);
        template<int dx, int dy, int dxi, int deta> double vertex_(int i, int j);
        template<int dx, int dy, int dxi, int deta> double edge_(int i, int j);


        // TODO template<int dx, int dy> double ghost_(int i, int j);

}; 

// ----------------------------------

// Definition of Center Scheme for 1st Derivative
template<int dx, int dy> double GridFunction::center_(int i, int j, Matrix const & fctn) {
    if (!(dx || dy)){return 0;}
    double const d = dy*std::fabs(domain->getSeta(i-dy,j)-domain->getSeta(i+dy,j)) + 2*dx/(domain->getSize().Y()-1.0);
    return(fctn.getElem(i+dy,j+dx) - fctn.getElem(i-dy,j-dx))/d;
}

// Definition of Asymmetric Scheme for 1st Derivative
template<int dx, int dy> double GridFunction::asym_(int i, int j, Matrix const & fctn) {
    if (!(dx || dy)){return 0;}
    double const d = 0.5*dy*std::fabs(domain->getSeta(i+2*dy,j)-domain->getSeta(i,j)) + dx/(domain->getSize().Y()-1.0);
    return (fctn.getElem(i+2*dy,j+2*dx) - 4*fctn.getElem(i+dy,j+dx) + 3*fctn.getElem(i,j))/(3*d);
}

// Vertex Derivative Calculation: 
// dx,dy either or 1,0 to determine derivative direction
// dxi,deta defines which directions have to be calculated as asym: 1/-1 pos/neg dir asym
template<int dx, int dy, int dxi, int deta> double GridFunction::vertex_(int i, int j){
    double const dudeta = asym_<0,deta>(i,j, u);
    double const dudxi = asym_<dxi,0>(i,j, u);
    double const dphixdxi = asym_<dxi,0>(i,j, phix);
    double const dphiydxi = asym_<dxi,0>(i,j, phiy);
    double const dphixdeta = asym_<0,deta>(i,j, phix);
    double const dphiydeta = asym_<0,deta>(i,j, phiy);
    double const detJ = dphixdxi * dphiydeta - dphixdeta*dphiydxi;
    return (1/detJ) * (dx*(dudxi*dphiydeta - dudeta*dphiydxi)+dy*(dudeta*dphixdxi - dudxi*dphixdeta));
}

// Edge Derivative Calculation: 
// dx,dy either or 1,0 to determine derivative direction
// dxi,deta defines which directions have to be calculated as asym: 1/-1 pos/neg dir asym, 0 center
template<int dx, int dy, int dxi, int deta> double GridFunction::edge_(int i, int j){
    double const dudeta = asym_<0,deta>(i,j,u) * abs(deta) + center_<0,1>(i,j,u) * (1-abs(deta));
    double const dudxi = asym_<dxi,0>(i,j, u)*abs(dxi) + center_<1,0>(i,j,u)* (1-abs(dxi));
    double const dphixdxi = asym_<dxi,0>(i,j, phix) * abs(dxi) + center_<1,0>(i,j, phix)* (1-abs(dxi));
    double const dphiydxi = asym_<dxi,0>(i,j, phiy) * abs(dxi) + center_<1,0>(i,j, phiy)* (1-abs(dxi));
    double const dphixdeta = asym_<0,deta>(i,j, phix)* abs(deta) + center_<0,1>(i,j, phix)* (1-abs(deta));
    double const dphiydeta = asym_<0,deta>(i,j, phiy)* abs(deta) + center_<0,1>(i,j, phiy)* (1-abs(deta));
    double const detJ = dphixdxi * dphiydeta - dphixdeta*dphiydxi;
    return (1/detJ) * (dx*(dudxi*dphiydeta - dudeta*dphiydxi)+dy*(dudeta*dphixdxi - dudxi*dphixdeta));
}

// Face Derivative Calculation: 
// dx,dy either or 1,0 to determine derivative direction with center
template<int dx, int dy> double GridFunction::face_(int i, int j){
    double const dudeta = center_<0,1>(i,j, u);
    double const dudxi = center_<1,0>(i,j, u);
    double const dphixdxi = center_<1,0>(i,j, phix);
    double const dphiydxi = center_<1,0>(i,j, phiy);
    double const dphixdeta = center_<0,1>(i,j, phix);
    double const dphiydeta = center_<0,1>(i,j, phiy);
    double const detJ = dphixdxi * dphiydeta - dphixdeta*dphiydxi;
    return (1/detJ) * (dx*(dudxi*dphiydeta - dudeta*dphiydxi)+dy*(dudeta*dphixdxi - dudxi*dphixdeta));
}

#endif
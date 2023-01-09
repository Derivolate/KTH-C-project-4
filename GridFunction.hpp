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

        template<int dx, int dy> double center_(int i, int j, Matrix const & fctn);
        template<int dx, int dy> double assym_(int i, int j, Matrix const & fctn);
        template<int dx, int dy> double face_(int i, int j);
        template<int dx, int dy, int dxi, int deta> double vertex_(int i, int j);
        template<int dx, int dy, int dxi, int deta> double edge_(int i, int j);


        // TODO template<int dx, int dy> double ghost_(int i, int j);

}; 

// ----------------------------------

template<int dx, int dy> double GridFunction::center_(int i, int j, Matrix const & fctn) {
    std::cout << "TEMP SETA" << std::endl;
    if (!(dx || dy)){return 0;}
    double const d = dy*std::fabs(domain->getSeta(i-dy,j)-domain->getSeta(i+dy,j)) + 2*dx/(domain->getSize().Y()-1.0);
    // double const d = dy/(domain->getSize().X()-1.0) + dx/(domain->getSize().Y()-1.0);
    return(fctn.getElem(i+dy,j+dx) - fctn.getElem(i-dy,j-dx))/d;
}

template<int dx, int dy> double GridFunction::assym_(int i, int j, Matrix const & fctn) {
    if (!(dx || dy)){return 0;}
    double const d = 0.5*dy*std::fabs(domain->getSeta(i+2*dy,j)-domain->getSeta(i,j)) + dx/(domain->getSize().Y()-1.0);
    // double const d = dy/(domain->getSize().X()-1.0) + dx/(domain->getSize().Y()-1.0);
    return (fctn.getElem(i+2*dy,j+2*dx) - 4*fctn.getElem(i+dy,j+dx) + 3*fctn.getElem(i,j))/(3*d);
}

template<int dx, int dy, int dxi, int deta> double GridFunction::vertex_(int i, int j){
    double const dudeta = assym_<0,deta>(i,j, u);
    double const dudxi = assym_<dxi,0>(i,j, u);
    double const dphixdxi = assym_<dxi,0>(i,j, phix);
    double const dphiydxi = assym_<dxi,0>(i,j, phiy);
    double const dphixdeta = assym_<0,deta>(i,j, phix);
    double const dphiydeta = assym_<0,deta>(i,j, phiy);
    double const detJ = dphixdxi * dphiydeta - dphixdeta*dphiydxi;
    // std::cout << "Vertex Element: [" << i << "," << j << "] direction: dxi,deta [" << dxi << "," << deta << "]" << std::endl;
    // std::cout << "dphixdxi: " << dphixdxi << " dphiydeta: " << dphiydeta << " dphixdeta: " << dphixdeta << " dphiydxi: "<< dphiydxi  << " detJ: " << detJ << std::endl;
    return (1/detJ) * (dx*(dudxi*dphiydeta - dudeta*dphiydxi)+dy*(dudeta*dphixdxi - dudxi*dphixdeta));
}

template<int dx, int dy, int dxi, int deta> double GridFunction::edge_(int i, int j){
    double const dudeta = assym_<0,deta>(i,j,u) * abs(deta) + center_<0,1>(i,j,u) * (1-abs(deta));
    double const dudxi = assym_<dxi,0>(i,j, u)*abs(dxi) + center_<1,0>(i,j,u)* (1-abs(dxi));
    double const dphixdxi = assym_<dxi,0>(i,j, phix) * abs(dxi) + center_<1,0>(i,j, phix)* (1-abs(dxi));
    double const dphiydxi = assym_<dxi,0>(i,j, phiy) * abs(dxi) + center_<1,0>(i,j, phiy)* (1-abs(dxi));
    double const dphixdeta = assym_<0,deta>(i,j, phix)* abs(deta) + center_<0,1>(i,j, phix)* (1-abs(deta));
    double const dphiydeta = assym_<0,deta>(i,j, phiy)* abs(deta) + center_<0,1>(i,j, phiy)* (1-abs(deta));
    double const detJ = dphixdxi * dphiydeta - dphixdeta*dphiydxi;
    return (1/detJ) * (dx*(dudxi*dphiydeta - dudeta*dphiydxi)+dy*(dudeta*dphixdxi - dudxi*dphixdeta));
}

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

// template<int dx, int dy> double GridFunction::center_(int i, int j){
//     if (dx*dy != 0){
//         std::cout << "Only x or y may be given as direction. No inbetween" << std::endl;
//         return(0);
//     }
//     const double deta = 1.0/(domain->getSize().X()-1);
//     const double dxi = 1.0/(domain->getSize().Y()-1);
//     const double dudeta = (u.getElem(i+1,j) - u.getElem(i-1,j))*0.5/deta;
//     const double dudxi = (u.getElem(i,j+1) - u.getElem(i,j-1))*0.5/dxi;
//     const double dphixdxi = (phix.getElem(i,j+1) - phix.getElem(i,j-1))*0.5/dxi;
//     const double dphiydxi = (phiy.getElem(i,j+1) - phiy.getElem(i,j-1))*0.5/dxi;
//     const double dphixdeta = (phix.getElem(i+1,j) - phix.getElem(i-1,j))*0.5/deta;
//     const double dphiydeta = (phiy.getElem(i+1,j) - phiy.getElem(i-1,j))*0.5/deta;
//     const double detJ = dphixdxi * dphiydeta - dphixdeta*dphiydxi;
//     std::cout << "Element: [" << i << "," << j << "] direction: center[" << dx << "," << dy << "] phiy.getElem: " << phiy.getElem(i+1,j) << std::endl;
//     std::cout<< "dxi: " << dxi <<"dphixdxi: " << dphixdxi << " dphiydeta: " << dphiydeta << " dphixdeta: " << dphixdeta << " dphiydxi: "<< dphiydxi  << " detJ: " << detJ << std::endl;
//     if (detJ == 0){
//         std::cout << "Matrix is singular, cannot calculate derivative" << std::endl;
//         return(0);
//     }
//     return((1/detJ) * (dx*(dudxi*dphiydeta - dudeta*dphiydxi)+dy*(dudeta*dphiydxi - dudxi*dphiydeta)));
// }

// template<int dx, int dy> double GridFunction::assym_(int i, int j){
//     if (dx*dy != 0){
//         std::cout << "Only x or y may be given as direction. No inbetween" << std::endl;
//         return(0);
//     }
//     const double deta = 1.0/(domain->getSize().X()-1);
//     const double dxi = 1.0/(domain->getSize().Y()-1);
//     const double dudeta = ((u.getElem(i+2*dy,j) - 4*u.getElem(i+dy,j)) + 3*u.getElem(i,j))/(3*deta);
//     const double dudxi = ((u.getElem(i,j+2*dx) - 4*u.getElem(i,j+dx)) + 3*u.getElem(i,j))/(3*dxi);
//     const double dphixdxi = ((phix.getElem(i,j+2*dx) - 4*phix.getElem(i,j+dx)) + 3*phix.getElem(i,j))/(3*dxi);
//     const double dphiydxi = ((phiy.getElem(i,j+2*dx) - 4*phiy.getElem(i,j+dx)) + 3*phiy.getElem(i,j))/(3*dxi);
//     const double dphixdeta = ((phix.getElem(i+2*dy,j) - 4*phix.getElem(i+dy,j)) + 3*phix.getElem(i,j))/(3*deta);
//     const double dphiydeta = ((phiy.getElem(i+2*dy,j) - 4*phiy.getElem(i+dy,j)) + 3*phiy.getElem(i,j))/(3*deta);
//     const double detJ = dphixdxi * dphiydeta - dphixdeta*dphiydxi;
//     std::cout << "Element: [" << i << "," << j << "] direction: center[" << dx << "," << dy << "] phiy.getElem: " << phiy.getElem(i+1,j) << std::endl;
//     std::cout<< "dxi: " << dxi <<"dphixdxi: " << dphixdxi << " dphiydeta: " << dphiydeta << " dphixdeta: " << dphixdeta << " dphiydxi: "<< dphiydxi  << " detJ: " << detJ << std::endl;
//     if (detJ == 0){
//         std::cout << "Matrix is singular, cannot calculate derivative" << std::endl;
//         return(0);
//     }
//     return(1/detJ * (dx*(dudxi*dphiydeta - dudeta*dphiydxi)+dy*(dudeta*dphiydxi - dudxi*dphiydeta)));
// }

#endif
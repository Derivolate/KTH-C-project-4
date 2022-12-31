#ifndef GRIDFUNCTION_HPP
#define GRIDFUNCTION_HPP
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

        GridFunction& operator+=(const GridFunction&) const;  //Pointwise grid addition assignment operator
        // GridFunction& operator+(const GridFunction&) const;  //Pointwise grid addition operator
        GridFunction& operator*=(const GridFunction&) const;  //Pointwise grid multiplication assignment operator
        // GridFunction& operator*(const GridFunction&) const;  //Pointwise grid multiplication operator
        
        GridFunction Dx();    //X derivative
        GridFunction Dy();    //Y derivative
        GridFunction Grad();  //Gradient function
        
        static GridFunction& Grad(const GridFunction&, const GridFunction&); //Overloaded option for the gradient in case the derivatives in X and Y are already available
        void printFkt(std::string) const;
        void fillGrid();
    private:
        std::shared_ptr<Domain> domain;
        Matrix u;
        Matrix phix, phiy;
        double u_function(Point<double>);

        template<int dx, int dy> double center_(int i, int j);
        template<int dx, int dy> double assym_(int i, int j);
        // TODO template<int dx, int dy> double ghost_(int i, int j);

}; 

// ----------------------------------

// template<int dx, int dy> double GridFunction::center_(int i, int j, Matrix const & fnctn) {
// }

// template<int dx> double GridFunction::center_(int i, int j, Matrix const & fnctn) {
// }

template<int dx, int dy> double GridFunction::center_(int i, int j){
    if (dx*dy != 0){
        std::cout << "Only x or y may be given as direction. No inbetween" << std::endl;
        return(0);
    }
    const double deta = 1.0/(domain->getSize().X()-1);
    const double dxi = 1.0/(domain->getSize().Y()-1);
    const double dudeta = (u.getElem(i+1,j) - u.getElem(i-1,j))*0.5/deta;
    const double dudxi = (u.getElem(i,j+1) - u.getElem(i,j-1))*0.5/dxi;
    const double dphixdxi = (phix.getElem(i,j+1) - phix.getElem(i,j-1))*0.5/dxi;
    const double dphiydxi = (phiy.getElem(i,j+1) - phiy.getElem(i,j-1))*0.5/dxi;
    const double dphixdeta = (phix.getElem(i+1,j) - phix.getElem(i-1,j))*0.5/deta;
    const double dphiydeta = (phiy.getElem(i+1,j) - phiy.getElem(i-1,j))*0.5/deta;
    const double detJ = dphixdxi * dphiydeta - dphixdeta*dphiydxi;
    std::cout << "Element: [" << i << "," << j << "] direction: center[" << dx << "," << dy << "] phiy.getElem: " << phiy.getElem(i+1,j) << std::endl;
    std::cout<< "dxi: " << dxi <<"dphixdxi: " << dphixdxi << " dphiydeta: " << dphiydeta << " dphixdeta: " << dphixdeta << " dphiydxi: "<< dphiydxi  << " detJ: " << detJ << std::endl;
    if (detJ == 0){
        std::cout << "Matrix is singular, cannot calculate derivative" << std::endl;
        return(0);
    }
    return((1/detJ) * (dx*(dudxi*dphiydeta - dudeta*dphiydxi)+dy*(dudeta*dphiydxi - dudxi*dphiydeta)));
}

template<int dx, int dy> double GridFunction::assym_(int i, int j){
    if (dx*dy != 0){
        std::cout << "Only x or y may be given as direction. No inbetween" << std::endl;
        return(0);
    }
    const double deta = 1.0/(domain->getSize().X()-1);
    const double dxi = 1.0/(domain->getSize().Y()-1);
    const double dudeta = ((u.getElem(i+2*dy,j) - 4*u.getElem(i+dy,j)) + 3*u.getElem(i,j))/(3*deta);
    const double dudxi = ((u.getElem(i,j+2*dx) - 4*u.getElem(i,j+dx)) + 3*u.getElem(i,j))/(3*dxi);
    const double dphixdxi = ((phix.getElem(i,j+2*dx) - 4*phix.getElem(i,j+dx)) + 3*phix.getElem(i,j))/(3*dxi);
    const double dphiydxi = ((phiy.getElem(i,j+2*dx) - 4*phiy.getElem(i,j+dx)) + 3*phiy.getElem(i,j))/(3*dxi);
    const double dphixdeta = ((phix.getElem(i+2*dy,j) - 4*phix.getElem(i+dy,j)) + 3*phix.getElem(i,j))/(3*deta);
    const double dphiydeta = ((phiy.getElem(i+2*dy,j) - 4*phiy.getElem(i+dy,j)) + 3*phiy.getElem(i,j))/(3*deta);
    const double detJ = dphixdxi * dphiydeta - dphixdeta*dphiydxi;
    std::cout << "Element: [" << i << "," << j << "] direction: center[" << dx << "," << dy << "] phiy.getElem: " << phiy.getElem(i+1,j) << std::endl;
    std::cout<< "dxi: " << dxi <<"dphixdxi: " << dphixdxi << " dphiydeta: " << dphiydeta << " dphixdeta: " << dphixdeta << " dphiydxi: "<< dphiydxi  << " detJ: " << detJ << std::endl;
    if (detJ == 0){
        std::cout << "Matrix is singular, cannot calculate derivative" << std::endl;
        return(0);
    }
    return(1/detJ * (dx*(dudxi*dphiydeta - dudeta*dphiydxi)+dy*(dudeta*dphiydxi - dudxi*dphiydeta)));
}

#endif
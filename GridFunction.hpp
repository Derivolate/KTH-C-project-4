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
        
        void Dx();    //X derivative
        void Dy();    //Y derivative
        void Grad();  //Gradient function
        
        static GridFunction& Grad(const GridFunction&, const GridFunction&); //Overloaded option for the gradient in case the derivatives in X and Y are already available
        void printFkt();
        void fillGrid();
    private:
        std::shared_ptr<Domain> domain;
        Matrix u;
        double u_function(Point<double>);
}; 
#endif
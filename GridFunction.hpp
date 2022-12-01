#ifndef GRIDFUNCTION_HPP
#define GRIDFUNCTION_HPP
#include "Matrix.hpp"
#include "Domain.hpp"
#include "GridFunction.cpp"
class GridFunction{
    public:
        
        GridFunction(const GridFunction&);
        GridFunction(Domain);
        virtual ~GridFunction();

        void fill_matrix(); //

        GridFunction& operator=(const GridFunction&);   //Assignment operator
        GridFunction operator+(const GridFunction&) const;  //Pointwise grid addition operator
        GridFunction operator*(const GridFunction&) const;  //Pointwise grid multiplication operator
        GridFunction Dx() const;    //X derivative
        GridFunction Dy() const;    //Y derivative
        GridFunction Grad() const;  //Gradient function
        GridFunction Grad(GridFunction*, GridFunction*) const; //Overloaded option for the gradient in case the derivatives in X and Y are already available
        
    private:
        Matrix u;
        Domain *grid; //TODO implement shared_pointer
        virtual void u_function(double, double) = 0;
}; 
#endif
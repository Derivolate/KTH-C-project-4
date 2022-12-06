#ifndef GRIDFUNCTION_HPP
#define GRIDFUNCTION_HPP
#include "Matrix.hpp"
#include "Domain.hpp"

class GridFunction {
    public:
        GridFunction(const GridFunction&);
        GridFunction(Domain*);
        virtual ~GridFunction();

        // void fill_matrix(); //fills Matrix u with values from u_function 
        // We don't need this if the constructor already fills up u

        GridFunction& operator=(const GridFunction&);   //Assignment operator
        GridFunction& operator+=(const GridFunction&) const;  //Pointwise grid addition assignment operator
        // GridFunction& operator+(const GridFunction&) const;  //Pointwise grid addition operator
        GridFunction& operator*=(const GridFunction&) const;  //Pointwise grid multiplication assignment operator
        // GridFunction& operator*(const GridFunction&) const;  //Pointwise grid multiplication operator
        
        GridFunction& Dx() const;    //X derivative
        GridFunction& Dy() const;    //Y derivative
        GridFunction& Grad() const;  //Gradient function
        
        static GridFunction& Grad(const GridFunction&, const GridFunction&); //Overloaded option for the gradient in case the derivatives in X and Y are already available
        void printFkt();

    protected:
        void fillGrid();

    private:
        Matrix u;
        Domain *grid; //TODO implement shared_pointer
        virtual double u_function(double, double) = 0;
}; 
#endif
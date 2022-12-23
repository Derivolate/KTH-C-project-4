#ifndef MATRIX_HPP
#define MATRIX_HPP
#include <iostream>
#include "Point.hpp"
using namespace std;
class Matrix {
    private:
        double *elems;
        Point<int> size;
        int getm() const;
        int getn() const;
    public:
        Matrix();
        Matrix(Point<int>);
        Matrix(int, int = -1);                  //Constructor of empty matrix with dimension m by m
        Matrix(const Matrix&);                  //Copy constrctor
        ~Matrix();

        Matrix& operator=(const Matrix&);       //Assignment operator
        Matrix& operator+=(const Matrix&);      //Element wise addition Assignment
        Matrix operator+(const Matrix&)const;   //Element wise addition
        Matrix& operator*=(const Matrix&);      //Matrix product assignment
        Matrix& operator%=(const Matrix&);      //Matrix element-wise product assignment
        Matrix& operator*=(const double);       //Scalar product assignment
        Matrix operator*(const double)const;    //Matrix product
        Matrix Mexp(int)const;                  //Matrix exponent function
        double operator()(int, int) const;      //Indexing
        void setElem(double,int,int); //et an element
        double getElem(int, int) const;
        double norm() const; //Two-norm
        double* getMat() const;
        static double norm(const Matrix&);
        void printMatrix() const;
        void fillMatrix(double); //Fill the matrix with random numbers
        Point<int> getSize() const;
};
#endif
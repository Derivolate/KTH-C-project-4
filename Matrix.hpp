#ifndef MATRIX_HPP
#define MATRIX_HPP
#include <iostream>
using namespace std;
class Matrix {
    private:
        double *elems;
        int size;
        // static int factorial(const int);
    public:
        Matrix();
        Matrix(int m);                          //Constructor of empty matrix with dimension m by m
        Matrix(const Matrix&);                  //Copy constrctor
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
        int getSize() const;
};
#include "Matrix.cpp"
#endif
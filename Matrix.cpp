#include <cmath>
#include <stdexcept>
#include "Matrix.hpp"

Matrix::Matrix() : size(Point<int>(-1,-1)),elems(nullptr){}

Matrix::Matrix(int m, int n=-1) : Matrix(Point<int>(m,n!=-1?n:m)){} //Initialise a square matrix if no value of n is provided

Matrix::~Matrix(){
    if(elems != nullptr)
    {
        delete[]elems;
    }
}

Matrix::Matrix(Point<int> _size) :size(_size) {//Initialise an identity matrix
    elems = new double[size.X()*size.Y()];
    for (int j(0); j < size.Y(); ++j){
        for (int i(0); i < size.X(); ++i){
            if (i==j){
                elems[i+j*size.X()] = 1;
            }
            else{
                elems[i+j*size.X()] = 0;
            }
        }
    }
}

Matrix::Matrix(const Matrix& M): size(M.size) { 
    elems = new double[size.Y()*size.X()];
    for (int j(0); j < size.Y(); ++j)
        for (int i(0); i < size.X(); ++i)
            setElem(M(i,j),i,j);
}

Matrix& Matrix::operator=(const Matrix& M) {
    if (this != &M) {
        size = M.size;
        elems = new double[size.Y()*size.X()];
        for (int j(0); j < size.Y(); ++j)
            for (int i(0); i < size.X(); ++i)
                setElem(M(i,j),i,j);
    }
    return *this; // dereferencing!
}

Matrix& Matrix::operator+=(const Matrix& M) { 
    if(size==M.size){
        for (int j(0); j < size.Y(); ++j)
            for (int i(0); i < size.X(); ++i)
                setElem(M(i,j)+getElem(i,j),i,j);
        return *this;
    }else{
        throw std::invalid_argument("Matrices of different sizes cannot be summed");
    }
}

Matrix Matrix::operator+(const Matrix& M) const{
    Matrix N = Matrix(size);
    return N+=M;    
}

Matrix& Matrix::operator%=(const Matrix& M) {
    if(size==M.size){
        for (int j(0); j < size.Y(); ++j)
            for (int i(0); i < size.X(); ++i)
                setElem(getElem(i,j)*M(i,j),i,j);
        return *this;
    }else{
        throw std::invalid_argument("Matrices of different sizes cannot be element-wise multiplied"); 
    }
}

Matrix& Matrix::operator*=(const double a) {
    for (int j(0); j < size.Y(); ++j)
        for (int i(0); i < size.X(); ++i)
            setElem(getElem(i,j)*a,i,j);
    return *this; // dereferencing!
}

Matrix Matrix::operator*(const double a) const{
    Matrix N = Matrix(size);
    return N*=a;
}

double Matrix::getElem(int i, int j) const{
    return elems[i+size.X()*j];
}

double Matrix::operator()(int i, int j) const{ //retrieve element
    return getElem(i,j);
}

double Matrix::norm() const{
    double norm(0);
    for (int j(0); j < size.Y(); ++j)
        for (int i(0); i < size.X(); ++i)
            norm+= getElem(i,j)*getElem(i,j);
    return std::sqrt(norm);
}

double Matrix::norm(const Matrix& M){
    return M.norm();
}

void Matrix::setElem(double val ,int i, int j){
    elems[i+j*size.X()] = val;
}

void Matrix::printMatrix() const { 
    cout << size.X() << ' ' << size.Y() << endl;
    for (int i(0); i < size.X(); ++i){
        for (int j(0); j < size.Y(); ++j){
            cout << elems[i + j*size.X()] <<" ";
        }
        cout << endl;
    }
    cout << endl;

}

void Matrix::fillMatrix(double maxVal = 10) {
    for (int j(0); j < size.Y(); ++j)
        for (int i(0); i < size.X(); ++i)
            setElem(((double)rand()/RAND_MAX)*maxVal, i, j);
}

double* Matrix::getMat() const{
    return elems;
}

Point<int> Matrix::getSize() const{
    return size;
}
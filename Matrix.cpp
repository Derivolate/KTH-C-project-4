#include <cmath>
#include <stdexcept>
#include "Matrix.hpp"

Matrix::Matrix() : size(Point<int>(-1,-1)),elems(nullptr){}

Matrix::Matrix(int m, int n) :Matrix(Point<int>(m,n!=-1?n:m)) {}

Matrix::~Matrix(){
    if(getm()>0 || getn()>0)
    {
        delete[]elems;
    }
}

Matrix::Matrix(Point<int> _size) :size(_size) {//Initialise an identity matrix
    elems = new double[getm()*getn()];
    for (int j(0); j < getn(); ++j){
        for (int i(0); i < getm(); ++i){
            if (i==j){
                elems[i+j*getm()] = 1;
            }
            else{
                elems[i+j*getm()] = 0;
            }
        }
    }
}

Matrix::Matrix(const Matrix& M): size(M.size) { 
    elems = new double[getn()*getm()];
    for (int j(0); j < getn(); ++j)
        for (int i(0); i < getm(); ++i)
            setElem(M(i,j),i,j);
}

Matrix& Matrix::operator=(const Matrix& M) {
    if (this != &M) {
        size = M.size;
        elems = new double[getn()*getm()];
        for (int j(0); j < getn(); ++j)
            for (int i(0); i < getm(); ++i)
                setElem(M(i,j),i,j);
    }
    return *this; // dereferencing!
}

Matrix& Matrix::operator+=(const Matrix& M) {
    if(size==M.size){
        for (int j(0); j < getn(); ++j)
            for (int i(0); i < getm(); ++i)
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

// Matrix& Matrix::operator*=(const Matrix& M) {
//     Matrix N = Matrix(size);
//     double elem(0);
//     if(size==M.size){
//         for (int i = 0; i < getn(); i++){
//             for (int j = 0; j < getm(); j++){
//                 for (int k = 0; k < size; k++){
//                     elem += getElem(i,k)*M(k,j);
//                 }
//                 N.setElem(elem,i,j);
//                 elem = 0;
//             }
//         }
//         *this = N;
//         return *this;
//     }else{
//         throw std::invalid_argument("Multiplication of different sized matrices is not implemented");
//     }
// }

Matrix& Matrix::operator%=(const Matrix& M) {
    if(size==M.size){
        for (int j(0); j < getn(); ++j)
            for (int i(0); i < getm(); ++i)
                setElem(getElem(i,j)*M(i,j),i,j);
        return *this;
    }else{
        throw std::invalid_argument("Matrices of different sizes cannot be element-wise multiplied"); 
    }
}

Matrix& Matrix::operator*=(const double a) {
    for (int j(0); j < getn(); ++j)
        for (int i(0); i < getm(); ++i)
            setElem(getElem(i,j)*a,i,j);
    return *this; // dereferencing!
}

Matrix Matrix::operator*(const double a) const{
    Matrix N = Matrix(size);
    return N*=a;
}

double Matrix::getElem(int i, int j) const{
    return elems[i+getm()*j];
}

double Matrix::operator()(int i, int j) const{ //retrieve element
    return getElem(i,j);
}

double Matrix::norm() const{
    double norm(0);
    for (int j(0); j < getn(); ++j)
        for (int i(0); i < getm(); ++i)
            norm+= getElem(i,j)*getElem(i,j);
    return std::sqrt(norm);
}

double Matrix::norm(const Matrix& M){
    return M.norm();
}

void Matrix::setElem(double val ,int i, int j){
    elems[i+j*getm()] = val;
}

void Matrix::printMatrix() const { 
    cout << endl;
    for (int i(0); i < getm(); ++i){
        for (int j(0); j < getn(); ++j){
            cout << elems[i + j*getm()] <<" ";
        }
        cout << endl;
    }
    cout << endl;

}

void Matrix::fillMatrix(double maxVal = 10) {
    for (int j(0); j < getn(); ++j)
        for (int i(0); i < getm(); ++i)
            setElem(((double)rand()/RAND_MAX)*maxVal, i, j);
}

double* Matrix::getMat() const{
    return elems;
}

Point<int> Matrix::getSize() const{
    return size;
}

int Matrix::getm() const{
    return (int)size.Y();
}

int Matrix::getn() const{
    return (int)size.X();
}

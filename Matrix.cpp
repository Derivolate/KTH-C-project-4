#include <cmath>
#include <stdexcept>
#include "Matrix.hpp"

Matrix::Matrix() : size(-1){}

Matrix::Matrix(int m) :size(m) {//Initialise an identity matrix
    elems = new double[size*size];
    for (int i = 0; i < size; i++) {//per row
        for (int j = 0; j < size; j++) {//per column
            if (i==j){
                elems[i+j*size] = 1;
            }
            else{
                elems[i+j*size] = 0;
            }
        }
    }
}

Matrix::Matrix(const Matrix& M): size(M.size) { 
    elems = new double[size*size];
    for (int i = 0; i < size; i++) //per row
        for (int j = 0; j < size; j++) //per column
            setElem(M(i,j),i,j);
}

Matrix& Matrix::operator=(const Matrix& M) {
    if (this != &M) {
        size = M.size;
        elems = new double[size*size];
        for (int i = 0; i < size; i++) //per row
            for (int j = 0; j < size; j++) //per column
                setElem(M(i,j),i,j);
    }
    return *this; // dereferencing!
}

Matrix& Matrix::operator+=(const Matrix& M) {
    if(size==M.size){
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
                setElem(M(i,j)+getElem(i,j),i,j);
        return *this;
    }else{
        throw std::invalid_argument("Matrices of different sizes cannot be summed");
    }
}

Matrix Matrix::operator+(const Matrix& M) const{
    Matrix N = Matrix(size);
    if(size==M.size){
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
                N.setElem(getElem(i,j)+M(i,j),i,j);
        return N;
    }else{
        throw std::invalid_argument("Matrices of different sizes cannot be summed");
    }
    
}

Matrix& Matrix::operator*=(const Matrix& M) {
    Matrix N = Matrix(size);
    double elem(0);
    if(size==M.size){
        for (int i = 0; i < size; i++){
            for (int j = 0; j < size; j++){
                for (int k = 0; k < size; k++){
                    elem += getElem(i,k)*M(k,j);
                }
                N.setElem(elem,i,j);
                elem = 0;
            }
        }
        *this = N;
        return *this;
    }else{
        throw std::invalid_argument("Multiplication of different sized matrices is not implemented");
    }
}

Matrix& Matrix::operator%=(const Matrix& M) {
    if(size==M.size){
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
                setElem(getElem(i,j)*M(i,j),i,j);
        return *this;
    }else{
        throw std::invalid_argument("Matrices of different sizes cannot be element-wise multiplied"); 
    }
}

Matrix& Matrix::operator*=(const double a) {
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            setElem(getElem(i,j)*a,i,j);
    return *this; // dereferencing!
}

Matrix Matrix::operator*(const double a) const{
    Matrix N = Matrix(size);
    for (int i = 0; i < size; i++) 
        for (int j = 0; j < size; j++)
            N.setElem(getElem(i,j)*a,i,j);
    return N; // dereferencing!
    //TODO lhs as well?
}

Matrix Matrix::Mexp(int tol=18)const{
    //Pick a scaling factor and series length based on the tolerance and norm from the pre-calculated table
    //Source: Nineteen Dubious Ways to Compute the Exponential of a Matrix, Twenty-Five Years Later; Cleve Moler and Charles Van Loan
    double norm2 = norm(); // Get the 2-norm
    int series_length;
    int scaling_pow;
    if (tol < 4){
        if(norm2 < 0.01){
            series_length = 1;
            scaling_pow = 0;
        }else if(norm2 < 0.1){
            series_length = 3;
            scaling_pow = 0;
        }else if(norm2 < 1){
            series_length = 5;
            scaling_pow = 1;
        }else if(norm2 < 10){
            series_length = 4;
            scaling_pow = 5;
        }else if(norm2 < 100){
            series_length = 4;
            scaling_pow = 8;
        }else{
            series_length = 5;
            scaling_pow = 11;
        }
    }else if(tol < 7){
        if(norm2 < 0.01){
            series_length = 2;
            scaling_pow = 1;
        }else if(norm2 < 0.1){
            series_length = 4;
            scaling_pow = 0;
        }else if(norm2 < 1){
            series_length = 7;
            scaling_pow = 1;
        }else if(norm2 < 10){
            series_length = 6;
            scaling_pow = 5;
        }else if(norm2 < 100){
            series_length = 5;
            scaling_pow = 9;
        }else{
            series_length = 7;
            scaling_pow = 11;
        }
    }else if(tol < 10){
        if(norm2 < 0.01){
            series_length = 3;
            scaling_pow = 1;
        }else if(norm2 < 0.1){
            series_length = 4;
            scaling_pow = 2;
        }else if(norm2 < 1){
            series_length = 6;
            scaling_pow = 3;
        }else if(norm2 < 10){
            series_length = 8;
            scaling_pow = 5;
        }else if(norm2 < 100){
            series_length = 7;
            scaling_pow = 9;
        }else{
            series_length = 6;
            scaling_pow = 13;
        }
    }else if(tol < 13){
        if(norm2 < 0.01){
            series_length = 4;
            scaling_pow = 1;
        }else if(norm2 < 0.1){
            series_length = 4;
            scaling_pow = 4;
        }else if(norm2 < 1){
            series_length = 8;
            scaling_pow = 3;
        }else if(norm2 < 10){
            series_length = 7;
            scaling_pow = 7;
        }else if(norm2 < 100){
            series_length = 9;
            scaling_pow = 9;
        }else{
            series_length = 8;
            scaling_pow = 13;
        }
    }else{
        if(norm2 < 0.01){
            series_length = 5;
            scaling_pow = 1;
        }else if(norm2 < 0.1){
            series_length = 5;
            scaling_pow = 4;
        }else if(norm2 < 1){
            series_length = 7;
            scaling_pow = 5;
        }else if(norm2 < 10){
            series_length = 9;
            scaling_pow = 7;
        }else if(norm2 < 100){
            series_length = 10;
            scaling_pow = 10;
        }else{
            series_length = 8;
            scaling_pow = 14;
        }
    }
    
    double factorial = 1; 
    Matrix A = *this;
    series_length += 5;
    scaling_pow+=5;
    A *= 1/pow(2,scaling_pow);
    Matrix scaled = A;                      //Temporary storage for matrix multiplications
    Matrix Mexp(size);                      //Initialized as identity matrix, will be returned as result
    Mexp = Mexp+A*(1/factorial);            //Taking this out of the loop saves one matrix multiplication
    for(int i = 2; i<=series_length;i++){
        factorial*=i;
        A%=scaled;
        Mexp = Mexp+ A*(1/factorial);
    }
    for(int i = 1; i<= scaling_pow;i++){
        Mexp*=Mexp;
    }
    return Mexp;
}

double Matrix::operator()(int i, int j) const{ //retrieve element
    return getElem(i,j);
}

double Matrix::getElem(int i, int j) const{
    return elems[i+j*size];
}

double Matrix::norm() const{
    double norm(0);
    for (int i = 0; i < size; i++) 
        for (int j = 0; j < size; j++)
            norm+= getElem(i,j)*getElem(i,j);
    return std::sqrt(norm);
}

double Matrix::norm(const Matrix& M){
    return M.norm();
}

void Matrix::setElem(double val ,int i, int j){
    elems[i+j*size] = val;
}

void Matrix::printMatrix() const { 
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            cout << elems[i + j*size] <<" ";
        }
        cout << endl;
    }
}

void Matrix::fillMatrix(double maxVal = 10) {
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            setElem(((double)rand()/RAND_MAX)*maxVal, i, j);
}

double* Matrix::getMat() const{
    return elems;
}

int Matrix::getSize() const{
    return size;
}
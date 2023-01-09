#include "Matrix.hpp"
#include "Domain.hpp"
#include "GridFunction.hpp"
#include <iostream>
#include <fstream>
#include <cmath>

GridFunction::GridFunction(std::shared_ptr<Domain> _dom) : domain(_dom),u(_dom->getSize()), phix(_dom->getSize()),phiy(_dom->getSize()){}

GridFunction::GridFunction(GridFunction const & B) : domain(B.domain), u(B.u), phix(B.phix), phiy(B.phiy){}

GridFunction& GridFunction::operator+=(const GridFunction& B) {
    u += B.u;
    return *this;
}

GridFunction GridFunction::operator+(const GridFunction& B) {
    GridFunction A(B);
    A.u += u;
    return A;
}

void GridFunction::printFkt(std::string file) const{
    u.printMatrix();
    std::string path = "../MatlabGrid/" + file;
    const char* str = path.c_str();
    double size = domain->getSize().X()*domain->getSize().Y();
	FILE *fu;
	fu = fopen(str,"wb");
	fwrite(u.getMat(),sizeof(double),size,fu);
	fclose(fu);
}

// This functionality cannot be in the constructor because it calls a virtual function
void GridFunction::fillGrid(){
    int const m(u.getSize().X()), n(u.getSize().Y());
    Point<double>* coords = domain->getGrid();
    Point<double> coord;
    std::cout << m << " " << n << std::endl;
    double val;
    for(int i(0);i<m;++i){
        for(int j(0);j<n;++j){
            coord = coords[i+m*j];
            val = u_function(coord);
            u.setElem(val,i,j);
            phix.setElem(coord.X(),i,j);
            phiy.setElem(coord.Y(),i,j);
        }
    }
    phix.printMatrix();
    phiy.printMatrix();
}

// Derivative in X direction returning new GridFunction Object with Derivatives
GridFunction GridFunction::Dx() {
    int m(u.getSize().X()), n(u.getSize().Y());
    GridFunction grdfktn(*this);
    // corners
    grdfktn.u.setElem(vertex_<1,0,1,1>(0,0),0,0);
    grdfktn.u.setElem(vertex_<1,0,-1,1>(0,n-1),0,n-1);
    grdfktn.u.setElem(vertex_<1,0,1,-1>(m-1,0),m-1,0);
    grdfktn.u.setElem(vertex_<1,0,-1,-1>(m-1,n-1),m-1,n-1);
    // boundary 1 and 3
    for(int i(1); i < m-1; i++){
        grdfktn.u.setElem(edge_<1,0,1,0>(i,0),i,0); 
        grdfktn.u.setElem(edge_<1,0,-1,0>(i,n-1),i,n-1); 
    }
    for(int j(1);j<n-1;j++){
        grdfktn.u.setElem(edge_<1,0,0,1>(0,j),0,j);
        grdfktn.u.setElem(edge_<1,0,0,-1>(m-1,j),m-1,j);
    }
    for(int i(1); i < m-1; i++){  
              for(int j(1); j < n-1; j++){
                    grdfktn.u.setElem(face_<1,0>(i,j),i,j);
        }
    }
    return(grdfktn);
}

GridFunction GridFunction::Dy() {
    int m(u.getSize().X()), n(u.getSize().Y());
    GridFunction grdfktn(*this);
    // corners
    grdfktn.u.setElem(vertex_<0,1,1,1>(0,0),0,0);
    grdfktn.u.setElem(vertex_<0,1,-1,1>(0,n-1),0,n-1);
    grdfktn.u.setElem(vertex_<0,1,1,-1>(m-1,0),m-1,0);
    grdfktn.u.setElem(vertex_<0,1,-1,-1>(m-1,n-1),m-1,n-1);
    // boundary 1 and 3
    for(int i(1); i<m-1; i++){
        grdfktn.u.setElem(edge_<0,1,1,0>(i,0),i,0); 
        grdfktn.u.setElem(edge_<0,1,-1,0>(i,n-1),i,n-1); 
    }
    for(int j(1);j<n-1;j++){
        grdfktn.u.setElem(edge_<0,1,0,1>(0,j),0,j);
        grdfktn.u.setElem(edge_<0,1,0,-1>(m-1,j),m-1,j);
    }
    for(int i(1); i < m-1; i++){  
        for(int j(1); j < n-1; j++){
            grdfktn.u.setElem(face_<0,1>(i,j),i,j);
        }
    }
    return(grdfktn);
}

// TODO: Implement Stensil Calculation of Laplacian. For now this suffices
GridFunction GridFunction::DDxy(GridFunction* Dx, GridFunction* Dy) {
    // TODO Write check if Dx and Dy come from the same domain and original gridfunction
    return(Dx->Dx() + Dy->Dy());
}

double GridFunction::u_function(Point<double> p){
    double x(p.X()),y(p.Y());
    return(sin(pow(x/10,2))*cos(x/10)+y);
}
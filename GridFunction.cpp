#include "Matrix.hpp"
#include "Domain.hpp"
#include "GridFunction.hpp"
#include <iostream>
#include <fstream>
#include <cmath>

GridFunction::GridFunction(std::shared_ptr<Domain> _dom) : domain(_dom),u(_dom->getSize()), phix(_dom->getSize()),phiy(_dom->getSize()){}

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
    Point<int> size =  domain->getSize();
    Point<double>* coords = domain->getGrid();
    Point<double> coord;
    double val;
    for(int i(0);i<u.getSize().X();++i){
        for(int j(0);j<u.getSize().Y();++j){
            coord = coords[i+u.getSize().X()*j];
            val = u_function(coord);
            u.setElem(val,i,j);
            phix.setElem(coord.X(),i,j);
            phiy.setElem(coord.Y(),i,j);
        }
    }
    std::cout << "------- PHIX AND PHIY -------" << std::endl;
    phix.printMatrix();
    phiy.printMatrix();
}

GridFunction GridFunction::Dx() {
    Point<double>* grid = domain->getGrid();
    int m(u.getSize().X()), n(u.getSize().Y());
    GridFunction grdfktn(domain);
    // corners
    std::cout<<"====== CORNERS DX ======" << std::endl;
    grdfktn.u.setElem(vertex_<1,0,1,1>(0,0),0,0);
    grdfktn.u.setElem(vertex_<-1,0,-1,1>(0,n),0,n);
    grdfktn.u.setElem(vertex_<1,0,1,-1>(n,0),m,0);
    grdfktn.u.setElem(vertex_<-1,0,-1,-1>(n,m),n,m);
    // boundary 1 and 3
    std::cout<<"====== BOUNDARIES ======" << std::endl;
    for(int i(1); i < m-1; i++){
        grdfktn.u.setElem(edge_<1,0,1,0>(i,0),i,0);
        grdfktn.u.setElem(edge_<1,0,-1,0>(i,n),i,n);
    }
    for(int j(1);j<n-1;j++){
        grdfktn.u.setElem(edge_<1,0,0,1>(0,j),0,j);
        grdfktn.u.setElem(edge_<-1,0,0,-1>(m,j),m,j);
    }
    std::cout<<"====== REST DX ======" << std::endl;
    for(int i(0); i < m; i++){  
              for(int j(1); j < n-1; j++){
                    grdfktn.u.setElem(face_<1,0>(i,j),i,j);
        }
    }
    return(grdfktn);
}

GridFunction GridFunction::Dy() {
    Point<double>* grid = domain->getGrid();
    int m(u.getSize().X()), n(u.getSize().Y());
    GridFunction grdfktn(domain);
    // corners
    std::cout<<"====== CORNERS DY ======" << std::endl;
    grdfktn.u.setElem(vertex_<0,1,1,1>(0,0),0,0);
    grdfktn.u.setElem(vertex_<0,1,-1,1>(m,0),m,0);
    grdfktn.u.setElem(vertex_<0,1,1,-1>(0,n),0,n);
    grdfktn.u.setElem(vertex_<0,1,-1,-1>(n,m),n,m);
    // boundary 0 and 2
    std::cout<<"====== BOUNDARIES ======" << std::endl;
       for(int i(1); i < m-1; i++){
        grdfktn.u.setElem(edge_<0,1,1,0>(i,0),i,0);
        grdfktn.u.setElem(edge_<0,1,-1,0>(i,n),i,n);
    }
    for(int j(1);j<n-1;j++){
        grdfktn.u.setElem(edge_<0,1,0,1>(0,j),0,j);
        grdfktn.u.setElem(edge_<0,1,0,-1>(m,j),m,j);
    }
    std::cout<<"====== REST DY ======" << std::endl;
    for(int i(0); i < m-1; i++){  
              for(int j(1); j < n; j++){
                    grdfktn.u.setElem(face_<0,1>(i,j),i,j);
        }
    }
    return(grdfktn);
}

double GridFunction::u_function(Point<double> p){
    double x(p.X()),y(p.Y());
    // return(sin(pow(x/10,2)*cos(x/10)+y));
    return(x);
}
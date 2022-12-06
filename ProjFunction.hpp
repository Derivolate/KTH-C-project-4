#ifndef PROJ_FUNCTION_HPP
#define PROJ_FUNCTION_HPP
#include "GridFunction.hpp"
#include "Domain.hpp"

class ProjFunction : public GridFunction {
    public:
        ProjFunction(Domain*);
        // ~ProjFunction();
    private:
        double u_function(double, double);
};
#endif
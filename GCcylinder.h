#ifndef GCCYLINDER_H
#define GCCYLINDER_H

#include "common.h"

//! Constraint Function Template
unsigned onCylinder_Z(double radius,unsigned index,LINCS &lincs)
//unsigned onCylinder (radius, index, lincs);
//The cylinder equation: sqrt(x*x + y*y) - radius = 0;
//where x,y are the components of the position vector of the arbitrary point in the surface.
//The function is about confining bead index of the chain of the lincs on the surface,
//and return the order of the confinement which will be used when calculating the constraining force.
{
    Vector3d grad = lincs.R[index];
    grad[2] = 0;
    const double len = grad.norm();
    const double g = len - radius;
    grad /= len;
    
    lincs.setGCerror(g);
    lincs.setGCdirection(index,grad);
    return lincs.endInsertOneGC();
}

#endif // GCCYLINDER_H

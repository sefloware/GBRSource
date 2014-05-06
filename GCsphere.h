#ifndef GCSPHERE_H
#define GCSPHERE_H

#include "common.h"

unsigned onSphere(double radius,unsigned index, LINCS &lincs)
//unsigned onSphere (radius, index, lincs);
//The sphere equation: r.norm() - radius = 0,
//where r is the position vector of the arbitrary point in the surface.
//and r.norm() means the norm of r.
//The function is about confining bead index of the chain of the lincs on the surface of the sphere,
//and return the order of the confinement which will be used when calculating the constraining force.
{
    const Vector3d &r = lincs.R[index];
    double rLen=r.norm();
    const double g = rLen -radius;
    Vector3d grad = r/rLen;
    
    lincs.setGCerror(g);
    lincs.setGCdirection(index,grad);
    return lincs.endInsertOneGC();
}

#endif // GCSPHERE_H

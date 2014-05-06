#ifndef RODLINK_H_
#define RODLINK_H_

#include "common.h"

unsigned rodLink(double length, unsigned i, unsigned j, LINCS &lincs)
//unsigned rodLink (length, i, j, lincs);
//The rod equation: (ri - rj).norm() - length = 0,
//where ri, rj are the position vectors of the two ends of the rod,
//and (ri - rj).norm() means the distance between the two ends.
//The function is about confining bead i and bead j of the chain of the lincs on the ends of the rod,
//and return the order of the confinement which will be used when calculating the constraining force.
{
    Vector3d tR =lincs.R[i]-lincs.R[j];
    double tRl = tR.norm();
    lincs.setGCerror(tRl-length);
    tR /= tRl;
    lincs.setGCdirection(i,tR);
    lincs.setGCdirection(j,-tR);
    return lincs.endInsertOneGC();
} 

unsigned rodLink(const Vector3d &r0,double length, unsigned index, LINCS &lincs)
//unsigned rodLink (r0, length, index, lincs);
//The rod equation: (r_index - r0).norm() - length = 0,
//where r_index is the position vectors of the bead index, r0 is the fixed point, 
//and (r_index - r0).norm() means the distance between the bead and the fixed point.
//The function is about confining bead index to the fixed point r0 with a rigid the rod length,
//and return the order of the confinement which will be used when calculating the constraining force.
{
    Vector3d tR =lincs.R[index]-r0;
    double tRl = tR.norm();
    lincs.setGCerror(tRl-length);
    tR /= tRl;
    lincs.setGCdirection(index,tR);
    return lincs.endInsertOneGC(); 
}

#endif // RODLINK_H_

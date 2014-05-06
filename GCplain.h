#ifndef GCPLAIN_H
#define GCPLAIN_H

#include "common.h"

unsigned onPlain (Vector3d n, double d, unsigned index, LINCS &lincs)
//unsigned onPlain (n, d, index, lincs);
//The plain equation: n*r + d = 0,
//where r is the position vector of the arbitrary point in the surface.
//The function is about confining bead index of the chain of the lincs on the surface,
//and return the order of the confinement which will be used when calculating the constraining force.
{
    lincs.setGCerror(n.dot(lincs.R[index])+d);
    lincs.setGCdirection(index,n);
    return lincs.endInsertOneGC();
}

#endif // GCPLAIN_H

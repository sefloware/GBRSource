#ifndef GCdistance_h
#define GCdistance_h
#include "common.h"

//! Constraint Function Template
unsigned seperateZdistance(double distance, unsigned i, unsigned j, LINCS &lincs)
{
    /*...*/
    double g=(lincs.R[j].z()-lincs.R[i].z())-distance;
    lincs.setGCerror(g);
    lincs.setGCdirection(j,Vector3d::UnitZ());
    lincs.setGCdirection(i,-Vector3d::UnitZ());
    return lincs.endInsertOneGC();
}

#endif


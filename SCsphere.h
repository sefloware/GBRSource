#ifndef SCSPHERE_H
#define SCSPHERE_H

#include "common.h"
#include "rebound_PB.h"

//! Roadblock Function Template
void inSphere(double radius, LINCS &lincs)
//void inSphere (radius, lincs);
//The sphere equation: r.norm() - radius = 0,
//where r is the position vector of the arbitrary point in the surface.
//and r.norm() means the norm of r.
//The function is about confining the chain of the lincs in the sphere.
//The rebounding function is rebound_PB() in file rebound_PB.h.
{
    const double sDT=std::sqrt(lincs.beadDiffusion()*lincs.timeStep);
    const double s5DT=std::sqrt(5.0)*sDT;

    for(int i=0; i<lincs.R.size(); ++i)
    {
        const Vector3d &r = lincs.R[i];
        //! @sphere{ f(r) = r.norm(); df(r) = r/f(r); f(r) > radius}
        //! #Assemble = sphere

        const double sphere_f = r.norm();
        
        if(radius-s5DT<sphere_f && sphere_f<radius)
            lincs.plusX(i, rebound_PB(sDT,lincs.generator,radius-sphere_f)*-r/sphere_f );
    }
}

#endif // SCSPHERE_H

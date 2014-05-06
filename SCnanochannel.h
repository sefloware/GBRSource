#ifndef SCNANOCHANNEL_H_
#define SCNANOCHANNEL_H_

#include "common.h"
#include "rebound_PB.h"

void inCylinder_Z(double radius, LINCS &lincs)
//void inCylinder_Z (radius, lincs);
//The cylinder equation: sqrt(x*x + y*y) - radius = 0;
//where x,y are the components of the position vector of the arbitrary point in the surface.
//The function is about confining the chain of the lincs in the cylinder.
//The rebounding function is rebound_PB() in file rebound_PB.h.
{
    const double sqrtDt=std::sqrt(lincs.beadDiffusion()*lincs.timeStep);

    for(int i=0;i<lincs.R.size();++i)
    {
        Vector3d R = lincs.R[i];
        R.z() = 0.0;
        double Rnorm = R.norm();
        double distance = radius - Rnorm;
        if(distance < std::sqrt(5.0)*sqrtDt)
        {
            Vector3d direction = R/(-Rnorm);

            lincs.plusX(i, rebound_PB(sqrtDt,lincs.generator,distance)*direction );
        }

    }
}

void inParallelPlate(unsigned component,double dist, LINCS &lincs)
//void inParallelPlate (component, dist, lincs);
//Argument component: 0->axis x; 1->axis y; 2->axis z.
//The parallel plate with distance dist and direction of axis component.
//The function is about confining the chain of the lincs between the plates.
//The rebounding function is rebound_PB() in file rebound_PB.h.
{
    const double sqrtDt=std::sqrt(lincs.beadDiffusion()*lincs.timeStep);

    for(int i=0;i<lincs.R.size();++i)
    {
        double tmp = lincs.R[i][component];
        double distance = dist/2-std::abs(tmp);
        if(distance < std::sqrt(5.0)*sqrtDt)
        {
            Vector3d direction = Vector3d::Unit(component)*( tmp>0.0 ? -1.0 : 1.0);
            lincs.plusX(i, rebound_PB(sqrtDt,lincs.generator,distance)*direction  );
        }
    }
}

void inSquareChannel_Z(double dist_x,double dist_y,LINCS &lincs)
//void inSquareChannel_Z (dist_x, dist_y, lincs);
//The square channel is constructed with z-direction and y-direction parallel plates.
//The function is about confining the chain of the lincs in the square channel.
//The rebounding function is rebound_PB() in file rebound_PB.h.
{
    inParallelPlate(0,dist_x,lincs);
    inParallelPlate(1,dist_y,lincs);
}

#endif // SCNANOCHANNEL_H_

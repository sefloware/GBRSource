#ifndef REBOUND_PB_H_
#define REBOUND_PB_H_

#include <cstdio>
#include <cmath>
#include <boost/math/special_functions/erf.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>
#include "common.h"

inline double rebound_PB(double sqtDt,base_generator_type &generator, double distance)
//double rebound_PB(sqtDt, generator, distance)
//Argument sqtDt is the square root of the multiplication between the diffusion of bead and the time step.
//Argument generator is the reference to the random generator.
//Argument distance is the distance to the refelcting wall.
//The function about the displacement of beads under the effect of the reflecting wall.
//Reference: Peters E, Barenbrug T M. Physical Review E, 2002, 66(5): 056701.
{
    distance/=sqtDt;
    double terf=boost::math::erf(distance/2.0);
    double f1=2.0/std::sqrt(M_PI)*std::exp(-distance*distance/4.0)-distance*(1.0-terf);
    double tmp = distance*terf+2.0/std::sqrt(M_PI)*std::exp(-distance*distance/4.0);
    double f2=std::sqrt(std::abs(2.0+distance*distance-tmp*tmp));
    boost::variate_generator<base_generator_type&,boost::uniform_smallint<> >
        DU01Gen(generator, boost::uniform_smallint<>(0,1));
    distance=std::abs(f1*sqtDt + f2*sqtDt*(DU01Gen()*2-1) );
    return distance;
}

#endif // REBOUND_PB_H_

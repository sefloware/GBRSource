#ifndef LANGEVIN_H
#define LANGEVIN_H
#include "common.h"
#include "rebound_PB.h"

Vector3d brownianMove(double step, double diffustion, base_generator_type &generator)
{
    typedef boost::random::normal_distribution<double> normal;        
    typedef boost::variate_generator<base_generator_type&, normal> gen_type;

    const double sigma = std::sqrt(2.0*step *diffustion);
    gen_type normalGen(generator, normal(0.0,sigma));
    
    return Vector3d(normalGen(),normalGen(),normalGen());
}

double interImpactForce(double epsilon, double sigma, double distance)
{
    if(distance > sigma || distance < sigma/1000.0)
        return 0.0;

    const double squareDist = distance * distance;
    const double squareSigma = sigma * sigma;
    const double squareRate = squareSigma/squareDist;

    const double term6 = squareRate * squareRate * squareRate;
    const double term12 = term6 * term6;
    return 24*epsilon/distance*(2*term12-term6);
}

Vector3d inSphere(double step, double diffusion, double beadRadius,const Vector3d &position, double sphereRadius, base_generator_type &generator)
{
    const double sDT=std::sqrt(diffusion*step);
    const double s5DT=std::sqrt(5.0)*sDT;

    //! @sphere{ f(r) = r.norm(); df(r) = r/f(r); f(r) > radius}
    //! #Assemble = sphere

    const double posNorm = position.norm();
    const double sphere_f = posNorm + beadRadius;

    if(sphereRadius-s5DT < sphere_f && posNorm > sphereRadius/1000.0)
        return rebound_PB(sDT,generator,sphereRadius-sphere_f) * (-position/posNorm);
    else
        return Vector3d::Zero();
}

#endif // LANGEVIN_H

#ifndef CHECKPARAMETER_H
#define CHECKPARAMETER_H
#include "common.h"

void checkRodlengthValidUnderStretching(double rodlength,double persistence,double T,double stretch)
//Reference: Wang J, Gao H. The Journal of chemical physics, 2005, 123(8): 084906-084906-13.
{
    double MaxRodlength = std::sqrt(M_KB*T*persistence/stretch);
    bool result=rodlength < MaxRodlength;
    if(!result)
    {
        std::cerr << "rod(under the streching force) should < " <<
            MaxRodlength << std::endl;
        exit(0);
    }
}
void checkRodlengthValidInCylinder(double rodlength,double persistence,double cylrds,double c=2.4)
//Reference: Wang J, Gao H. The Journal of chemical physics, 2005, 123(8): 084906-084906-13.
{
    double MaxRodlength=std::pow(4.0*cylrds*cylrds*persistence,1.0/3.0)/c;
    bool result=rodlength < MaxRodlength;
    if(!result)
    {
        std::cerr << "rod (in the cylinder) should < "
            << MaxRodlength << std::endl;
        exit(0);
    }
}

void checkStepValid(double beadRadius,double rodlength,double persistence, double T,double viscosity,double step)
//Reference: Wang J, Gao H. The Journal of chemical physics, 2005, 123(8): 084906-084906-13.
{
    double maxstep1 = ( 0.06*M_PI*viscosity*beadRadius*(rodlength*rodlength*rodlength) )/(persistence*M_KB*T);
    double Dtmp = (M_KB*T) /(6*M_PI*viscosity*beadRadius);
    double maxstep2= (0.06*rodlength)*(0.06*rodlength)/ (2*Dtmp);
    bool result=double(step) < std::min(maxstep1, maxstep2);
    if(!result)
    {
        std::cerr << "step(should) < " <<
            std::min(maxstep1,maxstep2) << std::endl;
        exit(0);
    }
}

#endif

#ifndef DIFFUSION_H
#define DIFFUSION_H

#include "type.h"

class Diffusion
{
public:
	inline void setParameter(double beadRadius,double T,double viscosity);
    bool operator() (const VectorXV3d &R);
private:
	inline double getBeadDiffusion();
	inline Matrix3d rotnePragerTensor(const Vector3d &p1,const Vector3d &p2);
public:
	MatrixXd diffusionMatrix;
	double beadDiffusion;
private:	
	double _beadRadius;
	double _T;
	double _viscosity;
};

void Diffusion::setParameter(double beadRadius,double T,double viscosity)
{
	_beadRadius = beadRadius;
	_T = T;
	_viscosity = viscosity;
}

double Diffusion::getBeadDiffusion ()
{
    return M_KB*_T/(6.0*M_PI*_viscosity*_beadRadius);
}

bool Diffusion::operator() (const VectorXV3d &R) 
{
	if(diffusionMatrix.size()!=R.size()*3)
		diffusionMatrix.resize(R.size()*3,R.size()*3);
	beadDiffusion=getBeadDiffusion();
    for (int i=0;i<R.size();++i)
		diffusionMatrix.block<3,3>(i*3,i*3)=beadDiffusion*Matrix3d::Identity();
    for (int i=0;i<R.size();++i)
        for (int j=0;j<i;++j)
			diffusionMatrix.block<3,3>(j*3,i*3) = diffusionMatrix.block<3,3>(i*3,j*3) = rotnePragerTensor(R[i],R[j]);
	return true;
}

//the Rotne-Prager tensor.
//JENS R, PRAGER S. Variation treatment of hydrodynamic interaction in polymers [J].
//The Journal of Chemical Physics, 1969, 50(11): 4831-4836.
Matrix3d Diffusion::rotnePragerTensor(const Vector3d &p1,const Vector3d &p2)
{
    Vector3d Utmp=p1-p2;
    double Dtmp=Utmp.norm();
    Utmp/=Dtmp;
    Matrix3d tmp,result;
    tmp.noalias()=Utmp*Utmp.transpose();
    if (Dtmp>2*_beadRadius){
        result.noalias() =
            (M_KB*_T)/(8.0*M_PI*_viscosity*Dtmp)
            *( (Matrix3d::Identity()+tmp)
            +(2.0*_beadRadius*_beadRadius)/(Dtmp*Dtmp)
            *(1.0/3.0*Matrix3d::Identity()-tmp) );
    }
    else {
        result.noalias() =
            (M_KB*_T)/(6.0*M_PI*_viscosity*_beadRadius)
            *( (1.0-9.0/32.0*Dtmp/_beadRadius)*Matrix3d::Identity()
            +(3.0/32.0*Dtmp/_beadRadius)*tmp );
    }   
	return result;
}

#endif // DIFFUSION_H

#ifndef FANGLESPRING2_H
#define FANGLESPRING2_H

#include <assert.h>
#include "common.h"

inline Vector3d Fanglespring2(double stiffness,double theta0,const Vector3d &u1,const Vector3d &u2, double d2)
//Spring 2: linear angle spring 
//Spring energy: U = 0.5*k*(theta-theta0)^2, where theta is the angle.
//1. Vector3d Fanglespring2 (stiffness, theta0, u1, u2, d2)
//   k=stiffness, the argument theta0 is the equilibrium angle size,
//   if the argument u1 is the direction of side 1,
//   the argument u2 is the direction of side 2
//   and the argument d2 is length of side 2,
//   then the returned value is the force on the far end of side 2.
{
    const double cosTheta = u1.dot(u2);
    const double sinTheta = std::sqrt(1.0 - cosTheta*cosTheta);

    //when sinTheta = 0, return zero vector.
    if(sinTheta < 1e-30) return Vector3d::Zero();

    const double cotTheta = cosTheta/sinTheta;

    const double cosTheta0 = std::cos(theta0);
    const double sinTheta0 = std::sin(theta0);

    const double cosDetlaTheta = cosTheta*cosTheta0+sinTheta*sinTheta0;
    const double absDetlaTheta = std::acos(cosDetlaTheta);

    double factor = 1.0;
    if(absDetlaTheta > 0.05)
    {
        const double absSinDetlaTheta = std::sqrt(1-cosDetlaTheta*cosDetlaTheta);
        factor = absDetlaTheta/absSinDetlaTheta;
    }

    return factor*stiffness/d2*(cosTheta0-sinTheta0*cotTheta)*(u1-cosTheta*u2);
}

inline void Fanglespring2(double stiffness,double theta0,int il,int ic,int ir, LINCS &lincs)
//2. void Fanglespring2 (stiffness, theta0, il, ic, ir, lincs)
//   stiffness, theta0 are defined above.
//   il is the index of the bead at the far end of one side,
//   ir to the far end of another side, and ic to the central bead.
//   lincs is type of LINCS which contain a chain.
{
    assert(ic<lincs.R.size() && il<lincs.R.size() && ir<lincs.R.size() );
    Vector3d Ul=lincs.R[il]-lincs.R[ic];
    double Dl=Ul.norm();
    Ul/=Dl;
    Vector3d Ur=lincs.R[ir]-lincs.R[ic];
    double Dr=Ur.norm();
    Ur/=Dr;

    const Vector3d fl=Fanglespring2(stiffness,theta0,Ur,Ul,Dl);
    lincs.plusF(il,fl);
    const Vector3d fr=Fanglespring2(stiffness,theta0,Ul,Ur,Dr);
    lincs.plusF(ir,fr);
    lincs.minusF(ic,fl+fr);
}

inline Vector3d Fanglespring2_pi(double stiffness,const Vector3d &u1, const Vector3d &u2,double d2)
//1. Vector3d Fanglespring2_pi (stiffness, u1, u2, d2).
//   It's the special case of the function Fanglespring2
//   when the theta0 = pi.
{
    const double cosTheta = u1.dot(u2);
    const double absSinTheta = std::sqrt(1.0-cosTheta*cosTheta);
    const double absTheta = std::acos(cosTheta);

    double factor = 1.0;
    if(absSinTheta > 0.05)
        factor = (M_PI-absTheta)/absSinTheta;

    return factor*stiffness/d2*(cosTheta*u2-u1);
}

inline void Fanglespring2_pi(double stiffness,int il,int ic,int ir,LINCS &lincs)
//2. void Fanglespring2_pi (stiffness, il, ic, ir, lincs).
//   It's the special case of the function Fanglespring2
//   when the theta0 = pi.
{
    assert(ic<lincs.R.size() && il<lincs.R.size() && ir<lincs.R.size() );
    Vector3d Utmpj=lincs.R[il]-lincs.R[ic];
    double Dtmpj=Utmpj.norm();
    Utmpj/=Dtmpj;
    Vector3d Utmpk=lincs.R[ir]-lincs.R[ic];
    double Dtmpk=Utmpk.norm();
    Utmpk/=Dtmpk;
    Vector3d jtmp=Fanglespring2_pi(stiffness,Utmpk,Utmpj,Dtmpj);
    lincs.plusF(il,jtmp);
    Vector3d ktmp=Fanglespring2_pi(stiffness,Utmpj,Utmpk,Dtmpk);
    lincs.plusF(ir,ktmp);
    lincs.minusF(ic,jtmp+ktmp);
}

#endif // FANGLESPRING2_H

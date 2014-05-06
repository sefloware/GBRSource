#ifndef LINCS_H_
#define LINCS_H_

#include <cmath>
#include <time.h>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include "type.h"
#include "geometricConstraint.h"
#include "diffusion.h"
extern "C" {
#include "lapacke.h"
#include "cblas.h"
}

class LINCS
{
public:
    LINCS(double temperature,double viscosity,double timeStep):
        temperature(temperature),timeStep(timeStep),viscosity(viscosity),
        _numGC(0),_normalGen(generator, normal(0.0,1.0))
    {
        generator.seed(static_cast<unsigned>(time(0)));

        fixBool.setConstant(false);
    }

    inline void plusF(unsigned i,Vector3d f) //load the force f on the ith bead.
    { F.segment<3>(i*3) += f; }
    inline void minusF(unsigned i, Vector3d f) //load the force -f on the ith bead.
    { F.segment<3>(i*3) -= f; }

    inline void plusX(unsigned i,Vector3d x) //load the displacement x on the ith bead.
    { X.segment<3>(i*3) += x; }
    inline void minusX(unsigned i,Vector3d x) //load the displacement -x on the ith bead.
    { X.segment<3>(i*3) -= x; }

    inline unsigned endInsertOneGC() // end of the inserting of one geometric constraint.
    { return geometricConstraint.nextRow()-1; }
    inline void setGCdirection(unsigned i, Vector3d u) // set the direction of the geometric constraint on the ith bead.
    {  geometricConstraint.pushGradient(i,u); }
    inline void setGCerror(double error) // set the error of the geometric constraint on the ith bead.
    {  geometricConstraint.pushError(error); }

    inline void setSeed(unsigned seed) // set the random seed.
	{ generator.seed(seed); }
    inline void setBeadRadius(double beadRadius) // set the bead Radius.
    {
        this->beadRadius = beadRadius;
        diffusion.setParameter(beadRadius,temperature,viscosity);
    }

    inline double beadDiffusion() const //return the diffusion of the beads.
    { return diffusion.beadDiffusion; }

    inline Vector3d forceGC(unsigned constraintOrder,unsigned beadId) const
    {
        // return the force of the geometric constraint 'constraintOrder' on the bead 'beadId'.
        const double lamda = -(_lamda[constraintOrder]/( timeStep/(M_KB*temperature)));
        Vector3d force = geometricConstraint.gradient.block(beadId*3,constraintOrder,3,1);
        force *= lamda;
        return force;
    }

    inline void reset()
    {
        //reset the F, X, geometricConstraint to zero and recalculate the diffusion matrix.
        F.setZero();
		X.setZero();
        if(_numGC) geometricConstraint.setZero();
		diffusion(R);
    }

    inline void conservativeResizeGC(unsigned numGC)
    {
        //Resizes geometric constraint number leaving old values untouched.
        _numGC=numGC;
        if(_numGC)
        {
            geometricConstraint.conservativeResize(NoChange,numGC);
            _DBt.conservativeResize(NoChange,numGC);
            _BDBt.conservativeResize(numGC,numGC);
            _lamda.conservativeResize(numGC);
        }
    }

    inline void resizeGC(unsigned numGC,unsigned reserveSizes = 6)
    {
        //Resizes geometric constraint number.
        _numGC=numGC;
        if(_numGC)
        {
            unsigned n3 = R.size() * 3;
            geometricConstraint.resize(n3,numGC);
            geometricConstraint.reserve(reserveSizes);
            _DBt.resize(n3,numGC);
            _BDBt.resize(numGC,numGC);
            _lamda.resize(numGC);
        }
    }

    inline void swapBeadPositions(VectorXV3d &R0)
    {
        //swaps R with R0. Import and export chain's positions.
        if( R0.size() && R.size() != R0.size())
        {
            const unsigned n3=R0.size()*3;
            F.resize(n3);
            X.resize(n3);
            if(_numGC) geometricConstraint.resize(n3,_numGC);
        }
        R.swap(R0);
    }

    inline void fixCoordinateValue(int axisIndice,double value)
    {
        //fix the coordinate component's value.
        fixBool[axisIndice] = true;
        fixValue[axisIndice] = value;
    }

    inline void freeCoordinates()  //free any fixation on the coordinate components.
    { fixBool.setConstant(false); }

	bool operator()()
	{
        //computing the new position with one step of the time based on the LINCS algorithm.
		_N01distGen();

        llt = diffusion.diffusionMatrix;
        LAPACKE_dpotrf(LAPACK_COL_MAJOR,'L',R.size()*3,llt.data(),R.size()*3);
        cblas_dtrmv(CblasColMajor,CblasLower,CblasNoTrans,CblasNonUnit,R.size()*3,llt.data(),R.size()*3,_randomMove.data(),1);
        cblas_daxpy(R.size()*3,std::sqrt(timeStep*2.0),_randomMove.data(),1,X.data(),1);

        const double tdivkT = timeStep/(M_KB*temperature);
        cblas_dgemv(CblasColMajor,CblasNoTrans,R.size()*3,R.size()*3,tdivkT,diffusion.diffusionMatrix.data(),R.size()*3,F.data(),1,1.0,X.data(),1);

        for (unsigned indice=0; indice<3; ++indice)
            if(fixBool[indice])
                for( int i=0; i<X.size(); i+=3)
                    X[i]=fixValue[indice];

        if(_numGC)
        {
            geometricConstraint.error+=(geometricConstraint.gradient.transpose()*X);

            if(_DBt.rows() != R.size()*3)
                _DBt.resize(R.size()*3, _numGC);

            _DBt.noalias()=geometricConstraint.gradient.transpose()*diffusion.diffusionMatrix;
            _BDBt.noalias()=_DBt*geometricConstraint.gradient;

            LAPACKE_dpotrf(LAPACK_COL_MAJOR,'L',_numGC,_BDBt.data(),_numGC);
            LAPACKE_dpotrs(LAPACK_COL_MAJOR,'L',_numGC,1,_BDBt.data(),_numGC,geometricConstraint.error.data(),_numGC);

            _lamda.swap(geometricConstraint.error);
            cblas_dgemv(CblasColMajor,CblasNoTrans,R.size()*3,_numGC,1.0,_DBt.data(),R.size()*3,_lamda.data(),1,-1.0,X.data(),1);
        }

        for(int i=0;i<R.size();++i)
			R[i]+=X.segment<3>(i*3);
		return true;
	}
private:
	inline void _N01distGen()
	{
		if(_randomMove.size()!=R.size()*3)
			_randomMove.resize(R.size()*3);
        for(int i=0;i<_randomMove.size();++i)
			_randomMove[i]=_normalGen();
	}
public:
    VectorXV3d R; // the position collection vector.
    base_generator_type generator; // the random generator.
    double temperature, timeStep, viscosity, beadRadius;
private:
    Vector3d fixValue;
    Vector3b fixBool;
    unsigned _numGC;

    VectorXd F,X;
    Diffusion diffusion;
    MatrixXd llt;

    GeometricConstraint geometricConstraint;
	VectorXd _randomMove;
	MatrixXd _BDBt,_DBt;
	VectorXd _lamda;

	typedef boost::random::normal_distribution<double> normal;
	typedef boost::variate_generator<base_generator_type&, normal> gen_type;
	gen_type _normalGen;
};

#endif // LINCS_H_

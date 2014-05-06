#ifndef LINCS_H_
#define LINCS_H_

#include <cmath>
#include <time.h>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include "type.h"
#include "geometricConstraint.h"
#include "diffusion.h"

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
        X.noalias()+=(llt.compute(diffusion.diffusionMatrix).matrixL()*_randomMove)*std::sqrt(timeStep*2.0);

        const double tdivkT = timeStep/(M_KB*temperature);
        X.noalias() += tdivkT*(diffusion.diffusionMatrix*F);

        for (unsigned indice=0; indice<3; ++indice)
            if(fixBool[indice])
                for( int i=0; i<X.size(); i+=3)
                    X[i]=fixValue[indice];

        if(_numGC)
        {
            geometricConstraint.error+=(geometricConstraint.gradient.transpose()*X);

            if(_DBt.rows() != R.size()*3)
                _DBt.resize(R.size()*3, _numGC);
            _DBt.noalias()=diffusion.diffusionMatrix*geometricConstraint.gradient;
            _BDBt.noalias()=geometricConstraint.gradient.transpose()*_DBt;

            _lamda=ldlt.compute(_BDBt).solve(geometricConstraint.error);

            X.noalias()-=_DBt*_lamda;
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
    GeometricConstraint geometricConstraint;
	VectorXd _randomMove;
	MatrixXd _BDBt,_DBt;
	VectorXd _lamda;
	LLT<MatrixXd> llt;
	LDLT<MatrixXd> ldlt;

	typedef boost::random::normal_distribution<double> normal;
	typedef boost::variate_generator<base_generator_type&, normal> gen_type;
	gen_type _normalGen;
};

#endif // LINCS_H_

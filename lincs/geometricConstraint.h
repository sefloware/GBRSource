#ifndef GEOMETRICCONSTRAINT_H_
#define GEOMETRICCONSTRAINT_H_

#include "type.h"

class GeometricConstraint
{
public:  
    inline void conservativeResize(unsigned rows, unsigned cols)
    {
        error.conservativeResize(cols);
        gradient.conservativeResize(rows,cols);
    }

    inline void resize(unsigned N3,unsigned NC)
    {
        error.resize(NC);
        gradient.resize(N3,NC);
    }

    inline void reserve(unsigned reserveSizes)
    { gradient.reserve(reserveSizes);         }

    inline void pushError(double x)
    { error[rowOffset] = x; }

    inline void pushGradient(unsigned beadId,const Vector3d &x) {
        beadId *= 3;
        gradient.insert(  beadId,rowOffset)=x[0];
        gradient.insert(++beadId,rowOffset)=x[1];
        gradient.insert(++beadId,rowOffset)=x[2];
    }

    inline void setZero() {
        gradient.setZero();
        rowOffset = 0;
    }

    inline unsigned nextRow()
    { return ++rowOffset; }

public:
    SparseMatrix<double> gradient; // the transpose matrix of B.
    VectorXd error; // the error vector.
private:
    unsigned rowOffset;
};

#endif // GEOMETRICCONSTRAINT_H_

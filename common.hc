class LINCS
//The linear constraint solver class.
{
public:    
    void plusF(unsigned i,Vector3d f);
    //Load the force f on the ith bead.

    void minusF(unsigned i, Vector3d f);
    //Load the force -f on the ith bead.

    void plusX(unsigned i,Vector3d x);
    //Load the displacement x on the ith bead.

    void minusX(unsigned i,Vector3d x);
    //Load the displacement -x on the ith bead.

    unsigned endInsertOneGC();
    //End the inserting of one geometric constraint.

    void setGCdirection(unsigned i, Vector3d u);
    //Set the direction of the geometric constraint on the ith bead.

    void setGCerror(double error);
    //Set the error of the geometric constraint on the ith bead.

    void setSeed(unsigned seed);
    //Set the random seed.

    void setBeadRadius(double beadRadius);
    //Set the bead Radius.

    double beadDiffusion();
    //Return the diffusion of the beads.

    Vector3d forceGC(unsigned constraintOrder,unsigned beadId);
    //Return the force of the geometric constraint 'constraintOrder' on the bead 'beadId'.

    void reset();
    //Set the F, X and geometricConstraint object to zero, and recalculate the diffusion matrix.

    void conservativeResizeGC(unsigned numGC);
    //Resizes geometric constraint leaving old values untouched.

    void resizeGC(unsigned numGC,unsigned reserveSizes = 6);
    //Resizes geometric constraint.

    void swapBeadPositions(VectorXV3d &R0);
    //Swaps R's datas with R0's datas.
	//Import and export beads' positions.

    void fixCoordinateValue(int axisIndice,double value);
    //Fix the coordinate component's value.

    void freeCoordinates();
    //Free any fixation on the coordinate components.

    bool operator()();
    //computing the new position by one step time based on the LINCS algorithm.

public:
    VectorXV3d R; //The collection vector of the bead's positions
    base_generator_type generator; //The random generator.
    double temperature; //The temperature of the solvent.
    double timeStep; //The time step.
    double viscosity; //The viscosity of the solvent.
    double beadRadius; //The radius of the beads.
};

//---Eigen---\
class Matrix3d
//The type of a 3*3 matrix.
{
public:
    double &operator()(int i,int j);

    static Matrix3d Identity(); //Return an identity 3*3 matrix.
    Matrix3d &setIdentity(); //identify the matrix.
    static Matrix3d Zero(); //Return an zero 3*3 matrix.
    Matrix3d &setZero(); //set the matrix to zero matrix.
    static Matrix3d Ones(); //Return an 3*3 matrix with elements 1.
    Matrix3d &setOnes(); //set the matrix's element to 1.
    static Matrix3d Constant(double value); //Return an 3*3 matrix with elements value.
    Matrix3d &setConstant(double value); //set the matrix's elements to value.

    Matrix3d transpose(); 
	//Return a transpose matrix of the matrix.
	//Note: m = m.transpose() is wrong! Use the "transposeInPlace()" method instead.
    void transposeInPlace(); 
	//Transpose the matrix.
};

class RowVector3d
//The type of a 3-D row vector
{
public:
    double &operator()(int index); //Return the reference to the index component.
    double &operator[](int index); //Return the reference to the index component.
    double &x(); //Return the reference to x component.
    double &y(); //Return the reference to y component.
    double &z(); //Return the reference to z component.

    double dot(const RowVector3d &x); //Return the dot product between the vector and the row vector x.
    RowVector3d cross(const RowVector3d &x);  //Return the cross product between the vector and the row vector x.

    double norm(); //Return the norm of the vector.
    double squaredNorm(); //Return the squared norm of the vector, e.g. norm()^2.

    RowVector3d normalized();  //Return a normalized vector of the vector.
    void normalize(); //Normalize the vector.

    static RowVector3d UnitX(); //Return a row vector [1,0,0]
    static RowVector3d UnitY(); //Return a row vector [0,1,0]
    static RowVector3d UnitZ(); //Return a row vector [0,0,1]
    static RowVector3d Unit(int i); //Return a row vector with elements of 1 in index i, and 0 in other indexs.

    Vector3d transpose();
	//Return a transpose vector(a column vector) of the vector.
	//Note: m = m.transpose() is wrong! Use the "transposeInPlace()" method instead.
    void transposeInPlace();
	//Transpose the vector.

    static RowVector3d Identity(); //Return a row vector [1,0,0].
    RowVector3d &setIdentity(); //Assign the vector to [1,0,0] and return the reference to the vector.
    static RowVector3d Zero(); //Return a row vector [0,0,0].
    RowVector3d &setZero(); //Assign the vector to [0,0,0] and return the reference to the vector.
    static RowVector3d Ones(); //Return a row vector [1,1,1].
    RowVector3d &setOnes(); //Assign the vector to [1,1,1] and return the reference to the vector.
    static RowVector3d Constant(double value); //Return a row vector [value,value,value].
    RowVector3d &setConstant(double value); //Assign the vector to [value,value,value] and return the reference to the vector.
};

class Vector3d
//The type of a 3-D column vector
{
public:
    double &operator()(int index); //Return the reference to the index component.
    double &operator[](int index); //Return the reference to the index component.
	double &x(); //Return the reference to x component.
    double &y(); //Return the reference to y component.
    double &z(); //Return the reference to z component.

    double dot(const Vector3d &x);  //Return the dot product between the vector and the column vector x.
    Vector3d cross(const Vector3d &x);  //Return the cross product between the vector and the column vector x.
    
    double norm(); //Return the norm of the vector.
    double squaredNorm(); //Return the squared norm of the vector, e.g. norm()^2.

    Vector3d normalized();  //Return a normalized vector of the vector.
    void normalize(); //Normalize the vector.
	
    static Vector3d UnitX(); //Return a column vector [1,0,0]
    static Vector3d UnitY(); //Return a column vector [0,1,0]
    static Vector3d UnitZ(); //Return a column vector [0,0,1]
    static Vector3d Unit(int i); //Return a column vector with elements of 1 in index i, and 0 in other indexs.

    RowVector3d transpose(); 
	//Return a transpose vector(a row vector) of the vector.
	//Note: m = m.transpose() is wrong! Use the "transposeInPlace()" method instead.
    void transposeInPlace(); //Transpose the vector.

    static Vector3d Identity(); //Return a column vector [1,0,0].
    Vector3d &setIdentity(); //Assign the vector to [1,0,0] and return the reference to the vector.
    static Vector3d Zero(); //Return a column vector [0,0,0].
    Vector3d &setZero(); //Assign the vector to [0,0,0] and return the reference to the vector.
    static Vector3d Ones(); //Return a column vector [1,1,1].
    Vector3d &setOnes(); //Assign the vector to [1,1,1] and return the reference to the vector.
    static Vector3d Constant(double value); //Return a column vector [value,value,value].
    Vector3d &setConstant(double value); //Assign the vector to [value,value,value] and return the reference to the vector.
};

class VectorXV3d
//The collect vector of datas of type Vector3d.
{
public:
    Vector3d &operator[](int index); //Return the reference to the index component.
    Vector3d &operator()(int index); //Return the reference to the index component.

    void resize(unsigned size); //Resize the vector.
    void conservativeResize(unsigned size); //Resize the vector leaving old values untouched.
    Vector3d *data() const; //Return a pointer to the first element of the vector.
    VectorXV3d &setConstant(const Vector3d &x); //Set the vector with elements x.
};
//Reference: http://eigen.tuxfamily.org/dox/group__QuickRefPage.html
//---------/

//---C++ Standards---\
namespace std
//A few of the c++ standards.
{
double sin(double arg); //Computes sine of arg.
double asin(double arg); //Computes arc sine of arg.
double cos(double arg); //Computes cosine of arg.
double acos(double arg); //Computes arc cosine of arg.
double tan(double arg); //Computes tangent of arg.
double atan(double arg); //Computes arc tangent of arg.
double sinh(double arg); //Computes hyperbolic sine of arg.
double cosh(double arg); //Computes hyperbolic cosine of arg.
double tanh(double arg); //Computes hyperbolic tangent of arg.
double exp(double arg); //Computes the e (Euler's number, 2.7182818) raised to the given power arg.
double log(double arg); //Computes the natural (base e) logarithm of arg.
double log10(double arg); //Computes the common (base 10) logarithm of arg.
double sqrt(double arg); //Computes square root of arg.
double ceil(double arg);  //Computes nearest integer not less than arg.
double floor(double arg); //Computes nearest integer not greater than arg.
double abs(double arg); //Computes the absolute value of an integer number.
double pow(double base, double exp); //Computes the value of base raised to the power exp.
class ofstream { }; 
ostream cout; 
ostream cerr;
ostream clog;
ostream& endl (ostream& os);

namespace ios
{ enum {app, ate , binary, in, out, trunc};}

}

enum { 
        M_PI =   3.14159265358979323846,
	    M_KB =   1.3806505E-23,
     };

//Reference: http://www.cplusplus.com/
//---------/

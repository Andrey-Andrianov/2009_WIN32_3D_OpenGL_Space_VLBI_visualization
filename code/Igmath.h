#ifndef _IGMATH_H_
#define _IGMATH_H_

#include <math.h>
#include "igtypes.h"
#include "ATemplate.h"


	#define M_E			2.7182818284590452354
	#define M_LOG2E		1.4426950408889634074
	#define M_LOG10E	0.43429448190325182765
	#define M_LN2		0.69314718055994530942
	#define M_LN10		2.30258509299404568402
	#define M_PI		3.14159265358979323846
	#define M_2PI		6.28318530717958647692
	#define M_PI_2		1.57079632679489661923
	#define M_PI_4		0.78539816339744830962
	#define M_1_PI		0.31830988618379067154
	#define M_2_PI		0.63661977236758134308
	#define M_2_SQRTPI	1.12837916709551257390
	#define M_SQRT2		1.41421356237309504880
	#define M_SQRT1_2	0.70710678118654752440

const	double ANearlyZero	=	1E-15; //1E-38;

void CalcHistogram(double *data, int n, TMnMx& mnmx,
				   double *hist, int nhist);

void CalcDispMean(double *data, int n, double& mean, double& disp);

Boolean IsGauss(double* y, int n, const double& mean, const double& disp);

class AGaussNoise
{
protected:
	Boolean	m_bUseLast;
	double	y2;
public:
	AGaussNoise();
	double Get(double mean, double sigma);
};

//double GradInRad(const double& angle);
//double RadInGrad(const double& angle);

double MaxSol3(double a, double b, double c);

int Just2Power(int val);

CString JulDateToDate(double jdate, int& day, int& month, int& year);
CString JulDateToDate1(double jdate, int& day, int& month, int& year);
double Date_JDate(int day, int month, int year, double time = 0);
double MapFunc(double gst, double raapp, double decapp, 
				 double antX, double antY, double antZ);
void RADECEpo_RADECApp(double raepo, double decepo, double jd, int epoch,
			   double& raapp, double& decapp);
double T_GST(double jdt);

int A_SLEQ_Gauss(ADMATRIX& A, ADMATRIX& B, double eps,
	ADMATRIX &X, int& n, double& det);
int A_SLEQ_Gauss_new(ADMATRIX& A, ADMATRIX& B, double eps,
	ADMATRIX &X, int& n, double& emax);

class ALinearMinSQR
{
protected:
	ADMATRIX	m_A, m_B, m_X;
public:
	ALinearMinSQR();
	ALinearMinSQR(int size);

	Boolean Create(int size);
	int Add(double *a, double b);
	int Calculate(double eps);
	ADMATRIX& X()		{ return (m_X); };
};


//double log2(double x);
//double AINT(const double& val);
//double FRAC(const double& val);


Boolean Cross2Line(const CPoint& pt1_1, const CPoint& pt2_1,
				   const CPoint& pt1_2, const CPoint& pt2_2,
				   CPoint& ptCross);

class ABAproximation
{
protected:
	typedef float afloat;

public:
	virtual ~ABAproximation()				{};
	virtual double fitX(afloat x)			{ return (0); };
	virtual int  Calculate(afloat* xData, afloat* yData, 
		int ptStart, int ptEnd)		{ return (0); };
};

class AAproximation : public ABAproximation
{
protected:
	typedef float afloat;

protected:
	int numTerms;
	afloat *pmBasis, *pmSolution;
	afloat multiplier, constant;
	
	AFMATRIX *coefficients;
	afloat*	constants;

public:
	AAproximation(int power);
	virtual ~AAproximation();
	afloat* GetSlolution()			{ return (pmSolution); };
	virtual double fitX(afloat x);
	virtual int  Calculate(afloat* xData, afloat* yData, 
					int ptStart, int ptEnd);
protected:
	virtual void GetBasis(afloat* basis, int n, afloat x) {};
	virtual void Transform(afloat* xData, afloat* yData, 
							int ptStart, int ptEnd)
			{ multiplier = 1;  constant = 0; }
	virtual void InverseTransform(afloat& yFit) {};
	virtual void CreateBasisFunctions(int numPoints) {};
	virtual void TransformSolution(afloat* oldSolution,	
									afloat* newSolution) {};
private:
	int InitializeAndFormBasisVectors(int ptStart, int ptEnd,
			afloat* xData, afloat* yData, afloat* solution);
	void Initial(int dimen, afloat* solution, int& error);
	void EROswitch(afloat* row1, afloat* row2);
	void EROmultAdd(afloat multiplier, int dimen, 
				afloat* referenceRow, afloat* changingRow)
		{ for (int i = 0; i < dimen; ++i)
			changingRow[i] += multiplier * referenceRow[i]; };
	void BackwardsSub(int dimen, afloat* solution);
	void ComputeNormalEquations(int ptStart, int ptEnd, afloat* xData, afloat* yData);
	void Pivot(int dimen, int referenceRow, int& error);
	void UpperTriangular(int dimen, int& error);
	void Partial_Pivoting(int dimen, afloat* solution, int& error);
	int CreateAndSolveEquations(int ptStart, int ptEnd,	afloat* solution, 
								afloat* xData, afloat* yData);
	void TransformSolutionAndFindResiduals(int ptStart, int ptEnd, 
								afloat* yData, afloat* solution)
			{  TransformSolution(solution, solution); } // My comment ?
};


class AAproxPoly : public AAproximation
{
public:
	AAproxPoly(int power) : AAproximation(power) {};
protected:
	virtual void GetBasis(afloat* basis, int n, afloat x);
	virtual void Transform(afloat* xData, afloat* yData, 
							int ptStart, int ptEnd);
	virtual void TransformSolution(afloat* oldSolution,	
							afloat* newSolution) {};
};
class AAproxUser : public AAproxPoly
{
public:
	AAproxUser(int power) : AAproxPoly(power) {};
protected:
	virtual void GetBasis(afloat* basis, int n, afloat x);
	virtual void Transform(afloat* xData, afloat* yData, 
							int ptStart, int ptEnd);
};

class AAproxFourier : public AAproximation
{
	double m_mean, m_mul;
public:
	AAproxFourier(int power) : AAproximation(power) {};

	virtual double fitX(afloat x);
	virtual int  Calculate(afloat* xData, afloat* yData, 
					int ptStart, int ptEnd);
protected:
	virtual void Transform(afloat* xData, afloat* yData, 
							int ptStart, int ptEnd);
private:
	void ReDFTLong(afloat* x, afloat* y, int ptStart, int ptEnd,
		 	 	   afloat* re, afloat* im, int k, int direct);
};

class AAproxPower : public ABAproximation
{
protected:
	double m_a, m_alpha;	
public:
	AAproxPower()							{ m_a = 1.0; m_alpha = 1.0; };
	double& Alpha()							{ return (m_alpha); };
	double& A()								{ return (m_a); };
	virtual double fitX(afloat x)			{ return (m_a * pow(double(x), m_alpha)); };
	virtual int  Calculate(afloat* xData, afloat* yData, 
				int ptStart, int ptEnd);
};

class AEstimateBeam
{
	double	m_x2, m_y2, m_xy,
			m_mx, m_my, m_sz,
			m_angle, m_sFWHM, m_lFWHM;
public:
	AEstimateBeam();
	void Add(const double& x, const double& y, const double& z = 1);
	void Calculate();
	double GetAngle()	{ return (m_angle); };
	double GetLFWHM()	{ return (m_lFWHM); };
	double GetSFWHM()	{ return (m_sFWHM); };
};

/*
 PLeastFourier = ^TLeastFourier;
 TLeastFourier = Object(TAprBase)
   procedure  GetBasis(var Basis: array of Float; N: Integer; X: Float); virtual;
 end;
*/

class AGaussFit
{
	double	m_sum_x2,
			m_sum_x4,
			m_sum_lny,
			m_sum_x2_lny,
			m_a,
			m_b;
	int		m_count;
public:
	AGaussFit();
	void add(double x, double y);
	void calculate();
	virtual double FitX(const double& x);
};

class AGaussTeylorFit
{
public:
	double	m_sum_x2,
			m_sum_x4,
			m_sum_x6,
			m_sum_x8,
			m_sum_x10,
			m_sum_x12,
			m_sum_y,
			m_sum_x2_y,
			m_sum_x4_y,
			m_sum_x6_y,
			m_a,
			m_b,
			m_c,
			m_d;
	int		m_count;
	TMnMx	m_mnmx;

public:
	AGaussTeylorFit();

	void add(const double& x, const double& y);
	void calculate();
	double FitX(double x);
};


class AB
{
	ADMATRIX	m_a,  
				m_b, 
				m_x;
public:
	AB(ADMATRIX& a, ADMATRIX& b);

	int Calculate(const double& epsilon);
	ADMATRIX& Result() { return (m_x); };

private:
	void Normalize();

};


class AMGauss
{
//public:
	ADMATRIX	m_a, 
				m_b, 
				m_x;
public:
	AMGauss(ADMATRIX& a, ADMATRIX& b);
	void Normalize();
	ADMATRIX* Calculate();
};

class AMGauss1
{
//public:
	ADMATRIX	m_a;
	ADARRAY		m_b, 
				m_x;
public:
	AMGauss1(ADMATRIX& a, ADARRAY& b);
	void Normalize();
	ADARRAY* Calculate();
};

class AMeanIntervals
{
	TMnMx m_mnmx;
	double	*m_yMean;
	double	*m_xMean;
	int m_count;
public:
	double *m_y, *m_y2, *m_n;

	AMeanIntervals(int n, TMnMx& mnmx);
	virtual ~AMeanIntervals();

	int count()	{ return (m_count); };
	double GetDX()	{return (m_mnmx.delta() / (count() - 1)); };
	void Add(const double& x, const double& y);
	void Squeeze();
	double* GetYMas()	{ return (m_yMean); };
	double* GetXMas()	{ return (m_xMean); };
};

class A3DBFunc
{
public:
	virtual ~A3DBFunc()				{};
	virtual double ValidXWorld()	{ return (0); };
	virtual double ValidYWorld()	{ return (0); };
	virtual double GetValue(double x, double y)		{ return (1); };
	virtual double GetFValue(double x, double y)	{ return (1); };
	virtual AFMATRIX* toMatrix(double dx, double dy);
	virtual AFMATRIX* toMatrix1(double dx, double dy);
};

class A3DPillBoxFunc : public A3DBFunc
{
	double m_wx, m_wy;
public:
	A3DPillBoxFunc(double wx, double wy);
	virtual double GetFValue(double x, double y);
};

class A3DGaussFunc : public A3DBFunc
{
	double	m_sigmaX,
			m_sigmaY,
			m_angle;
public:
	A3DGaussFunc(double FWHMx, double FWHMy, double angle);
	virtual double ValidXWorld();
	virtual double ValidYWorld();
	virtual double GetValue(double x, double y);
	virtual double GetFValue(double x, double y);
};

class AInterpLagrangeXCoeff
{
	double *m_xCoeff;
	int		m_count;
public:
	AInterpLagrangeXCoeff()		{ m_xCoeff = NULL; m_count = 0; };
	AInterpLagrangeXCoeff(double xc, double *x, int n);
	~AInterpLagrangeXCoeff()	{ if (m_xCoeff) delete [] m_xCoeff; };

	void SetNewXMas(double xc, double *x, int n);
	int count()				{ return (m_count); };
	double* GetXCoeff()		{ return (m_xCoeff); };
protected:
	void CalculateCoeff(double xc, double *x, int n);
};

class AInterpolate
{
public:
	virtual double InterpolateX(double x)		{ return (0); };
};

class ACubicSlineFree : public AInterpolate
{
protected:
	double	*m_b;
	double	*m_d;
	double	*m_c;

	double	*m_x;
	double	*m_y;
	int		m_npts;
public:
	ACubicSlineFree();
	ACubicSlineFree(double *x, double *y, int n);
	virtual ~ACubicSlineFree();

	void SetParameters(double *x, double *y, int n);
	virtual double InterpolateX(double x);
protected:
	void CalcCoeff();
};

class APhJustify
{
protected:
	double	m_a;	// a - prediduschee ne korrectirovanoe
	double	m_phi;	// phi - prediduschee korrectirovanoe
	Boolean	m_bFirstUse;

public:
	APhJustify()			{ m_bFirstUse = true; };
	double Correct(double phase);
};


//
//------------------------------------------------------------
/*
inline double GradInRad(const double& angle)
{
	return (angle * M_PI / 180.0);
};

inline double RadInGrad(const double& angle)
{
	return (angle * 180.0 / M_PI);
};
*/
inline int Just2Power(int val)
{
	return (1 << int(log(double(val + (val - 1))) / M_LN2));
}
/*
inline double log2(double x)
{
	return (log(x) / M_LN2);
}


inline double AINT(const double& val)
{
	double ret;
	modf(val, &ret);
	return (ret);
};
inline double FRAC(const double& val)
{
	double ret;
	return (modf(val, &ret));
};
*/
#endif

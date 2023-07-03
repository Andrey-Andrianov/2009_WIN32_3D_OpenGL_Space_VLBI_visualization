#include "igmath.h"
#include "ATemplate.cpp"
#include "astring.h"
#include "AMathQuick.h"
#include <stdlib.h>

const  double RAND_MAX_MUL = 2.0 / double(RAND_MAX);

double MaxSol3(double a, double b, double c)
{
	double d = - a/3.0;
	double p = b + a*d;
	double q = c + d*(b-2*d*d);
	double d0 = -0.5*q;
	double r0 = p/3.0;
	double r = d0*d0 + r0*r0*r0;
	double b0 = sqrt(fabs(r));
	double x;
	if (r < 0) 
	{
		double phi0 = atan2(b0, d0);
		x = 2*sqrt(-r0)*cos(phi0/3.0/*-M_2PI/3.0*/);
	}
	else
		x = SQRT3(d0+b0) + SQRT3(d0-b0);
	x += d;
	return (x);
};

AAproximation::AAproximation(int power)
{
	numTerms = power;
	pmSolution	= new afloat[numTerms];
	pmBasis		= new afloat[numTerms];
};

AAproximation::~AAproximation()
{
	delete [] pmSolution;
	delete [] pmBasis;
};
/*
constructor TAprBase.Load;
begin
 S.Read(NumTerms,SizeOf(Integer));
 GetMem(pmSolution,NumTerms*SizeOf(Float));
 GetMem(pmBasis,NumTerms*SizeOf(Float));
 S.Read(pmSolution^,NumTerms*SizeOf(Float));
 S.Read(Multiplier,SizeOf(Float));
 S.Read(Constant,SizeOf(Float));
end;
destructor TAprBase.Done;
begin
 FreeMem(pmSolution,NumTerms*SizeOf(Float));
 FreeMem(pmBasis,NumTerms*SizeOf(Float));
 inherited Done;
end;
procedure TAprBase.Store;
begin
 S.Write(NumTerms,SizeOf(Integer));
 S.Write(pmSolution^,NumTerms*SizeOf(Float));
 S.Write(Multiplier,SizeOf(Float));
 S.Write(Constant,SizeOf(Float));
end;
*/

double AAproximation::fitX(afloat x)
{
	x = multiplier * x + constant;
	GetBasis(pmBasis, numTerms, x);
	afloat y = 0;
	for (int i = 0; i < numTerms; ++i)
		y += pmSolution[i] * pmBasis[i];
	return (y);
};

//-------------------------------------------------------
int AAproximation::InitializeAndFormBasisVectors(int ptStart, int ptEnd,
			afloat* xData, afloat* yData, afloat* solution)
{
	int np = ptEnd - ptStart + 1;
	memset(solution, 0, numTerms * sizeof(afloat));
	if (np < 2) return (1);   // Less than 2 data points  
	if (numTerms < 1) return(2);   // Less than 1 coefficient in the fit  }
	if (numTerms > np) return (3);
    Transform(xData, yData, ptStart, ptEnd);
//    if Error = 0 then
	CreateBasisFunctions(np);
	return (0);
};
void AAproximation::Initial(int dimen, afloat* solution, int& error)
{
	error = 0;
	if (dimen < 1) error = 1;
	else
		if (dimen == 1)
		{
			if (fabs((*coefficients)(0, 0)) < ANearlyZero) error = 2;
			else
				solution[0] = constants[0] / (*coefficients)(0, 0);
		};
};
void AAproximation::BackwardsSub(int dimen, afloat* solution)
{
	int term = dimen - 1;
	while (term >= 0)
	{
		double sum = 0;
		for (int row = term + 1; row < dimen; ++row)
			sum += (*coefficients)(term, row) * solution[row];
		solution[term] = float((constants[term] - sum) / (*coefficients)(term, term));
		term--;
	};
};
void AAproximation::ComputeNormalEquations(int ptStart, int ptEnd, afloat* xData, afloat* yData)
{
	memset(constants, 0, numTerms * sizeof(afloat));
	for (int index = ptStart; index <= ptEnd; index++)
	{
		afloat x = multiplier * xData[index] + constant;
		GetBasis(pmBasis, numTerms, x);
		for (int column = 0; column < numTerms; ++column)
		{
			constants[column] += yData[index] * pmBasis[column];
			for (int row = column; row < numTerms; ++row)
			{
				afloat a = pmBasis[row] * pmBasis[column];
				(*coefficients)(row, column) += a;
				if (row != column)
					(*coefficients)(column, row) += a;
			};
		};
	};
};
void AAproximation::EROswitch(afloat* row1, afloat* row2)
{
	for (int i = 0; i < numTerms; ++i) SWAP(row1[i], row2[i]);
};
void AAproximation::Pivot(int dimen, int referenceRow, int& error)
{
// First, find the row with the largest element  
	int pivotRow = referenceRow;
	for (int row = referenceRow + 1; row < dimen; ++row)
		if (fabs((*coefficients)(row, referenceRow)) >
			fabs((*coefficients)(pivotRow, referenceRow)) ) pivotRow = row;
	if (pivotRow != referenceRow)
// Second, switch these two rows  
	{
		EROswitch(&(*coefficients)(pivotRow, 0), &(*coefficients)(referenceRow, 0)); //??
		afloat dummy = constants[pivotRow];
		constants[pivotRow] = constants[referenceRow];
		constants[referenceRow] = dummy;
	}
	else
		if (fabs((*coefficients)(referenceRow, referenceRow)) < ANearlyZero)
			error = 2;     // No solution  
};
void AAproximation::UpperTriangular(int dimen, int& error)
{
	int referenceRow = 0;
	while ((error == 0) && (referenceRow < dimen - 1))
	{
	// Find row with largest element in this column  
	// and switch this row with the ReferenceRow     
		Pivot(dimen, referenceRow, error);
		if (error == 0)
			for (int row = referenceRow + 1; row < dimen; ++row)
			{
			// Make the ReferenceRow element of these rows zero 
				if (fabs((*coefficients)(row, referenceRow)) > ANearlyZero)
				{
					afloat multiplier = -(*coefficients)(row, referenceRow) /
										(*coefficients)(referenceRow, referenceRow);
					EROmultAdd(multiplier, dimen, &(*coefficients)(referenceRow, 0), 
								&(*coefficients)(row, 0));
					constants[row] += multiplier * constants[referenceRow];
				};
			};
			referenceRow++;
	}; // while 
	if (fabs((*coefficients)(dimen - 1, dimen - 1)) < ANearlyZero)
		error = 2;		// No solution
};
void AAproximation::Partial_Pivoting(int dimen, afloat* solution, int& error)
{
	Initial(dimen, solution, error);
	if (dimen > 1)
	{
		UpperTriangular(dimen, error);
		if (error == 0)
			BackwardsSub(dimen, solution);
	};
};
int AAproximation::CreateAndSolveEquations(int ptStart, int ptEnd,
							afloat* solution, afloat* xData, afloat* yData)
{
	int error = 0;
	coefficients = new AFMATRIX(numTerms, numTerms);
	memset(coefficients->lpArray(), 0, 
		coefficients->GetiSize() * coefficients->GetjSize() * sizeof(afloat));
	constants = new afloat[numTerms];
	memset(constants, 0, numTerms * sizeof(afloat));

	ComputeNormalEquations(ptStart, ptEnd, xData, yData);
	Partial_Pivoting(numTerms, solution, error);
	if (error == 2) // Returned from Partial_Pivoting
		error = 4; // No solution
	delete [] constants;
	delete coefficients;
	return (error);
};
int AAproximation::Calculate(afloat* xData, afloat* yData, 
							int ptStart, int ptEnd)
{
	int error = InitializeAndFormBasisVectors(ptStart, ptEnd, xData,
						yData, pmSolution);
	if (error == 0)
		error = CreateAndSolveEquations(ptStart, ptEnd, pmSolution, xData, yData);
	if (error == 0)
		TransformSolutionAndFindResiduals(ptStart, ptEnd, yData, pmSolution);
	return (error);
};
//
//------------------------------------------------------
void AAproxPoly::GetBasis(afloat* basis, int n, afloat x)
{
	basis[0] = 1;
	basis[1] = x;
	for (int i = 2; i < n; ++i)
		basis[i] = 2 * x * basis[i - 1] - basis[i - 2];
};

void AAproxPoly::Transform(afloat* xData, afloat* yData, 
							int ptStart, int ptEnd)
{
	afloat xDataMin = xData[0];
	afloat xDataMax = xData[0];
	for (int row = ptStart; row <= ptEnd; ++row)
	{
		afloat x = xData[row];
		if (xDataMin > x) xDataMin = x;
		if (xDataMax < x) xDataMax = x;
	};
	multiplier = float(2.0 / (xDataMax - xDataMin));
	constant = -xDataMin * multiplier - 1;
//	multiplier = 1.0 / (xDataMax - xDataMin);
//	constant = -xDataMin * multiplier - 0.5;
};
//
//-------------------------------------------------------
void AAproxUser::GetBasis(afloat* basis, int n, afloat x)
{
	basis[0] = 1;
	for (int i = 1; i < n; ++i)
		basis[i] = x * basis[i - 1];
};
//???
void AAproxUser::Transform(afloat* xData, afloat* yData, 
							int ptStart, int ptEnd)
{
//	multiplier = 1;  constant = 0;
//	return;
	afloat xDataMin = xData[0];
	afloat xDataMax = xData[0];
	for (int row = ptStart; row <= ptEnd; ++row)
	{
		afloat x = xData[row];
		if (xDataMin > x) xDataMin = x;
		if (xDataMax < x) xDataMax = x;
	};
//	multiplier = 1.0 / (xDataMax - xDataMin);
//	constant = -xDataMin * multiplier;
	multiplier = float(2.0 / (xDataMax - xDataMin));
	constant = -xDataMin * multiplier - 1;
//	multiplier = float(1.0 / xDataMax);
//	constant = 0;
};
//???
//
//----------------------------------------------
double AAproxFourier::fitX(afloat x)
{
	double y = 0;
	for (int i = 0; i < numTerms; ++i)
	{
		double xx = constant + multiplier * x;
		xx = 1 - xx;
		double arg = M_2PI * i * xx;
		y += (pmBasis[i] * cos(arg) + pmSolution[i] * sin(arg));
	};
	return (y + m_mean);
};

void AAproxFourier::Transform(afloat* xData, afloat* yData, 
							int ptStart, int ptEnd)
{
	afloat xDataMin = xData[0];
	afloat xDataMax = xData[0];
	double sum = 0;
	for (int i = ptStart; i <= ptEnd; ++i)
	{
		afloat x = xData[i];
		if (xDataMin > x) xDataMin = x;
		if (xDataMax < x) xDataMax = x;
		sum += yData[i];
	};
//	multiplier = float(1.0 / (xDataMax - xDataMin));
	multiplier = float((1.0 - 1.0 / double(ptEnd - ptStart + 1)) / (xDataMax - xDataMin));
	constant = -xDataMin * multiplier;
//??
	m_mean = sum / (ptEnd - ptStart + 1);
//	m_mean = 0;
};

void AAproxFourier::ReDFTLong(afloat* x, afloat* y, int ptStart, int ptEnd,
							  afloat* re, afloat* im, int k, int direct)
{
	double cnst = 2.0 / double(ptEnd - ptStart + 1);
	for (int j = 0; j < k; ++j)
	{
		double sc = 0;
		double ss = 0;
		double temp = M_2PI * j;
		m_mul = 0;
		for (int i = ptStart; i <= ptEnd; ++i)
		{
			double angle = (constant + multiplier * x[i]) * temp;
			double yy = y[i] - m_mean;
			sc += yy * cos(angle);
			ss += yy * sin(angle);
		};
		if (direct < 0)
		{
			sc *=  cnst;
			ss *= -cnst;
		};
		re[j] = float(sc);
		im[j] = float(ss);
	};
};

int AAproxFourier::Calculate(afloat* xData, afloat* yData, 
					int ptStart, int ptEnd)
{
	Transform(xData, yData, ptStart, ptEnd);
	ReDFTLong(xData, yData, ptStart, ptEnd, pmBasis, pmSolution, numTerms, -1);
	return (0);
};
//
///////////////////////////////////////////////
int AAproxPower::Calculate(afloat* xData, afloat* yData, 
				int ptStart, int ptEnd)
{
	double sx=0, sy=0, sxy=0, sx2=0;
	int n=0;
	for (int i = ptStart; i <= ptEnd; ++i)
	{
		double y = yData[i];
		double x = xData[i];
/*		if (y < 0) return (-1);
		if (x < 0) return (-2);
		if (y == 0) y = 0;
		else y = log(y);
		if (x == 0) x = 0;
		else x = log(x);*/
		if (y == 0) y = -100.0;
		else y = log(y);
		if (x == 0) x = -100.0;
		else x = log(x);

		sy += y;
		sx += x;
		sxy += y*x;
		sx2 += x*x;
		n++;
	};
	m_alpha = (sy*sx/n - sxy) / (SQR(sx)/n - sx2);
	//m_alpha = (sxy) / (sx2);
	m_a = (sy - m_alpha*sx) / n;
	m_a = exp(m_a);
	return (0);
};

//
//----------------------------------------
//
AGaussFit::AGaussFit()
{
	m_sum_x2 = 0;
	m_sum_x4 = 0;
	m_sum_lny = 0;
	m_sum_x2_lny = 0;
	m_count = 0;
};

void AGaussFit::add(double x, double y)
{
	m_sum_x2_lny += log(y) * x * x;
	m_sum_lny += log(y);
	m_sum_x2 += x * x;
	m_sum_x4 += x * x * x * x;
	m_count++;
};

void AGaussFit::calculate()
{
	m_b = (m_sum_x2_lny - m_sum_lny * m_sum_x2 / m_count) /
		  (m_sum_x2 * m_sum_x2 / m_count - m_sum_x4);
	m_a = exp((m_b * m_sum_x2 + m_sum_lny) / m_count);
//	AfxMessageBox(toStr(m_b)+" "+toStr(m_a));
};

double AGaussFit::FitX(const double& x)
{
	double a = - SQR(x) * m_b;
	if (a < -1e4) a = -1e4;
	a = m_a * exp(a);
//	AfxMessageBox(toStr(x)+" "+toStr(a));
	return (a);
};
//
//-------------------------------------------------------
/*
procedure TLeastFourier.GetBasis;
{
	basis[0] = 1;
	basis[1] = cos(x);
	basis[2] = sin(X);
	for (int i = 3; i <n; ++i)
	{
		if odd(Column) then Basis[Column]:=Basis[1] * Basis[Column-2]
                                     -Basis[2] * Basis[Column-1]
		else
			Basis[Column]:=Basis[2] * Basis[Column-3]
                   +Basis[1] * Basis[Column-2];
	};
};
*/

void CalcHistogram(double *data, int n, TMnMx& mnmx,
					double *hist, int nhist)
{
	memset(hist, 0, nhist * sizeof(double));
	double c = 1.0 / n;
	double c1 = nhist / mnmx.delta();
	for (int i = 0; i < n; ++i)
	{
		int j = int((data[i] - mnmx.min) * c1);
		if (j >= nhist) j = nhist - 1;
		hist[j] = hist[j] + c;
	};
};

void CalcDispMean(double *data, int n, double& mean, double& disp)
{
	double sum1 = 0;
	double sum2 = 0;
	for (int i = 0; i < n; ++i)
	{
		sum1 += data[i];
		sum2 += SQR(data[i]);
	};
	mean = sum1 / n;
	disp = (sum2 - SQR(sum1) / n) / (n - 1);
};

double mom3(double* x, int n, const double& mean)
{
	double sum = 0.0;
	for (int i = 0; i < n; ++i)
	{
		double d = x[i] - mean;
		sum += SQR(d) * d;
	};
	return (sum / n);
};

double mom4(double* x, int n, const double& mean)
{
	double sum = 0.0;
	for (int i = 0; i < n; ++i)
		sum += SQR(x[i] - mean);
	return (sum / n);
};

Boolean SqewKur(const double& m2, const double& m3, 
				const double& m4, int n)
{
	double k2 = m2 / (1.0 - 1.0 / n);
	double k3 = m3 / ((1.0 - 1.0 / n) * (1.0 - 2.0 / n));
	double k4 = (m4 / ((1.0 - 2.0 / (n + 1.0)) * (1.0 - 2.0 / n) 
		* (1.0 - 3.0 / n))) - 3.0 * SQR(m2) / ((1.0 - 2.0 / n)
		* (1.0 - 3.0 / n));
	double g1 = k3 / (k2 * sqrt(k2));
	double g2 = k4 / SQR(k2);
	double dg1 = 6.0 * n * (n - 1.0) / ((n - 2.0) * (n + 1.0) * (n + 3.0));
	double dg2 = 24.0 * n * SQR(n - 1.0) / ((n - 3) * (n - 2.0) * (n + 3.0) * (n + 5.0));
	return ( (fabs(g1) <= sqrt(dg1)) && (fabs(g2) <= sqrt(dg2)) );
};

Boolean IsGauss(double* y, int n, const double& mean, const double& disp)
{
	double m2 = disp * (1.0 - 1.0 / n);
	double m3 = mom3(y, n, mean);
	double m4 = mom4(y, n, mean);
	return (SqewKur(m2, m3, m4, n));
};
//
//----------------------------------------------
//
AB::AB(ADMATRIX& a, ADMATRIX& b)
{
	m_a = a;
	m_b = b;
	Normalize();
	m_x = m_b;
	ASSERT ((a.GetiSize() == b.GetiSize()) && (b.GetjSize() == 1));
};

void AB::Normalize()
{
	for (int i = 0;  i < m_a.GetiSize(); ++i)
	{
		m_b(i, 0) /=  m_a(i, i);
		double a = m_a(i, i);
		for (int j = 0; j < m_a.GetjSize(); ++j)
		{
			m_a(i, j) = - m_a(i, j) / a;
		};
		m_a(i, i) = 0;
	};
}; 

int AB::Calculate(const double& epsilon)
{
	double max;
	int count = 0;
	do 
	{
		ADMATRIX tmp;
		tmp = m_x * m_a + m_b;
		max = -1e38; 
		for (int i = 0; i < tmp.GetiSize(); ++i)
		{
			max = max(max, fabs(m_x(i, 0) - tmp(i, 0)));
		};
		m_x = tmp;
		count++;
	} while (max > epsilon);

	return (count);
};
//
//----------------------------------
//
AGaussTeylorFit::AGaussTeylorFit() : m_mnmx(1e38, -1e38)
{
	m_sum_x2	= 0;
	m_sum_x4	= 0;
	m_sum_x6	= 0;
	m_sum_x8	= 0;
	m_sum_y		= 0;
	m_sum_x2_y	= 0;
	m_sum_x4_y	= 0;
	m_count		= 0;
};

void AGaussTeylorFit::add(const double& x, const double& y)
{
	m_mnmx.check(x);
	double xx = x * x;
	double x4 = xx * xx;
	m_sum_x2	+= xx;
	m_sum_x2_y	+= xx * y;
	m_sum_x4	+= x4;
	m_sum_x4_y	+= x4 * y;
	m_sum_x6	+= x4 * xx;
	m_sum_x8	+= x4 * x4;
	m_sum_y		+= y;
	m_count++;
};

void AGaussTeylorFit::calculate()
{
	ADMATRIX a(3, 3), b(3);
	a(0, 0) = m_count; a(0, 1) = m_sum_x2; a(0, 2) = m_sum_x4; 
	a(1, 0) = m_sum_x2; a(1, 1) = m_sum_x4; a(1, 2) = m_sum_x6;
	a(2, 0) = m_sum_x4; a(2, 1) = m_sum_x6; a(2, 2) = m_sum_x8; 
	b(0, 0) = m_sum_y; b(1, 0) = m_sum_x2_y; b(2, 0) = m_sum_x4_y;

	AMGauss ab(a, b);
	ADMATRIX* x = ab.Calculate();
	m_a = (*x)(0, 0);
	m_b = (*x)(1, 0);
	m_c = (*x)(2, 0);
};

double AGaussTeylorFit::FitX(double x)
{
	double xx = x * x;
	return (m_a + m_b * xx + m_c * xx * xx);
};

//
//------------------------
//
AMGauss::AMGauss(ADMATRIX& a, ADMATRIX& b)
{
	m_a = a;
	m_b = b;
	m_x = b;
};

void AMGauss::Normalize()
{
	for (int i = 0; i < m_a.GetiSize() - 1; ++i)
	{
		for (int j = i + 1; j < m_a.GetjSize(); ++j)
		{
			double a = -m_a(j, i) / m_a(i, i);
			for (int k = i; k < m_a.GetjSize(); ++k)
				m_a(j, k) += a * m_a(i, k);
			m_b(j, 0) += a * m_b(i, 0);
		};
	};
};

ADMATRIX* AMGauss::Calculate()
{
	Normalize();
	for (int i = m_a.GetiSize() - 1; i >= 0; --i)
	{
		double sum = 0;
		for (int j = i + 1; j < m_a.GetiSize(); ++j)
		{
			sum += m_a(i, j) * m_x(j, 0);
		};
		m_x(i, 0) = (m_b(i, 0) - sum) / m_a(i, i);
	};
	return (&m_x);
};

//////////////////////////
AMGauss1::AMGauss1(ADMATRIX& a, ADARRAY& b)
{
	m_a = a;
	m_b = b;
	m_x = b;
};

void AMGauss1::Normalize()
{
	for (int i = 0; i < m_a.GetiSize() - 1; ++i)
	{
		for (int j = i + 1; j < m_a.GetjSize(); ++j)
		{
			double a = -m_a(j, i) / m_a(i, i);
			for (int k = i; k < m_a.GetjSize(); ++k)
				m_a(j, k) += a * m_a(i, k);
			m_b(j) += a * m_b(i);
		};
	};
};

ADARRAY* AMGauss1::Calculate()
{
	Normalize();
	for (int i = m_a.GetiSize() - 1; i >= 0; --i)
	{
		double sum = 0;
		for (int j = i + 1; j < m_a.GetiSize(); ++j)
		{
			sum += m_a(i, j) * m_x(j);
		};
		m_x(i) = (m_b(i) - sum) / m_a(i, i);
	};
	return (&m_x);
};


CString JulDateToDate(double jdate, int& day, int& month, int& year)
{
	int  moff[13] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365};
	int moff1[13] = {0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366};
	//year = 2000 - int((2451545 - jdate) / 365.2422);
	//if (jdate < 2451545) year--;
	year = 2000 - int((2451544.5 - jdate) / 365.2422);
	if (jdate < 2451544.5) year--;
	int n = year - 1;
	int ic = n / 100;
	int nday = 1;
	if ( ( (((year - 1) % 4) == 0) && ((year % 100) != 0) ) ||
		((((year - 1) % 4) == 0) && ((year % 100) == 0) && 
		(((year / 100) % 4) == 0)) ) nday++;
	if((year > 2001) && (((year - 1) % 4) != 0 )) nday++;
	double jd_jan1 = 1721424 - ic + 365.0 * n + (n + ic) / 4 + nday - 0.5;
	double diff = jdate - jd_jan1;

	if (!((((year % 4) == 0) && ((year  % 100) != 0)) ||
		(((year % 4) == 0) && ((year % 100) == 0) && 
		(((year / 100) % 4) == 0))))
	{
	    for (int i = 0; i < 12; ++i)
		{
			if ((diff >= moff[i]) && (diff < moff[i + 1]))
			{
				month = i + 1;
				day = int(diff) - moff[i] + 1;
			};
		};
	}
	else
    for (int i = 0; i < 12; ++i)
	{
		if ((diff >= moff1[i]) && (diff < moff1[i + 1]))
		{
			month = i + 1;
			day = int(diff) - moff1[i] + 1;
		};
	};
	CTime time(year, month, day, 0, 0, 0, 0);

	return (time.Format("%d-%m-%Y"));
};

CString JulDateToDate1(double jdate, int& day, int& month, int& year)
{
	int  moff[13] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365};
	int moff1[13] = {0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366};
	//year = 2000 - int((2451545 - jdate) / 365.2422);
	//if (jdate < 2451545) year--;
	year = 2000 - int((2451544.5 - jdate) / 365.2422);
	if (jdate < 2451544.5) year--;
	int n = year - 1;
	int ic = n / 100;
	int nday = 1;
	if ( ( (((year - 1) % 4) == 0) && ((year % 100) != 0) ) ||
		((((year - 1) % 4) == 0) && ((year % 100) == 0) && 
		(((year / 100) % 4) == 0)) ) nday++;
	if((year > 2001) && (((year - 1) % 4) != 0 )) nday++;
	double jd_jan1 = 1721424 - ic + 365.0 * n + (n + ic) / 4 + nday - 0.5;
	double diff = jdate - jd_jan1;
	if (!((((year % 4) == 0) && ((year  % 100) != 0)) ||
		(((year % 4) == 0) && ((year % 100) == 0) && 
		(((year / 100) % 4) == 0))))
	    for (int i = 0; i < 12; ++i)
		{
			if ((diff >= moff[i]) && (diff < moff[i + 1]))
			{
				month = i + 1;
				day = int(diff) - moff[i] + 1;
			};
		}
		else
	    for (int i = 0; i < 12; ++i)
		{
			if ((diff >= moff1[i]) && (diff < moff1[i + 1]))
			{
				month = i + 1;
				day = int(diff) - moff1[i] + 1;
			};
		};

	CTime time(year, month, day, 0, 0, 0, 0);

	return (time.Format("%d/%m/%y"));
};


double MapFunc(double gst, double raapp, double decapp, 
				 double antX, double antY, double antZ)
{
	double s1 = gst - atan2(antY, antX);
	double t1 = s1 - raapp*M_PI/180.0;
	double phi = atan(antZ/sqrt(SQR(antX)+SQR(antY)));
	double cose = sin(phi)*sin(decapp*M_PI/180.0)+cos(phi)*cos(decapp*M_PI/180.0)*cos(t1);
	return (1.0/cose);
};

void RADECEpo_RADECApp(double raepo, double decepo, double jd, int epoch,
			   double& raapp, double& decapp) 
{
	double t = (jd - 2415020.5) / 36524.22;
	double m = 46.08506 + 0.027945*t + 0.00012*t*t;
	double n = 20.04685 - 0.0008533*t - 0.00037*t*t;

	double z_e = 0;
	if (epoch == 1950) z_e = 2433282.5;
	else
	if (epoch == 1900) z_e = 2415020.5;
	else 
		z_e = 2451544.5; // epoch == 2000

	raapp = raepo + ((jd-z_e)/365.2422) 
		* (m + n*sin(raepo)*tan(decepo))/3600.0;
	if (raapp == M_PI/2) 
		raapp = raepo + (jd-z_e)*m/(365.2422*3600.0);

	decapp = decepo + ((jd-z_e)/365.2422) * (n*cos(raepo))/3600.0;
};

double T_GST(double jdt)
{
	double d = FRAC(jdt);
	double jd;
	if (d > 0.5) jd = AINT(jdt) + 0.5;
	else jd = AINT(jdt) - 0.5;
	double jt = jdt - jd;

	double dt = (jd - 2451545.0) / 36525; //??
	double st = (24110.55 + (8640184.81 + 0.093 * dt) * dt)
				* M_PI / 43200.0;
	st = fmod(st, M_2PI);
	if (st < 0) st += M_2PI;
	return (st + jt * M_2PI * (366.25 / 365.25));
};
//
////////////////////////////////////////////////////
AEstimateBeam::AEstimateBeam()
{
	m_x2 = 0;
	m_y2 = 0;
	m_xy = 0;
	m_mx = 0;
	m_my = 0;
	m_sz = 0;
};
void AEstimateBeam::Add(const double& x, const double& y, const double& z)
{
	m_sz += z;
	m_x2 += x * x * z;
	m_y2 += y * y * z;
	m_xy += x * y * z;
	m_my += x * z;
	m_mx += y * z;
};

void AEstimateBeam::Calculate()
{

/*
	m_mx = (m_mx / m_sz);
	m_my = (m_my / m_sz);
	double m02 = (m_x2 / m_sz);
	double m20 = (m_y2 / m_sz);
	double m11 = (m_xy / m_sz);
	double sx2 = (m02 - SQR(m_mx));
	double sy2 = (m20 - SQR(m_my));
	double pxy = (m11 - m_mx * m_my);
	double beta = 0.5 * atan2(2 * pxy, sx2 - sy2);
	double cbt2 = SQR(cos(beta));
	double sbt2 = SQR(sin(beta));
	double sa2 = ((sx2 * cbt2 - sy2 * sbt2) / (cbt2 - sbt2));
	double sb2 = (sx2 * sbt2 - sy2 * cbt2) / (sbt2 - cbt2);
	double sx = 1.0 / (2 * sqrt(2.0) * M_PI * sqrt(sa2));
	double sy = 1.0 / (2 * sqrt(2.0) * M_PI * sqrt(sb2));
	m_sFWHM = 2.0 * sqrt(2 * M_LN2) * sx * (180.0 * 3600.0 / M_PI);
	m_lFWHM = 2.0 * sqrt(2 * M_LN2) * sy * (180.0 * 3600.0 / M_PI);
	m_angle = (0.5 * M_PI - (beta + 0.5 * M_PI)) * 180 / M_PI;
	if (m_angle > 90.0) m_angle -= 180.0;
	if (m_angle < -90.0) m_angle = 180.0 + m_angle;	
//	AfxMessageBox(toStr(sx)+" "+toStr(sy2));

*/
// ---------- My Code -----------------

	double m01 = m_mx / m_sz;
	double m10 = m_my / m_sz;
	double m02 = (m_x2 / m_sz);
	double m20 = (m_y2 / m_sz);
	double m11 = (m_xy / m_sz);
	double sigma_u2 = m02 - SQR(m01);
	double sigma_v2 = m20 - SQR(m10);
	double cov_uv = (m11 - m01 * m10);
	double rho_uv = cov_uv/(sqrt(sigma_u2*sigma_v2));
	double beta =-0.5 * atan2(2 * cov_uv , sigma_u2-sigma_v2);
	m_angle = beta * 180 / M_PI;
	double a = 1.0;
	double b = -(sigma_u2+sigma_v2)/(2*(1-rho_uv*rho_uv)*sigma_u2*sigma_v2);
	double c = 1/(4*(1-rho_uv*rho_uv)*sigma_u2*sigma_v2);
	double d = sqrt(b*b - 4*a*c);
		
//	AfxMessageBox("a, b, c, d "+toStr(a)+"  "+toStr(b)+"  "+toStr(c)+"  "+toStr(d));

	double lambda_1 = 0.5*(-b+d);
	double lambda_2 = 0.5*(-b-d);

	
//	AfxMessageBox("rho, l_1, l_2  "+toStr(rho_uv)+"  "+toStr(lambda_1)+"  "+toStr(lambda_2));

	double sigma_nu = 1/sqrt(2*lambda_1);
	double sigma_nv = 1/sqrt(2*lambda_2);
	
//	AfxMessageBox("s_nu, s_nv, beta "+toStr(sigma_nu)+" "+toStr(sigma_nv)+" "+toStr(m_angle));
	
	m_lFWHM = (180.0 * 3600.0 / M_PI) * 0.37478125 / sigma_nu;
	m_sFWHM = (180.0 * 3600.0 / M_PI) * 0.37478125 / sigma_nv;

//	AfxMessageBox("m_sFWHM, m_lFWHM "+toStr(m_sFWHM)+" "+toStr(m_lFWHM)+" "+toStr(m_angle));

	
	if (m_angle > 90.0) m_angle -= 180.0;
	if (m_angle < -90.0) m_angle = 180.0 + m_angle;


// ---------- End My Code -----------------

};

//
/////////////////////////////////////////////////
AMeanIntervals::AMeanIntervals(int n, TMnMx& mnmx)
{
	m_count = n;
	m_mnmx = mnmx;
	m_yMean = new double[n];
	memset(m_yMean, 0, n * sizeof(double));
	m_xMean = new double[n];
	memset(m_xMean, 0, n * sizeof(double));
	m_y = new double[n];
	memset(m_y, 0, n * sizeof(double));
	m_y2 = new double[n];
	memset(m_y2, 0, n * sizeof(double));
	m_n = new double[n];
	memset(m_n, 0, n * sizeof(double));
};

AMeanIntervals::~AMeanIntervals()
{
	delete [] m_yMean;
	delete [] m_xMean;
	delete [] m_y;
	delete [] m_y2;
	delete [] m_n;
};

void AMeanIntervals::Add(const double& x, const double& y)
{
	int i = int(0.5 + (x - m_mnmx.min) * (count() - 1) / m_mnmx.delta());
	m_yMean[i] += y;
	m_xMean[i] ++;
	m_y[i] += y;
	m_y2[i] += SQR(y);
	m_n[i] ++;
};

void AMeanIntervals::Squeeze()
{
	int ndrop = 0;
	for (int i = 0; i < count(); ++i)
	{
		if (m_xMean[i] == 0)
		{
			ndrop++;
		}
		else
		{
			m_yMean[i - ndrop] = m_yMean[i] / m_xMean[i];
			m_xMean[i - ndrop] = m_mnmx.min + i * m_mnmx.delta() / (count() - 1);
			m_y[i - ndrop] = m_y[i];
			m_y2[i - ndrop] = m_y2[i];
			m_n[i - ndrop] = m_n[i];

		};
	};
	m_count -= ndrop;
};
//
////////////////////////////////////////////////
AFMATRIX* A3DBFunc::toMatrix(double dx, double dy)
{	
	int nj = 2 * int(0.5 + ValidXWorld() / dx);
	nj = (nj & ~1) + 1;
	int ni = 2 * int(0.5 + ValidYWorld() / dy);
	ni = (ni & ~1) + 1;
	if (ni > 513) ni = 513;
	if (nj > 513) nj = 513;
	AFMATRIX *matr = new AFMATRIX(ni, nj);
	for (int i = 0; i < ni; ++i)
	{
		for (int j = 0; j < nj; ++j)
		{
			(*matr)(i, j) = float(GetValue((j - nj / 2) * dx,
								(i - ni / 2) * dy));
		};
	};
	return (matr);
};

AFMATRIX* A3DBFunc::toMatrix1(double dx, double dy)
{	
	int nj = 1 * int(0.5 + ValidXWorld() / dx);
	nj = (nj & ~1) + 1;
	int ni = 1 * int(0.5 + ValidYWorld() / dy);
	ni = (ni & ~1) + 1;
	if (ni > 513) ni = 513;
	if (nj > 513) nj = 513;
	AFMATRIX *matr = new AFMATRIX(ni, nj);
	for (int i = 0; i < ni; ++i)
	{
		for (int j = 0; j < nj; ++j)
		{
			(*matr)(i, j) = float(GetValue((j - nj / 2) * dx,
								(i - ni / 2) * dy));
		};
	};
	return (matr);
};


A3DPillBoxFunc::A3DPillBoxFunc(double wx, double wy)
{
	m_wx = wx;
	m_wy = wy;
};

double A3DPillBoxFunc::GetFValue(double x, double y)
{
	if ((x == 0) || (y == 0)) return (1);
	return ((sin(x) / x) * sin(y) / y);
};

A3DGaussFunc::A3DGaussFunc(double FWHMx, double FWHMy, double angle)
{
	m_sigmaX = 0.5 * FWHMx / sqrt(2.0 * M_LN2);
	m_sigmaY = 0.5 * FWHMy / sqrt(2.0 * M_LN2);
	m_angle = -angle;
};

double A3DGaussFunc::ValidXWorld()
{
	double x = 4 * m_sigmaX;
	double x1 = cos(m_angle) * x;
	double y = 4 * m_sigmaY;
	double y1 = -sin(m_angle) * y;
	if (fabs(x1) > fabs(y1)) x = x1;
	else x = y1;
	return (fabs(x));
};

double A3DGaussFunc::ValidYWorld()
{
	double x = 4 * m_sigmaX;
	double x1 = sin(m_angle) * x;
	double y = 4 * m_sigmaY;
	double y1 = cos(m_angle) * y;
	if (fabs(x1) > fabs(y1)) x = x1;
	else x = y1;
	return (fabs(x));
};

double A3DGaussFunc::GetValue(double x, double y)
{
	double x1 = cos(m_angle) * x - sin(m_angle) * y;
	double y1 = sin(m_angle) * x + cos(m_angle) * y;
	double a = -0.5 * (SQR(x1 / m_sigmaX) + SQR(y1 / m_sigmaY));
	if (a < -6e2) a = -6e2;
	return (exp(a));
};

double A3DGaussFunc::GetFValue(double x, double y)
{
	double a = -0.25 * (SQR(x) * SQR(m_sigmaX)
					  + SQR(y) * SQR(m_sigmaY));
	if (a < -6e2) a = -6e2;
	return (exp(a));
};
//
//------------------------------------------------
AInterpLagrangeXCoeff::AInterpLagrangeXCoeff(double xc, double *x, int n)
{
	m_count = n; 
	m_xCoeff = new double[n]; 
	CalculateCoeff(xc, x, n);
};
void AInterpLagrangeXCoeff::SetNewXMas(double xc, double *x, int n)
{
	if (n != count())
	{
		if (m_xCoeff) delete [] m_xCoeff;
		m_xCoeff = new double[n];
	};
	m_count = n;
	CalculateCoeff(xc, x, n);
};


void AInterpLagrangeXCoeff::CalculateCoeff(double xc, double *x, int n)
{
	for (int i = 0; i < n; ++i)
	{
		double mul1 = 1.0;
		double mul2 = 1.0;
		for (int j = 0; j < n; ++j)
		{
			if (i == j) continue;
			mul1 *= (xc - x[j]);
			mul2 *= (x[i] - x[j]);
		}
		m_xCoeff[i] = mul1 / mul2;
	};
};




ACubicSlineFree::ACubicSlineFree()
{
	m_npts = 0;
	m_b = NULL;
	m_d = NULL;
	m_c = NULL;
};
ACubicSlineFree::ACubicSlineFree(double *x, double *y, int n)
{
	m_npts = n;
	m_b = new double[n];
	m_d = new double[n];
	m_c = new double[n];
	m_x = x;
	m_y = y;

	CalcCoeff();
};
ACubicSlineFree::~ACubicSlineFree()
{
	if (m_b) delete [] m_b;
	if (m_d) delete [] m_d;
	if (m_c) delete [] m_c;
};
void ACubicSlineFree::SetParameters(double *x, double *y, int n)
{
	if (m_npts != n)
	{
		if (m_b) delete [] m_b;
		if (m_d) delete [] m_d;
		if (m_c) delete [] m_c;
	};
	m_npts = n;
	m_b = new double[n];
	m_d = new double[n];
	m_c = new double[n];
	m_x = x;
	m_y = y;

	CalcCoeff();
};

void ACubicSlineFree::CalcCoeff()
{
	double x;
	double dl = 0;
	m_d[0] = m_b[0] = 0;

	int i;
    for (i = 1; i < m_npts - 1; ++i)
	{
		double x = m_x[i+1]-m_x[i];
		double xx = m_x[i]-m_x[i-1];
		double alpha = 3 * ((m_y[i+1] * xx)
                  - (m_y[i] * (m_x[i+1]-m_x[i-1]))
                  + (m_y[i-1] * x))
                  / (xx * x);

		dl = 2 * (m_x[i+1]-m_x[i-1]) - xx*m_d[i-1];
		m_d[i] = x / dl;
		m_b[i] = (alpha - xx*m_b[i-1]) / dl;
	};
    
    m_c[m_npts-1] = 0;
    for (i = m_npts-2; i >= 0; i--)
      m_c[i] = m_b[i] - m_d[i] * m_c[i+1];

    for (i = m_npts-2; i >= 0; i--)
    {
		x = m_x[i+1] - m_x[i];
		m_b[i] = (m_y[i+1] - m_y[i])/x
				- x*(m_c[i+1] + 2*m_c[i]) / 3;
		m_d[i] = (m_c[i+1] - m_c[i]) / (3*x);
	};
};

double ACubicSlineFree::InterpolateX(double x)
{
	int location = 0;
//	for (int i = 0; i <  m_npts-2; ++i)
	for (int i = 0; i <  m_npts; ++i)
	{
//		location = i;
		if (x <= m_x[i]) break;
		location = i;
	};
//	cout << location<<endl;
	x -= m_x[location];
//	cout<<x<<" "<<m_d[location]<<" "<<
//		m_c[location]<<" "<<
//		m_b[location]<<" "<<" "<<m_y[location]<<endl;
    return (((m_d[location]*x + m_c[location])*x +
			m_b[location])*x + m_y[location]);
}

/*
function CubicSplineFree;
var
 pmC,pmB,pmD: PMFloat0;


var
 Index: Integer;
 Dx: Double;
 Error: Integer;
begin { procedure CubicSplineFree }
  Error := 0;
  if NumPoints < 2 then Error := 3
  else
  begin
   for Index := 0 to NumPoints - 2 do
   begin
    Dx:=XData[Index+1] - XData[Index];
    if Dx < 0 then  Error := 2
    else
     if ABS(Dx) < TNNearlyZero then Error := 1;
   end;
  end;
  if Error = 0 then
  begin
    GetMem(pmC,SizeOf(Float)*NumPoints);
    GetMem(pmB,SizeOf(Float)*NumPoints);
    GetMem(pmD,SizeOf(Float)*NumPoints);
    CalculateCoefficients;
    Interpolate;
    FreeMem(pmC,SizeOf(Float)*NumPoints);
    FreeMem(pmB,SizeOf(Float)*NumPoints);
    FreeMem(pmD,SizeOf(Float)*NumPoints);
  end;
  CubicSplineFree:=Error;
end; { procedure CubicSplineFree }
*/


Boolean Cross2Line(const CPoint& pt1_1, const CPoint& pt2_1,
				   const CPoint& pt1_2, const CPoint& pt2_2,
				   CPoint& ptCross)
{
	int dx1 = pt2_1.x - pt1_1.x;
	int dy1 = pt2_1.y - pt1_1.y;
	int a1 = dy1;
	int b1 = -dx1;
	int c1 = dx1*pt1_1.y - dy1*pt1_1.x;

	int dx2 = pt2_2.x - pt1_2.x;
	int dy2 = pt2_2.y - pt1_2.y;
	int a2 = dy2;
	int b2 = -dx2;
	int c2 = dx2*pt1_2.y - dy2*pt1_2.x;
	//cout<<"a1 = "<<a1<<", b1 = "<<b1<<", a1 = "<<a2<<", b2 = "<<b2<<endl;
	if (a1*b2 == a2*b1) return (false);

	double del = a1*b2 - a2*b1;
	double x = (b1*c2 - b2*c1) / del;
	double y = (c1*a2 - c2*a1) / del;

	ptCross = CPoint(x, y);
	return (true);
};
//--------------------------------------------
double APhJustify::Correct(double phase)
{
	if (m_bFirstUse)
	{
		m_a = phase;
		m_phi = m_a;
		m_bFirstUse = false;
		return (phase);
	};
	double a1 = phase;
	double d = a1-m_a;
	m_a = a1;
	double d1 = d-M_2PI;
	double d2 = d+M_2PI;
	if (fabs(d1) < fabs(d)) d = d1;
	if (fabs(d2) < fabs(d)) d = d2;
	m_phi += d;
	return (m_phi);
};
// correctirovka fazi konetz

//
// ¬ходные параметры
// a - матрица дл€ которой ищетс€ X
// b - вектор (матрица) правой части
// eps - критерий выбора np.
// ¬ыходные параметры
// X - матрица решений 
// n - ранг матрицы ј
// det - детерминант обратимой части матрицы ј
// Ќа выходе вектор (матрица) ¬ - есть искомый вектор ’
//
// возвращаемое значение:
//		признак 0 - норма (решение единственное, n = N2)
//			    1 - решений много
//				2 - система несовместна
//

int A_SLEQ_Gauss(ADMATRIX& A, ADMATRIX& B, double eps,
	ADMATRIX &X, int& n, double& det)
{
	X.ZeroFill();

	int N1 = A.GetiSize();
	int N2 = A.GetjSize();
	int N3 = B.GetjSize();

	AIARRAY J(N2);

	for (int k = 0; k < N2; ++k) J(k) = k;
	n = -1;
	int L = 0;
	det = 1.0;
	double e = eps*1e8;
	double r;
	int mm, km;//??
	goto L2;

L1:		
	SWAP(J(L), J(mm));
	for (int k = L; k < N2; ++k)
	{
		SWAP(A(L,k), A(km, k));
	};
	for (int k = 0; k < N3; ++k)
	{
		SWAP(B(L,k),B(km, k));
	};
	for (int k = 0; k < N1; ++k)
	{
		SWAP(A(k, L),A(k, mm));
	};
	n = L;
	r = 1.0 + 1.0/(N2-L);
	eps *= r * sqrt(r);
	e = eps*1e8;
	L += 1;  
	r = 1.0/A(n,n);
	for (int k = L; k < N2; ++k) A(n,k) *= r;
	for (int k = 0; k < N3; ++k) B(n,k) *= r;
	for (int k = 0; k < N1; ++k)
	{
		if (k != n) 
		{
			r = A(k, n);
			for (int m = L; m < N2; ++m)
				A(k, m) -= r*A(n,m);
			for (int m = 0; m < N3; ++m)
				B(k, m) -= r*B(n,m);
		}; 
	};  
//пока не используетс€	det *= A(n,n);
//	if (mm != n) det = -det;
//	if (km != n) det = -det;
 L2:  
	r = 0.0;  
	km = L;  
	mm = L;
	for (int k = L; k < N1; ++k)
	{
		for (int m = L; m < N2; ++m)
		{
			double s = fabs(A(k, m));
			if (s > r) 
			{
				r = s; 
				km = k;  
				mm = m;
			};
		};
	}; 
	if (r > eps) goto L1;
	n += 1;
	for (int m = 0; m < n; ++m)
	{
		mm = J(m);
		for (int k = 0; k < N3; ++k)
		{
			X(mm, k) = B(m, k);
		};
	};

	for (int k = L; k < N1; ++k)
	{
		for (int m = 0; m < N3; ++m)
		{
			if (fabs(B(k, m)) > e) return (2);
		};
	};
	if (n != N2) return (1);
	return (0);
};


//
// ¬ходные параметры
// a - матрица дл€ которой ищетс€ X
// b - вектор (матрица) правой части
// eps - критерий выбора np.
// ¬ыходные параметры
// X - матрица решений 
// n - ранг матрицы ј
// emax - точность решени€
// Ќа выходе вектор (матрица) ¬ - есть искомый вектор ’
//
// возвращаемое значение:
//		признак 0 - норма (решение единственное, n = N2)
//			    1 - решений много
//				2 - система несовместна
//

int A_SLEQ_Gauss_new(ADMATRIX& A, ADMATRIX& B, double eps,
	ADMATRIX &X, int& n, double& emax)
{
	ADMATRIX A2(A), B2(B);

	X.ZeroFill();

	int N1 = A.GetiSize();
	int N2 = A.GetjSize();
	int N3 = B.GetjSize();

	AIARRAY J(N2);

	for (int k = 0; k < N2; ++k) J(k) = k;

	int q = 0;
	n = 0;
	for (int i = 0; i < min(N1,N2); ++i)
	{
		double r = 0.0;  
		int km = i;
		int mm = i;
		for (int k = i; k < N1; ++k)
		{
			for (int m = i; m < N2; ++m)
			{
				double s = fabs(A(k, m));
				if (s > r) 
				{
					r = s; 
					km = k;  
					mm = m;
				};
			}; // m
		};  // k
		if (r < eps)
		{
			double s = r * emax * (N2 - i);
			for (int k = i; k < N1; ++k)
			{
				for (int m = 0; m < N3; ++m)
				{
					if (fabs(B(k,m)) > s) q = 2;
				};// m
			}; // k
			break;
		}; // if
		if (i != km)
		{
			for (int k = i; k < N2; ++k)
				SWAP(A(i,k), A(km, k));
			for (int k = 0; k < N3; ++k)
				SWAP(B(i,k),B(km, k));
		};
		if (i != mm)
		{
			SWAP(J(i), J(mm));
			for (int k = 0; k < N1; ++k)
				SWAP(A(k, i),A(k, mm));
		};
		n = i+1;
		r = 1.0/A(i,i);
		for (int k = n; k < N2; ++k) A(i,k) *= r;
		for (int k = 0; k < N3; ++k) B(i,k) *= r;
		for (int k = 0; k < N1; ++k)
		{
			if (k == i) continue;
			r = A(k, i);
			for (int m = n; m < N2; ++m)
				A(k, m) -= r*A(i,m);
			for (int m = 0; m < N3; ++m)
				B(k, m) -= r*B(i,m);
		};  //k
	}; // i
	if ((n < N2) && (q != 2)) q = 1;
	
	for (int k = 0; k < N3; ++k)
	{
		for (int m = 0; m < n; ++m)
		{
			double r = 0;
			for (int j = n; j < N2; ++j)
				r += A(m,j)*X(J(j),k);
			X(J(m), k) = B(m, k) - r;
		}; // m
	}; // k
	emax = 0;
	for (int i = 0; i < N1; ++i)
	{
		for (int j = 0; j < N3; ++j)
		{
			double r = B2(i,j); 
			for (int k = 0; k < N2; ++k)
				r -= A2(i,k)*X(k,j);
			emax = max(emax, fabs(r));
		}; // j
	}; // i
	return (q);
};



//////////////////////////////
/* Include requisits */
#include <cstdlib>
#include <ctime>

/* Generate a new random seed from system time - do this once in your constructor */
/*
double GaussNoise(TMnMx range)
{
const int q = 15;
const double c1 = (1 << q) - 1;
const double c2 = ((int)(c1 / 3)) + 1;
const double c3 = 1.0 / c1;

	double random = range.min + range.delta()
		*double(rand()) / (RAND_MAX + 1);
	double noise = (2.0 * ((random * c2) + (random * c2) 
		+ (random * c2)) - 3.0 * (c2 - 1.0)) * c3;
	return (noise);
}
*/



AGaussNoise::AGaussNoise()
{
	m_bUseLast = false;
	y2 = 0;
};

double	AGaussNoise::Get(double mean, double sigma)
{
	double y1, v, x1, x2;
	if (m_bUseLast)
	{
		y1 = y2;
		m_bUseLast = false;
	}
	else
	{
		do 
		{
			x1 = rand() * RAND_MAX_MUL - 1.0;
			x2 = rand() * RAND_MAX_MUL - 1.0;
			v = x1*x1 + x2*x2;
		} while (v >= 1.0);
		v = sqrt(-2.0 * log10(v) /  v);
		y1 = x1 * v;
		y2 = x2 * v;
		m_bUseLast = true;
	};
	return (mean + y1 * sigma);
};

ALinearMinSQR::ALinearMinSQR() : m_A(), m_B(), m_X()
{
};
ALinearMinSQR::ALinearMinSQR(int size) : m_A(), m_B(), m_X()
{
	Create(size);
};

Boolean ALinearMinSQR::Create(int size)
{
	m_A.GrowBy(size, size);
	m_B.GrowBy(size, 1);
	m_X.GrowBy(size, 1);
	m_A.ZeroFill();
	m_B.ZeroFill();
	m_X.ZeroFill();
	return (true);
};

int ALinearMinSQR::Add(double *a, double b)
{
	for (int i = 0; i < m_A.GetiSize(); ++i)
	{
		m_B(i) += b*a[i];
		for (int j = 0; j < m_A.GetjSize(); ++j)
		{
			m_A(i,j) += a[i]*a[j];
		};
	};
	return (0);
};

fstream fs("d:\\ALinearMinSQR.txt", ios::out);

int ALinearMinSQR::Calculate(double eps)
{
	int rang;
	double emax = eps*1e-2;
	int ret = A_SLEQ_Gauss_new(m_A, m_B, eps, m_X, rang, emax);
if (ret != 0)
{
	fs<<"---------------------------"<<endl;
	for (int i = 0; i < m_A.GetiSize(); ++i)
	{
		for (int j = 0; j < m_A.GetjSize(); ++j)
		{
			fs<<"A("<<i<<","<<j<<")="<<m_A(i,j)<<"  ";
		};
		fs<<"B("<<i<<")="<<m_B(i)<<endl;
	};
};
	return (ret);
//	return (A_SLEQ_Gauss(rang, det));
};




double Date_JDate(int day, int month, int year, double time)
//; Var SEC,DJ,IDJ,FDJ: Float);
{
	int MoN[12] = 
	{0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};

	year -= 1900; //??????
	int jd = int(2433282+365.0*(year-50.0)+(year-49.0)/4.0
		+MoN[month-1]+day);
	if (((year % 4) == 0) && (month > 2)) jd++;
	double fjd = time/8.64e4 - 0.5;
	double a = jd + fjd;
	//??
	if (fjd != 0)
	{
		jd--;
		fjd += 1;
	};
	return (a);
};

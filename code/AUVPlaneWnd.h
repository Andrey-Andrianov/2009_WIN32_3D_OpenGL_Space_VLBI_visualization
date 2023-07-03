#pragma once
#include "ASLWnds.h"

struct AUV
{
	double	u;
	double	v;
	double	t;
};

class AUVPlaneGraph : public AGraph
{
protected:
	AUV	*m_mas;
	int	m_count;
	double	m_umax;
	double m_vmax;

public:
	AUVPlaneGraph(ABAxis* aParent, Boolean asItem = false) 
		: AGraph(aParent, asItem) { m_mas = NULL; m_umax = 1; m_vmax = 1;};
	void SetParameters(AUV* mas, int n, double umax, double vmax);
					  

	virtual Boolean hasHeaderText(Boolean postView) { return (!postView); };
	virtual CString getHeaderText(Boolean postView)	{ return ("RadioAstron-Earth"); };
	virtual Boolean isColors() { return (false); };

	virtual void getMinMax(TMinMax& mnmx);
	virtual void draw(CDC& dc, CRect& r, TMinMax& mnmx);
};


// AUVPlaneWnd

class AUVPlaneWnd : public AAxisWnd
{
	DECLARE_DYNAMIC(AUVPlaneWnd)

public:
	AUVPlaneWnd();
	virtual ~AUVPlaneWnd();

	void SetParameters(AUV *mas, int n, double umax, double vmax);
protected:
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnDestroy();
};



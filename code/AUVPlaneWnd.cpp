// AUVPlaneWnd.cpp : implementation file
//

#include "stdafx.h"
//#include "test_uvplane.h"
#include "AUVPlaneWnd.h"


void AUVPlaneGraph::getMinMax(TMinMax& mnmx)
{
	mnmx.x = TMnMx(-m_umax, m_umax);
	mnmx.y = TMnMx(-m_vmax, m_vmax);
};

void AUVPlaneGraph::draw(CDC& dc, CRect& r, TMinMax& mnmx)
{
	if (m_mas == NULL) return;
	int tmp = GetPartCrossHeight();
	SetPartCrossHeight(16);
	for (int i = 0; i < m_count; ++i)
	{
		int x = parent->RPXtoGPX(m_mas[i].u);
		int y = parent->RPYtoGPY(m_mas[i].v);
		parent->pointCross(x, y, parent->workColor(7));
	};
	SetPartCrossHeight(tmp);
};

void AUVPlaneGraph::SetParameters(AUV* mas, int n, double umax, double vmax)
{
	m_mas = mas;
	m_count = n;
	m_umax = umax;
	m_vmax = vmax;
};

// AUVPlaneWnd

IMPLEMENT_DYNAMIC(AUVPlaneWnd, AAxisWnd)

AUVPlaneWnd::AUVPlaneWnd() : AAxisWnd()
{
	SetFlag(WNF_NOFLICKER|WNF_NOWAIT);
}

AUVPlaneWnd::~AUVPlaneWnd()
{
}


void AUVPlaneWnd::SetParameters(AUV *mas, int n, double umax, double vmax)
{
	A2_3DAxis *pa = new A2_3DAxis(-1, 1, NULL, AXIS_DRAWNET
			| AXIS_USEDIB | AXIS_SCALE/* | AXIS_WANTRESTORE*/);
	pa->setNameOX("U [%sWaveLengths]");
	pa->setNameOY("V [%sWaveLengths]");
	AUVPlaneGraph *pi = new AUVPlaneGraph(pa);
	pi->SetParameters(mas, n, umax, vmax);
	SetAxis(pa);
	if(IsWindow(m_hWnd)) InvalidateRect(NULL);
};

BEGIN_MESSAGE_MAP(AUVPlaneWnd, AAxisWnd)
	ON_WM_DESTROY()
END_MESSAGE_MAP()



// AUVPlaneWnd message handlers



void AUVPlaneWnd::OnDestroy()
{
	AAxisWnd::OnDestroy();
    PostQuitMessage(0);
	// TODO: Add your message handler code here
}

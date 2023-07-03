// VisualizationDlg.h : header file
//

#pragma once
#include "afxcmn.h"


// CVisualizationDlg dialog
class CVisualizationDlg : public CDialog
{
// Construction
public:
	CVisualizationDlg(CWnd* pParent = NULL);	// standard constructor

// Dialog Data
	enum { IDD = IDD_VISUALIZATION_DIALOG };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support


// Implementation
protected:
	HICON m_hIcon;

	// Generated message map functions
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnDestroy();
	afx_msg void OnBnClickedOk();
	double m_dec;
	double m_RA;
	double m_UVtime;
	BOOL m_Lines;
	afx_msg void OnDeltaposSpin1(NMHDR *pNMHDR, LRESULT *pResult);
	bool m_Time_scale;
	CSpinButtonCtrl m2_Time_Scale;
	double m_TimeScale;
	afx_msg void OnBnClickedButton1();
	double m_Time;
	double m_TIME;
	CTime start_time;
	COleDateTime start_time2;
	COleDateTime start_time1;
	COleDateTime UV_time2;
	afx_msg void OnBnClickedButton2();
	afx_msg void OnDtnDatetimechangeDatetimepicker4(NMHDR *pNMHDR, LRESULT *pResult);
	int m_UV_time_days;
};

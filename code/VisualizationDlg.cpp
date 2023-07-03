// VisualizationDlg.cpp : implementation file
//

#include "stdafx.h"
#include "Visualization.h"
#include "VisualizationDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif
#include <windows.h>
#include <stdio.h>

#include <gl/gl.h>

#include "constants.h"
#include "3d_view.h"
#include "math.h"

#include <io.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <share.h>
#include "Igmath.h"

typedef  GLuint (* _stdcall MYPROC)(GLenum); 
PROC glaaa;
GLuint (__stdcall *glCreateShader)(GLenum);
#define GL_VERTEX_SHADER                     0x8B31


int nCmdShow;
int TimerID;

int Draw_Timer_Counter = 0;

BOOL CreateConsole(void)
{
 FreeConsole(); 
//на всякий случай 
if ( AllocConsole() ) 
{ 
int hCrt = _open_osfhandle((long)
 GetStdHandle(STD_OUTPUT_HANDLE), _O_TEXT);
 *stdout = *(::_fdopen(hCrt, "w")); 
::setvbuf(stdout, NULL, _IONBF, 0);
 *stderr = *(::_fdopen(hCrt, "w"));
 ::setvbuf(stderr, NULL, _IONBF, 0);
 return TRUE; }
return FALSE;
}


int SystemInitializate()
{
	STARTUPINFO StInfo;
	GetStartupInfo(&StInfo);
	if ((StInfo.dwFlags && 1) != 0){
		nCmdShow = StInfo.wShowWindow;
	}else{
		nCmdShow = 10;
	}
	hInstance = AfxGetInstanceHandle();
//	hInstance = AfxGetModuleHandle(NULL);
	return 0;
}




LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{

	if (message == WM_CHAR){
		if (HIWORD(lParam) == 1){
		    PostQuitMessage(0);
		}
	}
	if (message == WM_TIMER){
		KillTimer(wnd,100);
		Draw_Timer_Counter++;
		InvalidateRect(wnd, NULL, 0); 
		UpdateWindow(wnd);
		InvalidateRect(wnd_earth, NULL, 0); 
		UpdateWindow(wnd_earth);
/*		if (Draw_Timer_Counter > 10){
//			InvalidateRect(wnd_cap, NULL, 0); 
			InvalidateRect(wnd_earth, NULL, 0); 
		}
		if (Draw_Timer_Counter > 20){
			InvalidateRect(wnd_graph, NULL, 0); 
			Draw_Timer_Counter = 0;
		} */
		SetTimer (wnd, 100, 20, NULL);
		Angle = Angle+0.1;
	}

	BOOL b = 1;
	if (hWnd == wnd){
		return WndProc_3d_view(hWnd, message, wParam, lParam);
		b = 0;
//		exit(0);
	}
	if (hWnd == wnd_earth){
		return WndProc_3d_view_earth(hWnd, message, wParam, lParam);
		b = 0;
//		exit(0);
	}
//	if (hWnd == wnd_graph){
//		return WndProc_graph(hWnd, message, wParam, lParam);
//		b = 0;
//		exit(0);
//	}
	if (b == 1){
		return DefWindowProc(hWnd, message, wParam, lParam);
	}
	return 0;
}

// CAboutDlg dialog used for App About

class CAboutDlg : public CDialog
{
public:
	CAboutDlg(CWnd* parent=NULL);

// Dialog Data
	enum { IDD = IDD_ABOUTBOX };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

// Implementation
protected:
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg(CWnd* parent) : CDialog(CAboutDlg::IDD,parent)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialog)
END_MESSAGE_MAP()


// CVisualizationDlg dialog




CVisualizationDlg::CVisualizationDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CVisualizationDlg::IDD, pParent)
	, m_dec(40)
	, m_RA(18.27)
//	, m_UVtime(10000)
	, m_Lines(TRUE)
	, m_Time_scale(false)
	, m_TimeScale(10000.0)
	, m_Time(0)
//	, m_TIME(212135843659.8)
//	, start_time(0)
	, start_time2(COleDateTime::GetCurrentTime())
	, start_time1(COleDateTime::GetCurrentTime())
	, UV_time2(COleDateTime::GetCurrentTime())
	, m_UV_time_days(0)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
//	m_hIcon = AfxGetApp()->LoadIcon(IDI_ICON1);

}

void CVisualizationDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_EDIT1, m_dec);
	DDV_MinMaxDouble(pDX, m_dec, -90, 90);
	DDX_Text(pDX, IDC_EDIT2, m_RA);
	DDV_MinMaxDouble(pDX, m_RA, 0, 24);
	//	DDX_Text(pDX, IDC_EDIT3, m_UVtime);
	DDX_Check(pDX, IDC_CHECK1, m_Lines);
	DDX_Control(pDX, IDC_SPIN1, m2_Time_Scale);
	DDX_Text(pDX, IDC_EDIT4, m_TimeScale);

	//	DDX_Text(pDX, IDC_EDIT5, m_TIME);
	//	DDX_DateTimeCtrl(pDX, IDC_DATETIMEPICKER1, start_time);
	DDX_DateTimeCtrl(pDX, IDC_DATETIMEPICKER2, start_time2);
	DDX_DateTimeCtrl(pDX, IDC_DATETIMEPICKER3, start_time1);
	DDX_DateTimeCtrl(pDX, IDC_DATETIMEPICKER4, UV_time2);
	DDX_Text(pDX, IDC_EDIT3, m_UV_time_days);
	DDV_MinMaxInt(pDX, m_UV_time_days, 0, 10);
}

BEGIN_MESSAGE_MAP(CVisualizationDlg, CDialog)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	//}}AFX_MSG_MAP
	ON_WM_DESTROY()
	ON_BN_CLICKED(IDOK, &CVisualizationDlg::OnBnClickedOk)
	ON_NOTIFY(UDN_DELTAPOS, IDC_SPIN1, &CVisualizationDlg::OnDeltaposSpin1)
	ON_BN_CLICKED(IDC_BUTTON1, &CVisualizationDlg::OnBnClickedButton1)
	ON_BN_CLICKED(IDC_BUTTON2, &CVisualizationDlg::OnBnClickedButton2)
	ON_NOTIFY(DTN_DATETIMECHANGE, IDC_DATETIMEPICKER4, &CVisualizationDlg::OnDtnDatetimechangeDatetimepicker4)
END_MESSAGE_MAP()


// CVisualizationDlg message handlers

BOOL CVisualizationDlg::OnInitDialog()
{


	CDialog::OnInitDialog();

	// IDM_ABOUTBOX must be in the system command range.
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != NULL)
	{
		CString strAboutMenu;
		strAboutMenu.LoadString(IDS_ABOUTBOX);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// Set the icon for this dialog.  The framework does this automatically
	//  when the application's main window is not a dialog
	SetIcon(m_hIcon, TRUE);			// Set big icon
	SetIcon(m_hIcon, FALSE);		// Set small icon

	// TODO: Add extra initialization here

	CreateConsole();
    SystemInitializate();

	MSG msg;
	WNDCLASS wcex;

	wcex.style			= CS_HREDRAW | CS_VREDRAW | CS_OWNDC;//CS_HREDRAW | CS_VREDRAW;
	wcex.lpfnWndProc	= WndProc;
	wcex.cbClsExtra		= 0;
	wcex.cbWndExtra		= 0;
	wcex.hInstance		= hInstance;
	wcex.hIcon			= LoadIcon(NULL, IDI_APPLICATION);
	wcex.hCursor		= LoadCursor(NULL, IDC_ARROW);
	wcex.hbrBackground	= 0;//(HBRUSH)(COLOR_WINDOW+1);
	wcex.lpszMenuName	= NULL;
	wcex.lpszClassName	= "MyClass";

	RegisterClass(&wcex);

	printf("*\n");

	int screenW=GetSystemMetrics(SM_CXSCREEN);//Получить ширину экрана
	int  screenH=GetSystemMetrics(SM_CYSCREEN);//Получить высоту  экрана

	wnd = CreateWindow("MyClass", "3D_View", WS_OVERLAPPEDWINDOW,0,0, screenW*2/3, screenH, 0, 0, hInstance, NULL);
    if (!wnd)
    {
		return FALSE;
	    exit(1);
    } 
	wnd_earth = CreateWindow("MyClass", "Earth", WS_OVERLAPPEDWINDOW,screenW*2/3,0, screenW/3, screenH/2, 0, 0, hInstance, NULL);
    if (!wnd_earth)
    {
		return FALSE;
	    exit(1);
    } 

	InitWndGraphSettings();
//	::ShowWindow(wnd, SW_SHOWNORMAL);
	::ShowWindow(wnd, nCmdShow);
	::UpdateWindow(wnd);

//	::ShowWindow(wnd_earth, SW_SHOWNORMAL);
	::ShowWindow(wnd_earth, nCmdShow);
	::UpdateWindow(wnd_earth);


    ::SetTimer (wnd, 100, 1, NULL);


	double t4,t5,t6,t7;
	int t1,t2,t3;
	JulDateToDate(Time/24.0/60.0/60.0, t3, t2, t1);
	t4 = (Time/24.0/60.0/60.0-0.5 -(int)(Time/24.0/60.0/60.0 -0.5))*24;
	t5 = (t4 - (int)(t4))*60;
	t6 = (t5 - (int)(t5))*60;
	start_time1.SetDate(t1,t2,t3);
	start_time2.SetTime(t4,t5,t6);
	t7 = UVtime/60/60/24;
//	t4 = UVtime/60/60;
	t4 = (t7 - (int)(t7))*24;
	t5 = (t4 - (int)(t4))*60;
	t6 = (t5 - (int)(t5))*60;
	m_UV_time_days = t7;
	UV_time2.SetTime(t4,t5,t6);
	UpdateData(false);

	return TRUE;  // return TRUE  unless you set the focus to a control
}

void CVisualizationDlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialog::OnSysCommand(nID, lParam);
	}
}

// If you add a minimize button to your dialog, you will need the code below
//  to draw the icon.  For MFC applications using the document/view model,
//  this is automatically done for you by the framework.

void CVisualizationDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // device context for painting

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// Center icon in client rectangle
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// Draw the icon
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialog::OnPaint();
	}
}

// The system calls this function to obtain the cursor to display while the user drags
//  the minimized window.
HCURSOR CVisualizationDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}


void CVisualizationDlg::OnDestroy()
{
	::KillTimer(wnd,100);
	CDialog::OnDestroy();
}

void CVisualizationDlg::OnBnClickedOk()
{
	if (!UpdateData(true) )  {
		return ;
	}

	printf("%f \n",sigma_source);
	sigma_source = m_dec/180*Pi;
	H_source = (24-m_RA) / 24 * 2.0*Pi;
	printf("%f \n",sigma_source);
	UVtime = m_UV_time_days*24*3600 + UV_time2.GetHour()*3600+UV_time2.GetMinute()*60+UV_time2.GetSecond(); //m_UVtime;
//	UV_time2.GetHour()*3600+UV_time2.GetMinute()*60+UV_time2.GetSecond();
//	m_UV_time_days
	if (m_Lines)
	{
		DRAW_EKL = 1;
	}else {
		DRAW_EKL = 0;
	}
	double t1,t2,t3,t4,t5,t6;
	t1 = start_time1.GetYear();
	t2 = start_time1.GetMonth();
	t3 = start_time1.GetDay();
	t4 = start_time2.GetHour();
	t5 = start_time2.GetMinute();
	t6 = start_time2.GetSecond();
//	printf("fjgfdjklgnfjgnflgndflkjgnlkfdgnlds %f %f %f %f %f %f \n",t1, t2, t3, t4, t5, t6);
	Time = Date_JDate(t3, t2, t1, (t4*60+t5)*60+t6) *24.0*60.0*60.0;;
//	Time = m_TIME;
//	n_max = UVtime/dTime*1000;  //&&&&&&
	n = 1;
	m_TimeScale = dTime*10.0;
	UpdateData(false);
}

void CVisualizationDlg::OnDeltaposSpin1(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMUPDOWN pNMUpDown = reinterpret_cast<LPNMUPDOWN>(pNMHDR);
	printf("%i \n",pNMUpDown->iDelta);
	if (pNMUpDown->iDelta == -1){
		dTime = dTime*2.0;
	}else {
		dTime = dTime/2.0;
	}
	m_TimeScale = dTime*10.0;
	UpdateData(false);

//	n_max = UVtime/dTime*1000; //&&&&&&
	n = 1;
//	mas[n]
	printf("N %i \n", n_max);
	// TODO: Add your control notification handler code here
	*pResult = 0;
}

void CVisualizationDlg::OnBnClickedButton1()
{
	if (dTime < 0.001){
		dTime = Time_old;
	}else{
		dTime = 0;
	}
	// TODO: Add your control notification handler code here
}

void CVisualizationDlg::OnBnClickedButton2()
{
		CAboutDlg dlgAbout(this);
		dlgAbout.DoModal();
//		dlgAbout.;

	// TODO: Add your control notification handler code here
}

void CVisualizationDlg::OnDtnDatetimechangeDatetimepicker4(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMDATETIMECHANGE pDTChange = reinterpret_cast<LPNMDATETIMECHANGE>(pNMHDR);
	// TODO: Add your control notification handler code here
	*pResult = 0;
}

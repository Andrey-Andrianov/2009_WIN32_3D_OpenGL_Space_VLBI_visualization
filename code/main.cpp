#include <windows.h>
#include <stdio.h>

#include <gl/gl.h>
//#include <gl/glu.h>
//#include <gl/glaux.h>

#include "constants.h"
#include "interface.h"
#include "3d_view.h"
#include "cap.h"
#include "graph.h"
#include "logic.h"
#include "math.h"

#include <io.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <share.h>

//HINSTANCE hInstance;
typedef  GLuint (* _stdcall MYPROC)(GLenum); 
//extern  int _stdcall   glArrayElement (GLint index); 
//extern  int _stdcall   glCreateShader (GLint index); 
//WINGDIAPI int APIENTRY glCreateShader (int type);
//PROC glCreateShader;
PROC glaaa;
//MYPROC glCreateShader;
GLuint (__stdcall *glCreateShader)(GLenum);

//void (*glaaa)();
#define GL_VERTEX_SHADER                     0x8B31


int nCmdShow;
int TimerID;

int Draw_Timer_Counter = 0;


int SystemInitializate()
{
	STARTUPINFO StInfo;
	GetStartupInfo(&StInfo);
	if ((StInfo.dwFlags && 1) != 0){
		nCmdShow = StInfo.wShowWindow;
	}else{
		nCmdShow = 10;
	}
	hInstance = GetModuleHandle(NULL);
	return 0;
}

void CALLBACK MyTimerFunc(UINT uID, UINT uMsg, DWORD dwUser, DWORD dw1, DWORD dw2)
{
	if (BUSY == 0){
		for (int i=0;i<100;i++){
			TimerProc();         
   		}
	}
}


LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
//	if (message == WM_PAINT){
//		if (BUSY == 0){
//			for (int i=0;i<1000;i++){
//				TimerProc();         // ��� ������, �� ��� ������
//    		}
//		}
//	}
	if (message == WM_CHAR){
		KeysCommand((char)wParam);
 
		if (HIWORD(lParam) == 1){
		    PostQuitMessage(0);
		}
	}
	if (message == WM_TIMER){
		KillTimer(wnd,100);
/*		if (BUSY == 0){
			for (int i=0;i<100;i++){
				TimerProc();         
    		}
		}*/
		Draw_Timer_Counter++;
		InvalidateRect(wnd, NULL, 0); 
		if (Draw_Timer_Counter > 10){
			InvalidateRect(wnd_cap, NULL, 0); 
		}
		if (Draw_Timer_Counter > 20){
			InvalidateRect(wnd_graph, NULL, 0); 
			Draw_Timer_Counter = 0;
		}
		SetTimer (wnd, 100, 1, NULL);
		Angle = Angle+0.1;
	}

	BOOL b = 1;
	if (hWnd == wnd){
		return WndProc_3d_view(hWnd, message, wParam, lParam);
		b = 0;
//		exit(0);
	}
	if (hWnd == wnd_cap){
		return WndProc_cap(hWnd, message, wParam, lParam);
		b = 0;
//		exit(0);
	}
	if (hWnd == wnd_graph){
		return WndProc_graph(hWnd, message, wParam, lParam);
		b = 0;
//		exit(0);
	}
	if (b == 1){
		return DefWindowProc(hWnd, message, wParam, lParam);
	}
	return 0;
}
//LPTIMECALLBACK
//void TimerProc(UINT uTimerid, UINT umessage, DWORD_PTR dwuser, DWORD_PTR dw1, DWORD_PTR dw2) 

//{//TIMECALLBACK
//	return 0;
//	int b = 10;
//}


//int _tmain(int argc, _TCHAR* argv[])
int main(int argc, char* argv[])
//int WINAPI WinMain(HINSTANCE hInstance,HINSTANCE hPrevInstance,LPTSTR lpCmdLine,int nCmdShow)

{
    SystemInitializate();
	draw_pixels_light = (DRAW_PIXELS*)malloc(sizeof(DRAW_PIXELS));
	rays = (RAYS*)malloc(sizeof(RAYS));


	InitPulsarParameters();
	ResetPulsarParameters();

	PrintPulsarParameters();

	for(int i=0;i<2000;i++){
			double r_sm = abs((1000.0-i)/1000.0)*2.0;
			draw_pulse2[i].x = 1.0 * exp(-r_sm*r_sm) / (1.0+pow(0.0/r_sm,5.0));
			draw_pulse2[i].y = 1.0 * exp(-r_sm*r_sm) / (1.0+pow(0.5/r_sm,5.0));
//			draw_pulse2[i].x = 0.33;// printf("%i%e\n",mm,draw_pulse2[mm].x);
	}
	draw_pulse2[1000].x = 1.0;
//	TimerProc& a = TimerProc;
//	a =  TimerProc;
//	void (*a)(UINT uTimerid, UINT umessage, DWORD_PTR dwuser, DWORD_PTR dw1, DWORD_PTR dw2);
//	a = TimerProc;
    //id_timer = timeSetEvent(1,5,(LPTIMECALLBACK)a,0,TIME_PERIODIC);  

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

	wnd_graph = CreateWindow("MyClass", "Graph", WS_OVERLAPPEDWINDOW,0,512, 200, 200+28, 0, 0, hInstance, NULL);
    if (!wnd_graph)
    {
		return FALSE;
	    exit(1);
    } 
	wnd_cap = CreateWindow("MyClass", "Cap", WS_OVERLAPPEDWINDOW,400,512, 200, 200+28, 0, 0, hInstance, NULL);
    if (!wnd_cap)
    {
		return FALSE;
	    exit(1);
    } 
	wnd = CreateWindow("MyClass", "3D_View", WS_OVERLAPPEDWINDOW,0,0, 512, 500, 0, 0, hInstance, NULL);
    if (!wnd)
    {
		return FALSE;
	    exit(1);
    } 

    mainDlg = CreateDialog(hInstance, MAKEINTRESOURCE(101), 0, MainDlgProc);
    ShowWindow(mainDlg, nCmdShow);
	UpdateWindow(mainDlg);

	InitWndGraphSettings();
    ShowWindow(wnd, nCmdShow);
	UpdateWindow(wnd);

    ShowWindow(wnd_cap, nCmdShow);
	UpdateWindow(wnd_cap);

    //ShowWindow(wnd_graph, nCmdShow);
	//UpdateWindow(wnd_graph);

    SetTimer (wnd, 100, 1, NULL);
    TimerID = timeSetEvent(10,5,MyTimerFunc,0,TIME_PERIODIC);

	while (GetMessage(&msg, NULL, 0, 0))
	{
		TranslateMessage(&msg);
		DispatchMessage(&msg);
	}

	KillTimer(wnd,100);
    timeKillEvent(TimerID);

	free(rays);
	free(draw_pixels_light);
	return (int) msg.wParam;
}

#include <windows.h>
#include "constants.h"

double Angle;                 // ”гол поворота пульсара относительго оси вращени€ (градусы)
double Angle_naklona;         // ”гол между осью вращени€ и магнитной осью (градусы)

HINSTANCE hInstance;

HDC dc;			
HWND wnd;
HGLRC HRC;

HDC dc_earth;			
HWND wnd_earth;
HGLRC HRC_earth;



HWND mainDlg;


int WndWidth(HWND window)
{
    RECT r;
    long nWindowWidth = 100;
	GetClientRect(window, &r);
	nWindowWidth = r.right;
	return nWindowWidth;
}

int WndHeight(HWND window)
{
	RECT r;
	long nWindowHeight = 100;
	GetClientRect(window, &r);
	nWindowHeight = r.bottom;
	return nWindowHeight;
}



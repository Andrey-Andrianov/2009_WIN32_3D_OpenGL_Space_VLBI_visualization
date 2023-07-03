extern double Angle;                 // ”гол поворота пульсара относительго оси вращени€ (градусы)
extern double Angle_naklona;         // ”гол между осью вращени€ и магнитной осью (градусы)
#define Pi 3.1415926535897932384626433832795

struct _VECTOR {
	double x;
	double y;
	double z;
};
typedef struct  _VECTOR VECTOR;

struct _VECTOR_single {
	double x;
	double y;
	double z;
};
typedef struct  _VECTOR_single VECTOR_single;

extern HINSTANCE hInstance;

extern HDC dc;			
extern HWND wnd;
extern HGLRC HRC;

extern HDC dc_earth;			
extern HWND wnd_earth;
extern HGLRC HRC_earth;

extern HWND mainDlg;

int WndWidth(HWND window);
int WndHeight(HWND window);


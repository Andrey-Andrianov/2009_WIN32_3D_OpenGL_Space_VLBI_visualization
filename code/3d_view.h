extern HGLRC HRC;
extern double H_source, sigma_source;
extern double UVtime;
extern int n_max;
extern int n;
//extern AUV *mas;

extern double Time, dTime,Time_old;
extern int DRAW_EKL;


int SetDcPixelFormat(HDC dc);
int InitWndGraphSettings();
LRESULT CALLBACK WndProc_3d_view(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam);
LRESULT CALLBACK WndProc_3d_view_earth(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam);


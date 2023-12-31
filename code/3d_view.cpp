#include "AUVPlaneWnd.h"

#include <windows.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include <gl/gl.h>
#include <gl/glu.h>
#include <vector>
//#include <gl/glut.h>
//#include <gl/glaux.h>

#include "constants.h"
#include "3d_view.h"
//#include "logic.h"
//#include "interface.h"
#include "ARandom.h"
#include "Igmath.h"
//#include "AUVPlaneWnd.h"

//HGLRC HRC;

	GLuint lists[1000];
	GLfloat LightPos[4];
	GLfloat global_ambient[4];
	unsigned char tex[512][512][4];
	

double angle_x_old = 0;//30;
double angle_y_old = 0;
double angle_z_old = 0;

#define V_light 299792458.0
#define GLF_START_LIST 500
#define T0  40459.8
double H_source = 1.5;
double sigma_source = 0.7;
double UVtime = 100000;


double CamX, CamY, CamZ, CamAngle, CamAngleNaklona,CamAlpha, CamTeta;
double CamX_Earth, CamY_Earth, CamZ_Earth, CamAngle_Earth, CamAngleNaklona_Earth,CamAlpha_Earth, CamTeta_Earth;
double EarthX, EarthY, EarthZ, EarthAngle, EarthAngleNaklona;
double LunaX, LunaY, LunaZ, LunaAngle, LunaAngleNaklona;
double SunX, SunY, SunZ, SunAngle, SunAngleNaklona;
double RadioAstronX, RadioAstronY, RadioAstronZ, RadioAstronAngle, RadioAstronAngleNaklona;
double R_earth, R_luna;

double Time,dTime,Time_old;
double scale_km;
int orbita_last;
int orbita_current;
double MouseX,MouseY;
int MouseKey;
double MouseX_Earth,MouseY_Earth;
int MouseKey_Earth;

struct _ORB {
	VECTOR r;
	double time;
};
typedef struct  _ORB ORB;
struct _VLINE {
	int v1;
	int v2;
};
typedef struct  _VLINE VLINE;
struct _VTRIANGL {
	int v1;
	int v2;
	int v3;
	VECTOR Normal;
};
typedef struct  _VTRIANGL VTRIANGL;
struct _VSQUARE {
	int v1;
	int v2;
	int v3;
	int v4;
	VECTOR Normal;
};
typedef struct  _VSQUARE VSQUARE;
ORB orbita[15000];
VECTOR Vertexs[91000];
VLINE Lines[10000];
VTRIANGL Triangles[30000];
VSQUARE Squares[30000];
int Num_Triangles, Num_Squares, Num_Lines;
int DRAW_EKL;


#pragma pack(push,1)
struct _vert {
	float x;
	float y;
	float z;
};
typedef struct  _vert vert;
struct _pol {
	__int16 v1;
	__int16 v2;
	__int16 v3;
	__int16 f;
};
typedef struct  _pol pol;
struct _tex2 {
	float u;
	float v;
};
typedef struct  _tex2 tex2;
struct _obj {
	int v;
	int t;
	int p;
	int parents;
	char name[30];
	float x_rot;
	float y_rot;
	float z_rot;
};
typedef struct  _obj obj;
vert  vertixs[80000];
vert  normals[80000];
tex2  textures[80000];
pol   poligons[80000];
obj   objekts[10000];
#pragma pack(pop)


int num_objekts_vert;
int num_objekts_polig;
int num_objekts_tex;
int num_objekts;
int num_vertixs;
int num_poligons;
int num_textures;

struct _CTELESCOP {
	char Name[40];
	char Name2[5];
	double X;
	double Y;
	double Z;
	double Diam;
	double AZmin;
	double AZmax;
	double Zmin;
	double Zmax;
	int flag;
};
typedef struct  _CTELESCOP CTELESCOP;

CTELESCOP CoordTelescops[100];
int NumTelescops;


int n = 1;
int n_max = 10000000;
AUV *mas = new AUV[n_max];
double umax = 500e6;
double vmax = 500e6;
AUVPlaneWnd	m_wnd;
//double temp_time;


int windparamX,windparamY,windparamXold,windparamYold,windparamwindW,windparamwindH;
int  NewCount, LastCount, FrameCount;

char szString[17];

vector<VECTOR_single> Telescope_Vert;
vector<VECTOR_single> Telescope_Triangles;
vector<VECTOR_single> Telescope_Normals;

int SetDcPixelFormat(HDC dc)
{
//Telescope_Vert.
   //
   // Fill in the pixel format descriptor.
   //
   PIXELFORMATDESCRIPTOR pfd ;
   memset(&pfd, 0, sizeof(PIXELFORMATDESCRIPTOR)) ;
   pfd.nSize      = sizeof(PIXELFORMATDESCRIPTOR); 
   pfd.nVersion   = 1 ; 
   pfd.dwFlags    = PFD_DOUBLEBUFFER |
                    PFD_SUPPORT_OPENGL |
                    PFD_DRAW_TO_WINDOW ;
   pfd.iPixelType = PFD_TYPE_RGBA;
   pfd.cColorBits = 32;
   pfd.cDepthBits = 32; 
   pfd.iLayerType = PFD_MAIN_PLANE ;

   int nPixelFormat = ChoosePixelFormat(dc, &pfd);
   BOOL bResult = SetPixelFormat (dc, nPixelFormat, &pfd);

   return 0;
}

int OutText (char Litera[40])
{
	glPushAttrib(GL_LIST_BIT);
	glListBase(GLF_START_LIST);
	glCallLists((unsigned short int)Litera[0], GL_UNSIGNED_BYTE, Litera+1);
		//glListBase(0);
	glPopAttrib();
	return 0;
}

int read_telescopes() // 0 - radioastron
{

	FILE* file;
	char str[256];// = "123456";
	printf("1111111\n");
//	POLYCO Pol;

	char str2[256];
	char str3[256];
	file = fopen("data\\telescopes.dat","r");

	int N;
	for (int i = 1; i < 11; i++){

		fscanf(file,"%s",str);
		N = atoi(str);
		fscanf(file,"%s",CoordTelescops[i].Name2);
		fscanf(file,"%s",CoordTelescops[i].Name);
		fscanf(file,"%s",str);
		CoordTelescops[i].X = atof(str);
		fscanf(file,"%s",str);
		CoordTelescops[i].Y = -atof(str);
		fscanf(file,"%s",str);
		CoordTelescops[i].Z = atof(str);
		fscanf(file,"%s",str);
		CoordTelescops[i].Diam = atof(str);
		fscanf(file,"%s",str);
		CoordTelescops[i].AZmin = atof(str);
		fscanf(file,"%s",str);
		CoordTelescops[i].AZmax = atof(str);
		fscanf(file,"%s",str);
		CoordTelescops[i].Zmin = atof(str);
		fscanf(file,"%s",str);
		CoordTelescops[i].Zmax = atof(str);
		fscanf(file,"%s",str);
		N = atoi(str);

	printf("%i %s %s %e %e %e %e %e %e %e %e %i \n",N,CoordTelescops[i].Name2,CoordTelescops[i].Name,CoordTelescops[i].X,CoordTelescops[i].Y,CoordTelescops[i].Z,
			CoordTelescops[i].Diam,CoordTelescops[i].AZmin,CoordTelescops[i].AZmax,CoordTelescops[i].Zmin,CoordTelescops[i].Zmax,N);
	//CoordTelescops[i].X = random(2.0)-1.0;
	//CoordTelescops[i].Y = random(2.0)-1.0;
	//CoordTelescops[i].Z = random(2.0)-1.0;
	}



	fclose(file);

	NumTelescops = 11;
	return 0;
}

VECTOR CalculateUVPoint(int T1, int T2, double H, double sigma, double Freq)
{
	VECTOR uv;
	double u,v;
	VECTOR D,v_u,v_v;
	D.x = CoordTelescops[T1].X - CoordTelescops[T2].X;
	D.y = CoordTelescops[T1].Y - CoordTelescops[T2].Y;
	D.z = CoordTelescops[T1].Z - CoordTelescops[T2].Z;
	v_u.x = sin(H); v_u.y = cos(H); v_u.z = 0;
	v_v.x = -cos(H)*sin(sigma); v_v.y = sin(H)*sin(sigma); v_v.z = cos(sigma);
	u = Freq / V_light * (D.x*v_u.x+D.y*v_u.y+D.z*v_u.z);
	v = Freq / V_light * (D.x*v_v.x+D.y*v_v.y+D.z*v_v.z);
	uv.x = u; uv.y =v;
	return uv;
}

int read_textures()
{
  FILE *id_file;
  int i,j;
  unsigned char buffer [512][512][3];
  unsigned char buffer3 [256][256][3];
  unsigned char* buffer2 =  (unsigned char*)malloc(2048*2048*3);
  unsigned char* tex2 =  (unsigned char*)malloc(2048*2048*4);
  unsigned char* tex3 =  (unsigned char*)malloc(256*256*4);
 
  id_file=fopen("data\\nebo5.bmp","r+t");
  fseek(id_file,54,0);
  fread(buffer2,2048*2048*3,1,id_file);
  fclose(id_file);
  for(i = 0;i < 2048;i++)
  {
     for(j = 0;j < 2048;j++)
     {
		 tex2[(j*2048+i)*4+0] = buffer2[(j*2048+i)*3+2];
	     tex2[(j*2048+i)*4+1] = buffer2[(j*2048+i)*3+1];
	     tex2[(j*2048+i)*4+2] = buffer2[(j*2048+i)*3+0];
	 }
  }
	lists[0] = glGenLists(1);
	glNewList(lists[0],GL_COMPILE);

	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
//	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
//	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
	glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA,2048,2048,0,GL_RGBA,GL_UNSIGNED_BYTE,tex2);
	glEndList();
  free(tex2);
  free(buffer2);
  id_file=fopen("data\\earth4.bmp","r+t");
  fseek(id_file,54,0);
  fread(buffer,sizeof(buffer),1,id_file);
  fclose(id_file);
  for(i = 0;i < 512;i++)
  {
     for(j = 0;j < 512;j++)
     {
	   tex[j][i][0] = buffer[j][i][2];
	   tex[j][i][1] = buffer[j][i][1];
	   tex[j][i][2] = buffer[j][i][0];
	 }
  }
	lists[1] = glGenLists(1);
	glNewList(lists[1],GL_COMPILE);

	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
//	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
//	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
	glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA,512,512,0,GL_RGBA,GL_UNSIGNED_BYTE,&tex);
	glEndList();

  id_file=fopen("data\\luna.bmp","r+t");
  fseek(id_file,54,0);
  fread(buffer,sizeof(buffer),1,id_file);
  fclose(id_file);
  for(i = 0;i < 512;i++)
  {
     for(j = 0;j < 512;j++)
     {
	   tex[j][i][0] = buffer[j][i][2];
	   tex[j][i][1] = buffer[j][i][1];
	   tex[j][i][2] = buffer[j][i][0];
	 }
  }
	lists[2] = glGenLists(1);
	glNewList(lists[2],GL_COMPILE);

	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
//	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
//	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
	glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA,512,512,0,GL_RGBA,GL_UNSIGNED_BYTE,&tex);
	glEndList();

  id_file=fopen("data\\sun14.bmp","r+t");
  fseek(id_file,54,0);
  fread(buffer3,sizeof(buffer3),1,id_file);
  fclose(id_file);
  for(i = 0;i < 256;i++)
  {
     for(j = 0;j < 256;j++)
     {
	   tex3[(j*256+i)*4+0] = buffer3[j][i][2];
	   tex3[(j*256+i)*4+1] = buffer3[j][i][1];
	   tex3[(j*256+i)*4+2] = buffer3[j][i][0];
//	   tex3[(j*256+i)*4+0] = buffer3[j][i][2];
//	   if ( (double)(buffer3[j][i][1]+100) < 255 ){
//	   tex3[(j*256+i)*4+1] = buffer3[j][i][1]+100;
//	   }else{
//	   tex3[(j*256+i)*4+1] = 255;
//	   }
//	   tex3[(j*256+i)*4+2] = buffer3[j][i][0];
//	   tex3[(j*256+i)*4+0] = (buffer3[j][i][2]+buffer3[j][i][1])*0.8;
//	   tex3[(j*256+i)*4+1] = (buffer3[j][i][1]+buffer3[j][i][2])*0.8;
//	   tex3[(j*256+i)*4+2] = 0;//buffer3[j][i][0];

//	   tex3[(j*256+i)*4+1] = tex3[(j*256+i)*4+0];
	  // if ((tex3[(j*256+i)*4+2]+tex3[(j*256+i)*4+1]+tex3[(j*256+i)*4+0]) > 230*3 ){
		//   tex3[(j*256+i)*4+3] = 0;
	   //} else {
		//   tex3[(j*256+i)*4+3] = 255;
	  // }
	   tex3[(j*256+i)*4+3] = (tex3[(j*256+i)*4+2]+tex3[(j*256+i)*4+1]+tex3[(j*256+i)*4+0])/(double)(3);
	 }
  }
	lists[3] = glGenLists(1);
	glNewList(lists[3],GL_COMPILE);

	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
//	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
//	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
	glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA,256,256,0,GL_RGBA,GL_UNSIGNED_BYTE,tex3);
	glEndList();

  id_file=fopen("data\\text1.bmp","r+t");
  fseek(id_file,54,0);
  fread(buffer3,sizeof(buffer3),1,id_file);
  fclose(id_file);
  for(i = 0;i < 256;i++)
  {
     for(j = 0;j < 256;j++)
     {
	   tex3[(j*256+i)*4+0] = buffer3[j][i][2];
	   tex3[(j*256+i)*4+1] = buffer3[j][i][1];
	   tex3[(j*256+i)*4+2] = buffer3[j][i][0];
	   if ((tex3[(j*256+i)*4+2]+tex3[(j*256+i)*4+1]+tex3[(j*256+i)*4+0]) > 0 ){
		   tex3[(j*256+i)*4+3] = 255;
	   } else {
		   tex3[(j*256+i)*4+3] = 0;
	   }
	 }
  }
	lists[4] = glGenLists(1);
	glNewList(lists[4],GL_COMPILE);

	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
//	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
//	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
	glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA,256,256,0,GL_RGBA,GL_UNSIGNED_BYTE,tex3);
	glEndList();

  id_file=fopen("data\\text2.bmp","r+t");
  fseek(id_file,54,0);
  fread(buffer3,sizeof(buffer3),1,id_file);
  fclose(id_file);
  for(i = 0;i < 256;i++)
  {
     for(j = 0;j < 256;j++)
     {
	   tex3[(j*256+i)*4+0] = buffer3[j][i][2];
	   tex3[(j*256+i)*4+1] = buffer3[j][i][1];
	   tex3[(j*256+i)*4+2] = buffer3[j][i][0];
	   if ((tex3[(j*256+i)*4+2]+tex3[(j*256+i)*4+1]+tex3[(j*256+i)*4+0]) > 0 ){
		   tex3[(j*256+i)*4+3] = 255;
	   } else {
		   tex3[(j*256+i)*4+3] = 0;
	   }
	 }
  }
	lists[5] = glGenLists(1);
	glNewList(lists[5],GL_COMPILE);

	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
//	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
//	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
	glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA,256,256,0,GL_RGBA,GL_UNSIGNED_BYTE,tex3);
	glEndList();

  id_file=fopen("data\\text3.bmp","r+t");
  fseek(id_file,54,0);
  fread(buffer3,sizeof(buffer3),1,id_file);
  fclose(id_file);
  for(i = 0;i < 256;i++)
  {
     for(j = 0;j < 256;j++)
     {
	   tex3[(j*256+i)*4+0] = buffer3[j][i][2];
	   tex3[(j*256+i)*4+1] = buffer3[j][i][1];
	   tex3[(j*256+i)*4+2] = buffer3[j][i][0];
	   if ((tex3[(j*256+i)*4+2]+tex3[(j*256+i)*4+1]+tex3[(j*256+i)*4+0]) > 0 ){
		   tex3[(j*256+i)*4+3] = 255;
	   } else {
		   tex3[(j*256+i)*4+3] = 0;
	   }
	 }
  }
	lists[6] = glGenLists(1);
	glNewList(lists[6],GL_COMPILE);

	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
//	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
//	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
	glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA,256,256,0,GL_RGBA,GL_UNSIGNED_BYTE,tex3);
	glEndList();

  id_file=fopen("data\\text4.bmp","r+t");
  fseek(id_file,54,0);
  fread(buffer3,sizeof(buffer3),1,id_file);
  fclose(id_file);
  for(i = 0;i < 256;i++)
  {
     for(j = 0;j < 256;j++)
     {
	   tex3[(j*256+i)*4+0] = buffer3[j][i][2];
	   tex3[(j*256+i)*4+1] = buffer3[j][i][1];
	   tex3[(j*256+i)*4+2] = buffer3[j][i][0];
	   if ((tex3[(j*256+i)*4+2]+tex3[(j*256+i)*4+1]+tex3[(j*256+i)*4+0]) > 0 ){
		   tex3[(j*256+i)*4+3] = 255;
	   } else {
		   tex3[(j*256+i)*4+3] = 0;
	   }
	 }
  }
	lists[7] = glGenLists(1);
	glNewList(lists[7],GL_COMPILE);

	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
//	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
//	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
	glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA,256,256,0,GL_RGBA,GL_UNSIGNED_BYTE,tex3);
	glEndList();

  id_file=fopen("data\\source.bmp","r+t");
  fseek(id_file,54,0);
  fread(buffer3,sizeof(buffer3),1,id_file);
  fclose(id_file);
  for(i = 0;i < 256;i++)
  {
     for(j = 0;j < 256;j++)
     {
	   tex3[(j*256+i)*4+0] = buffer3[j][i][2];
	   tex3[(j*256+i)*4+1] = buffer3[j][i][1];
	   tex3[(j*256+i)*4+2] = buffer3[j][i][0];
	   if ((tex3[(j*256+i)*4+2]+tex3[(j*256+i)*4+1]+tex3[(j*256+i)*4+0]) > 240*3 ){
		   tex3[(j*256+i)*4+3] = 0;
	   } else {
		   tex3[(j*256+i)*4+3] = 255;
	   }
	 }
  }
	lists[8] = glGenLists(1);
	glNewList(lists[8],GL_COMPILE);

	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
//	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
//	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
	glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA,256,256,0,GL_RGBA,GL_UNSIGNED_BYTE,tex3);
	glEndList();

  id_file=fopen("data\\antenna3.bmp","r+t");
  fseek(id_file,54,0);
  fread(buffer,sizeof(buffer),1,id_file);
  fclose(id_file);
  for(i = 0;i < 512;i++)
  {
     for(j = 0;j < 512;j++)
     {
	   tex[j][i][0] = buffer[j][i][2];
	   tex[j][i][1] = buffer[j][i][1];
	   tex[j][i][2] = buffer[j][i][0];
	 }
  }
	lists[9] = glGenLists(1);
	glNewList(lists[9],GL_COMPILE);

	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
//	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
//	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
	glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA,512,512,0,GL_RGBA,GL_UNSIGNED_BYTE,&tex);
	glEndList();


	return 0;
}

int	read_orbita()
{
	FILE* file;
	char str[256];// = "123456";
	printf("1111111\n");
//	POLYCO Pol;

	char str2[256];
	char str3[256];
	file = fopen("data\\r2011_07_18.txt","r");
//	file = fopen("data\\RA2003.txt","r");
	fgets(str,250,file);

	int N = 0;
	double t_y,t_m,t_d,t_h,t_min,t_s;
  char* end_file = (char*)1;
  do{


	end_file = fgets(str,250,file);
//	printf("adfafsdfgsd %i\n",N);
	if (end_file != NULL){
		//printf("%s",str);

		for (int i = 0; i < 23;i++){str2[i] = str[i];} str2[24] = '\0'; //str2[0] = '4';
		str3[0] = str2[0];
		str3[1] = str2[1];
		str3[2] = str2[2];
		str3[3] = str2[3];
		str3[4] = '\0';
		t_y = atof(str3);
		str3[0] = str2[5];
		str3[1] = str2[6];
		str3[2] = '\0';
		t_m = atof(str3);
		str3[0] = str2[8];
		str3[1] = str2[9];
		str3[2] = '\0';
		t_d = atof(str3);
		str3[0] = str2[11];
		str3[1] = str2[12];
		str3[2] = '\0';
		t_h = atof(str3);
		str3[0] = str2[14];
		str3[1] = str2[15];
		str3[2] = '\0';
		t_min = atof(str3);
		str3[0] = str2[17];
		str3[1] = str2[18];
		str3[2] = str2[19];
		str3[3] = str2[20];
		str3[4] = str2[21];
		str3[5] = str2[22];
		str3[6] = '\0';
		t_s = atof(str3);
		orbita[N].time = Date_JDate(t_d, t_m, t_y, (t_h*60+t_min)*60+t_s) *24.0*60.0*60.0;
		//printf("%s!!!! %f %f %f %f %f %f   time %f TIME %f\n",str2,t_y,t_m,t_d,t_h,t_min,t_s,orbita[N].time,((((t_y*12+t_m)*30+t_d)*24+t_h)*60+t_min)*60+t_s);

		for (int i = 0; i < 11;i++){str2[i] = str[i+23];} str2[11] = '\0'; 
			orbita[N].r.x = atof(str2);
		//printf("%f!!!!\n",orbita[N].r.x);
		for (int i = 0; i < 11;i++){str2[i] = str[i+35];} str2[11] = '\0'; 
			orbita[N].r.y = atof(str2);
		//printf("%f!!!!\n",orbita[N].r.y);
		for (int i = 0; i < 11;i++){str2[i] = str[i+47];} str2[11] = '\0'; 
			orbita[N].r.z = atof(str2);
		//printf("%f!!!!\n",orbita[N].r.z);
		//	printf("%f %f %f \n",orbita[N].r.x,orbita[N].r.y,orbita[N].r.z);

//		Sleep(1000);
				  N++;
	}

  }while (end_file != NULL);
  orbita_last = N;

  printf("\n****\n");


	fclose(file);


	return 0;
}

VECTOR CalcNormal(VECTOR v1, VECTOR v2, VECTOR v3)
{
	VECTOR norm,dv1,dv2;
	double r;
	dv1.x = v1.x-v2.x;
	dv1.y = v1.y-v2.y;
	dv1.z = v1.z-v2.z;
	dv2.x = v2.x-v3.x;
	dv2.y = v2.y-v3.y;
	dv2.z = v2.z-v3.z;
	norm.x = dv1.y*dv2.z - dv1.z*dv2.y;
	norm.y = dv1.z*dv2.x - dv1.x*dv2.z;
	norm.z = dv1.x*dv2.y - dv1.y*dv2.x;
	r = sqrt(norm.x*norm.x+norm.y*norm.y+norm.z*norm.z);
	norm.x = norm.x / r;
	norm.y = norm.y / r;
	norm.z = norm.z / r;
	return norm;
}



int	read_model()
{

//VECTOR Vertexs[76000];
//VLINE Lines[10000];
//VTRIANGL Triangles[30000];
//VSQUARE Squares[30000];
//int Num_Triangles, Num_Squares, Num_Lines;

	for (int i = 0; i < 76000;i++){
		Vertexs[i].x = 0;
		Vertexs[i].y = 0;
		Vertexs[i].z = 0;
	}

	FILE* file;
	char str[256];// = "123456";
	printf("1111111\n");
//	POLYCO Pol;

	char str2[256];
	file = fopen("data\\full_KRT.bdf","r");

	int N,N2 = 0;

	Num_Squares = 0;
	Num_Lines = 0;
	Num_Triangles = 0;
  char* end_file = (char*)1;
  do{

	end_file = fgets(str,250,file);
//	printf("adfafsdfgsd %i\n",end_file);
	if (end_file != NULL){
		if (  (str[0] == 'G')&&(str[1] == 'R')&&(str[2] == 'I')&&(str[3] == 'D') ){
		//	printf("%s",str);
			//printf("%i\n", N);
			for (int i = 0; i < 8;i++){str2[i] = str[i+8];} str2[8] = '\0'; 
			N2 = atoi(str2);
			for (int i = 0; i < 8;i++){str2[i] = str[i+24];} str2[8] = '\0'; 
			Vertexs[N2].x = atof(str2);
			for (int i = 0; i < 8;i++){str2[i] = str[i+32];} str2[8] = '\0'; 
			Vertexs[N2].y = atof(str2);
			for (int i = 0; i < 8;i++){str2[i] = str[i+40];} str2[8] = '\0'; 
			Vertexs[N2].z = atof(str2);
		    // printf("%i %f %f %f!!!!\n",N2,Vertexs[N2].x,Vertexs[N2].y,Vertexs[N2].z);

			N++;
		}
		if (  (str[0] == 'C')&&(str[1] == 'B')&&(str[2] == 'A')&&(str[3] == 'R') ){
		//	printf("%s",str);
			for (int i = 0; i < 8;i++){str2[i] = str[i+24];} str2[8] = '\0'; 
			Lines[Num_Lines].v1 = atoi(str2);
			for (int i = 0; i < 8;i++){str2[i] = str[i+32];} str2[8] = '\0'; 
			Lines[Num_Lines].v2 = atoi(str2);
		  //   printf("%i %i %i!!!!\n",Num_Lines,Lines[Num_Lines].v1,Lines[Num_Lines].v2);
			Num_Lines++;

		}
		if (  (str[0] == 'C')&&(str[1] == 'Q')&&(str[2] == 'U')&&(str[3] == 'A')&&(str[4] == 'D')&&(str[5] == '4') ){
//			printf("%s",str);
			for (int i = 0; i < 8;i++){str2[i] = str[i+24];} str2[8] = '\0'; 
			Squares[Num_Squares].v1 = atoi(str2);
			for (int i = 0; i < 8;i++){str2[i] = str[i+32];} str2[8] = '\0'; 
			Squares[Num_Squares].v2 = atoi(str2);
			for (int i = 0; i < 8;i++){str2[i] = str[i+40];} str2[8] = '\0'; 
			Squares[Num_Squares].v3 = atoi(str2);
			for (int i = 0; i < 8;i++){str2[i] = str[i+48];} str2[8] = '\0'; 
			Squares[Num_Squares].v4 = atoi(str2);
//		     printf("%i %i %i %i %i!!!!\n",Num_Squares,Squares[Num_Squares].v1,Squares[Num_Squares].v2,Squares[Num_Squares].v3,Squares[Num_Squares].v4);
			Num_Squares++;
		}
		if (  (str[0] == 'C')&&(str[1] == 'T')&&(str[2] == 'R')&&(str[3] == 'I')&&(str[4] == 'A')&&(str[5] == '3') ){
//			printf("%s",str);
			for (int i = 0; i < 8;i++){str2[i] = str[i+24];} str2[8] = '\0'; 
			Triangles[Num_Triangles].v1 = atoi(str2);
			for (int i = 0; i < 8;i++){str2[i] = str[i+32];} str2[8] = '\0'; 
			Triangles[Num_Triangles].v2 = atoi(str2);
			for (int i = 0; i < 8;i++){str2[i] = str[i+40];} str2[8] = '\0'; 
			Triangles[Num_Triangles].v3 = atoi(str2);
//		     printf("%i %i %i %i %i!!!!\n",Num_Triangles,Triangles[Num_Triangles].v1,Triangles[Num_Triangles].v2,Triangles[Num_Triangles].v3);
			Num_Triangles++;
		}
	}


  }while (end_file != NULL);

  printf("\n****\n");


	fclose(file);

	double max = 0;
	double r;
	for (int i = 0; i < 76000;i++){
		r = sqrt(Vertexs[i].x*Vertexs[i].x+Vertexs[i].y*Vertexs[i].y+Vertexs[i].z*Vertexs[i].z);
		if (r>max){max=r;}
	}
	for (int i = 0; i < 76000;i++){
		Vertexs[i].x = Vertexs[i].x/max*0.03;
		Vertexs[i].y = Vertexs[i].y/max*0.03;
		Vertexs[i].z = Vertexs[i].z/max*0.03;
//		Vertexs[i].x = Vertexs[i].x*0.1;
//		Vertexs[i].y = Vertexs[i].y*0.1;
//		Vertexs[i].z = Vertexs[i].z*0.1;
	}


		for(int i = 0; i<Num_Squares;i++){
			Squares[i].Normal = CalcNormal(Vertexs[Squares[i].v1],Vertexs[Squares[i].v2],Vertexs[Squares[i].v3]);
		} 
		for(int i = 0; i<Num_Triangles;i++){
			Triangles[i].Normal = CalcNormal(Vertexs[Triangles[i].v1],Vertexs[Triangles[i].v2],Vertexs[Triangles[i].v3]);
		}


	return 0;
}


int read_telescope_model() // 0 - radioastron
{

	char str[256];// = "123456";
  int num_objekt_animation = 0;
num_objekts_vert = 0;
num_objekts_polig = 0;
num_objekts_tex = 0;
num_objekts = 0;
num_vertixs = 0;
num_poligons = 0;
num_textures = 0;

    unsigned __int16 chunk_id;
    int chunk_leng;
    int poz = 0;
    char h;
    unsigned __int16 num;
    int i,j;	
	char Name[40] = "data\\tdrs_ant_right.3ds";
//	char Name[40] = "data\\tall_dish1.3ds";
//	char Name[40] = "data\\70m_rev_e.3ds";
	int handle,bytes;
	if((handle = open(Name,O_RDONLY|O_BINARY))== -1)
	{
		printf("Error open 3ds file.\n");
		exit(1);
	}

//	if((bytes = _read(handle,&Header,sizeof(Header))) != sizeof(Header))
//	{
//		printf("Error read in 3ds file\n");
		//gets(&c);
//		exit(1);
//	}

	do
	{
		if ( !eof(handle) ){
			if((bytes = _read(handle,&chunk_id,sizeof(chunk_id))) != sizeof(chunk_id))
			{
				printf("%i Error read in 3ds file!!!!!\n",bytes);
			}
			if((bytes = _read(handle,&chunk_leng,sizeof(chunk_leng))) != sizeof(chunk_leng))
			{
				printf("%i Error read in 3ds file!!!!!\n",bytes);
			}
		}
		printf("** %i %i \n", poz, chunk_leng);
		switch (chunk_id) {
		    case 0x4D4D:
                poz = poz + 6;
			    break;
		    case 0x3D3D:
                poz = poz + 6;
			    break;
		    case 0x4000:  //obj
                poz = poz + 6;
				i = 0;
				do
				{
					if((bytes = _read(handle,&str[i],sizeof(str[i]))) != sizeof(str[i]))
					{
						printf("%i Error read in 3ds file\n",bytes);
					}
					i++;
					poz++;
				}while ( str[i-1] != 0 );
				printf(str);
				
			    break; 
		    case 0x4100:  //trimesh-obj
                poz = poz + 6;
			    break;
		    case 0x4110:  //Vertixs
					if((bytes = _read(handle,&num,sizeof(num))) != sizeof(num))
					{
						printf("%i 11111Error read in 3ds file\n",bytes);
					}
					num_objekts_vert++;
					objekts[num_objekts_vert].v = num_vertixs;
					if((bytes = _read(handle,&vertixs[num_vertixs],sizeof(vert)*num)) != sizeof(vert)*num)
					{
						printf("%i Err555or read in 3ds file\n",bytes);
					}
		            num_vertixs = num_vertixs + num;
		            poz = poz + chunk_leng;
	                _lseek(handle,poz,0);
					printf("Vert %i All %i \n",num,num_vertixs);
			    break; 
		    case 0x4120: //Polygons
//					printf("Pol******************* %i \n",poz);
					if((bytes = _read(handle,&num,sizeof(num))) != sizeof(num))
					{
						printf("%i 111111Error read in 3ds file\n",bytes);
					}
					num_objekts_polig++;
					objekts[num_objekts_polig].p = num_poligons;
					if((bytes = _read(handle,&poligons[num_poligons],sizeof(pol)*num)) != sizeof(pol)*num)
					{
						printf("%i %i %i Erro444r read in 3ds file\n",bytes,poz,sizeof(pol));
					}
		            num_poligons = num_poligons + num;
					objekts[num_objekts_polig+1].p = num_poligons;
		            poz = poz + chunk_leng;
	                _lseek(handle,poz,0);
					printf("Pol %i All %i \n",num,num_poligons);
			    break; 
		    case 0x4140:  //Tex coord
//					printf("T1******************* %i \n",poz);
					if((bytes = _read(handle,&num,sizeof(num))) != sizeof(num))
					{
						printf("%i 111Error read in 3ds file\n",bytes);
					}
					num_objekts_tex++;
					objekts[num_objekts_tex].t = num_textures;
					if ( (bytes = _read(handle,&textures[num_textures],sizeof(tex2)*num)) != sizeof(tex2)*num )
					{
						printf("%i %i %i Error read in 3ds file\n",bytes,poz,num);
					}
		            num_textures = num_textures + num;
		            poz = poz + chunk_leng;
	                _lseek(handle,poz,0);
					printf("Tex %i All %i \n",poz,num_textures); 
//		        poz = poz + chunk_leng;
//					printf("T2******************* %i \n",poz);
  //              _lseek(handle,poz,0);
			    break;  
			default: 
		        poz = poz + chunk_leng;
                _lseek(handle,poz,0);
		}


	}while ( !eof(handle) &&( bytes != 0)  );
	//	Sleep(10000);


//	printf(str);
	_close(handle);

	double max = 0;
	double r;
	for (int i = 0; i < num_vertixs;i++){
		r = sqrt(vertixs[i].x*vertixs[i].x+vertixs[i].y*vertixs[i].y+vertixs[i].z*vertixs[i].z);
		if (r>max){max=r;}
	}
	for (int i = 0; i < num_vertixs;i++){
		vertixs[i].x = vertixs[i].x/max*0.03;
		vertixs[i].y = vertixs[i].y/max*0.03;
		vertixs[i].z = vertixs[i].z/max*0.03;
	}

	for(int i = 0; i<num_poligons;i++){
		VECTOR tt,v1,v2,v3;
		v1.x = vertixs[poligons[i].v1].x; v1.y = vertixs[poligons[i].v1].y; v1.z = vertixs[poligons[i].v1].z;
		v2.x = vertixs[poligons[i].v2].x; v2.y = vertixs[poligons[i].v2].y; v2.z = vertixs[poligons[i].v2].z;
		v3.x = vertixs[poligons[i].v3].x; v3.y = vertixs[poligons[i].v3].y; v3.z = vertixs[poligons[i].v3].z;
		tt = CalcNormal(v1,v2,v3);
		normals[i].x = tt.x; 
		normals[i].y = tt.y; 
		normals[i].z = tt.z; 
	}

	return 0;
}


/*
int read_telescope_model() // 0 - radioastron
{

	FILE* file;
	char str[256];// = "123456";
	printf("1111111\n");
//	POLYCO Pol;
	int N = 0;

	char str2[256];
	char str3[256];
    char* end_file = (char*)1;
	int group;
	VECTOR_single Temp_vector1,Temp_vector2,Temp_vector3;
	
	file = fopen("data\\telescope.dxf","r");
	do
	{
		end_file = fgets(str,250,file);
//		if(strcmp(str,"ENTITIES\n") == 0){printf(str);}
//		printf(str);
//		Sleep(1000);
	}while ( (strcmp(str,"ENTITIES\n") != 0)&&(end_file != NULL) );
	printf(str);
	while (end_file != NULL)
	{
//		end_file = (char*)fscanf(file,"%s",str);
		end_file = fgets(str,250,file);
		group = atoi(str);
		end_file = fgets(str,250,file);

//vector<VECTOR_single> Telescope_Vert;
//vector<VECTOR_single> Telescope_Triangles;
//vector<VECTOR_single> Telescope_Normals;

		switch (group) {
		    case 0:
				Telescope_Vert.push_back(Temp_vector3);
				Telescope_Vert.push_back(Temp_vector2);
				Telescope_Vert.push_back(Temp_vector1); N++;
			    break;
		    case 10:
				Temp_vector1.x = atof(str); 
//				printf("%f\n",Temp_vector1.x);
			    break;
		    case 20:
				Temp_vector1.y = atof(str); 
//				printf("%f\n",Temp_vector1.x);
			    break;
		    case 30:
				Temp_vector1.z = atof(str); 
//				printf("%f\n",Temp_vector1.x);
			    break;
		    case 11:
				Temp_vector2.x = atof(str);
//				printf("%f\n",Temp_vector1.x);
			    break;
		    case 21:
				Temp_vector2.y = atof(str);
				printf("%f\n",Temp_vector2.y);
			    break;
		    case 31:
				Temp_vector2.z = atof(str);
//				printf("%f\n",Temp_vector1.x);
			    break;
		    case 12:
				Temp_vector3.x = atof(str);
//				printf("%f\n",Temp_vector1.x);
			    break;
				Temp_vector3.y = atof(str);
//				printf("%f\n",Temp_vector1.x);
			    break;
		    case 32:
				Temp_vector3.z = atof(str);
//				printf("%f\n",Temp_vector3.z);
			    break;
		}

	}
	printf(str);
//		end_file = fgets(str,250,file);
//		printf(str);
	printf("DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD      %i",N);
	fclose(file);
	return 0;
}
*/

int DrawAntenna(float R, float G, float B)
{
/***********************************************************************************************************************/
//				glCullFace(GL_FRONT);
//				glCullFace(GL_BACK);
				VECTOR Normal;

//	if (R+G+B < 0.0001){
//		printf("fwdffgregererergv\n");
//		lists[15] = glGenLists(1);
//		glNewList(lists[15],GL_COMPILE);
				
		global_ambient[0] = 0;
		global_ambient[1] = 0;
		global_ambient[2] = 0;
		global_ambient[3] = 1;
		glLightModelfv(GL_LIGHT_MODEL_AMBIENT,global_ambient);
		LightPos[0] = 0;
		LightPos[1] = -1000;
		LightPos[2] = -1000;
		LightPos[3] = 0.0;
		glLightfv(GL_LIGHT0,GL_POSITION,LightPos);
				GLfloat MatColor[4]; 
		MatColor[0] = R; MatColor[1] = G; MatColor[2] = B; MatColor[3] = 0;
//		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, MatColor);
//		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, MatColor);
//		glEnable(GL_COLOR_MATERIAL);
		glLightfv(GL_LIGHT0, GL_AMBIENT, MatColor);
		glLightfv(GL_LIGHT0, GL_DIFFUSE, MatColor);
//		glLightfv(GL_LIGHT0, GL_AMBIENT_AND_DIFFUSE, MatColor);
        glEnable(GL_LIGHTING);

/*********************************************************************************************************************************************************
																														//Антенна с модельки радиоастрона
		for(int i = 0; i<Num_Squares-8500;i++){
//			Normal = CalcNormal(Vertexs[Squares[i].v1],Vertexs[Squares[i].v2],Vertexs[Squares[i].v3]);
//			glBegin(GL_TRIANGLE_STRIP);
//	        glNormal3fv((GLfloat*)&Squares[i].Normal);
//			glVertex3fv((GLfloat*)&Vertexs[Squares[i].v1]);
//			glEnd();
			glBegin(GL_POLYGON);
//			glBegin(GL_QUADS);
	        glNormal3f(Squares[i].Normal.x, Squares[i].Normal.y, Squares[i].Normal.z);
			glVertex3f(Vertexs[Squares[i].v1].x,Vertexs[Squares[i].v1].y,Vertexs[Squares[i].v1].z);
			glVertex3f(Vertexs[Squares[i].v2].x,Vertexs[Squares[i].v2].y,Vertexs[Squares[i].v2].z);
			glVertex3f(Vertexs[Squares[i].v3].x,Vertexs[Squares[i].v3].y,Vertexs[Squares[i].v3].z);
			glVertex3f(Vertexs[Squares[i].v4].x,Vertexs[Squares[i].v4].y,Vertexs[Squares[i].v4].z);
			glEnd(); 
		} 
		
		for(int i = 0; i<Num_Triangles;i++){
//			Normal = CalcNormal(Vertexs[Triangles[i].v1],Vertexs[Triangles[i].v2],Vertexs[Triangles[i].v3]);
//			glBegin(GL_TRIANGLE_STRIP);
//	        glNormal3fv((GLfloat*)&Triangles[i].Normal);
//			glVertex3fv((GLfloat*)&Vertexs[Triangles[i].v1]);
//			glEnd();
			glBegin(GL_TRIANGLES);
	        glNormal3f(Triangles[i].Normal.x, Triangles[i].Normal.y, Triangles[i].Normal.z);
			glVertex3f(Vertexs[Triangles[i].v1].x,Vertexs[Triangles[i].v1].y,Vertexs[Triangles[i].v1].z);
			glVertex3f(Vertexs[Triangles[i].v2].x,Vertexs[Triangles[i].v2].y,Vertexs[Triangles[i].v2].z);
			glVertex3f(Vertexs[Triangles[i].v3].x,Vertexs[Triangles[i].v3].y,Vertexs[Triangles[i].v3].z);
			glEnd();  
		}  
		
		for(int i = 0; i<Num_Lines-500;i++){
			glBegin(GL_LINES);
			glVertex3f(Vertexs[Lines[i].v1].x,Vertexs[Lines[i].v1].y,Vertexs[Lines[i].v1].z);
			glVertex3f(Vertexs[Lines[i].v2].x,Vertexs[Lines[i].v2].y,Vertexs[Lines[i].v2].z);
			glEnd();
		} 
/***********************************************************************************************************************************************************/

			glCallList(lists[9]); // надпись млечный путь
	        glEnable(GL_TEXTURE_2D);
		for(int i = 0; i<num_poligons;i++){
//			Normal = CalcNormal(Vertexs[Triangles[i].v1],Vertexs[Triangles[i].v2],Vertexs[Triangles[i].v3]);
//			glBegin(GL_TRIANGLE_STRIP);
//	        glNormal3fv((GLfloat*)&Triangles[i].Normal);
//			glVertex3fv((GLfloat*)&Vertexs[Triangles[i].v1]);
//			glEnd();
			glBegin(GL_TRIANGLES);
	        glNormal3f(normals[i].x, normals[i].y, normals[i].z);
			glTexCoord2f(textures[poligons[i].v1].u,textures[poligons[i].v1].v);
			glVertex3f(vertixs[poligons[i].v1].x,vertixs[poligons[i].v1].y,vertixs[poligons[i].v1].z);
			glTexCoord2f(textures[poligons[i].v2].u,textures[poligons[i].v2].v);
			glVertex3f(vertixs[poligons[i].v2].x,vertixs[poligons[i].v2].y,vertixs[poligons[i].v2].z);
			glTexCoord2f(textures[poligons[i].v3].u,textures[poligons[i].v3].v);
			glVertex3f(vertixs[poligons[i].v3].x,vertixs[poligons[i].v3].y,vertixs[poligons[i].v3].z);
			glEnd();  
		}  
		    glDisable(GL_TEXTURE_2D); 


		MatColor[0] = 1; MatColor[1] = 1; MatColor[2] = 1; MatColor[3] = 0;
		glLightfv(GL_LIGHT0, GL_AMBIENT, MatColor);
		glLightfv(GL_LIGHT0, GL_DIFFUSE, MatColor);
		global_ambient[0] = 2;
		global_ambient[1] = 2;
		global_ambient[2] = 2;
		global_ambient[3] = 1;
		glLightModelfv(GL_LIGHT_MODEL_AMBIENT,global_ambient);
//	glEndList();
//	}else{
//		glPushMatrix();
//		glCallList(15);
//		glPopMatrix();
//	}
/***********************************************************************************************************************/
 /*       glEnable(GL_LIGHTING);
	glColor3f(R,G,B);
			GLUquadric * quadObj;
			quadObj = gluNewQuadric();
			gluQuadricTexture(quadObj,TRUE);
			gluQuadricOrientation(quadObj,GLU_INSIDE);
//		gluQuadricOrientation(quadObj,GLU_OUTSIDE);
//			glCallList(lists[1]);
			double eqn[4];
			eqn[0] = 0; eqn[1]=-1; eqn[2]=0;eqn[3]=0;
			glClipPlane(GL_CLIP_PLANE0,eqn);
			glEnable(GL_CLIP_PLANE0);
			glTranslatef(0.000,0.011,0.00);
			gluSphere(quadObj,R_earth/4,30,30);
			glDisable(GL_CLIP_PLANE0);
			//DrawAntenna();
			gluDeleteQuadric(quadObj); */

	return 0;
}


int InitWndGraphSettings()
{
	dc = GetDC(wnd);
	SetDcPixelFormat(dc);

	dc_earth = GetDC(wnd_earth);
	SetDcPixelFormat(dc_earth);

//	dc_cap = GetDC(wnd_cap);

//	dc_graph = GetDC(wnd_graph);
//	SetDcPixelFormat(dc_graph);
//    HRC_graph = wglCreateContext(dc_graph);

    // Make a GL Context
    HRC = wglCreateContext(dc);
    wglMakeCurrent(dc, HRC);
    // Clear background color to black
    glClearColor(0.0, 0.0, 0.0, 0.0);
    // Clear the depth buffer
    glClearDepth(2.0);
    // Type of depth test
    glDepthFunc(GL_LESS);
    // Smooth color shading
    glShadeModel(GL_SMOOTH);
    // Depth test
    glEnable(GL_DEPTH_TEST);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
//    glOrtho(-10,10,-10,10,-10,10);
			gluPerspective(90,1.0,0.05,25);
//	gluPerspective(90,(double)WndWidth(hWnd)/WndHeight(hWnd),0.05,25);
	glViewport(0, 0, WndWidth(wnd), WndHeight(wnd));
    glLoadIdentity();  
    glMatrixMode(GL_MODELVIEW);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);


	CamX = 1.7427; CamY  = 0.93; CamZ = -0.44; CamAngle = -0.65; CamAngleNaklona = -0.35; CamAlpha = 0.0; CamTeta = 0.0;
	CamX_Earth = 0.055; CamY_Earth  = 0.104; CamZ_Earth = 0.0389; CamAngle_Earth = -1.09; CamAngleNaklona_Earth = 0.32; CamAlpha_Earth = 0.0; CamTeta_Earth = 0.0;
	EarthX = 0.0; EarthY = 0.00; EarthZ = 0.00; EarthAngle = 0.0; EarthAngleNaklona = 0.41;
	LunaX = 0.3; LunaY = 0.4; LunaZ = 0.2; LunaAngle = 0.0; LunaAngleNaklona = 0.1;
	SunX = -1; SunY = 1; SunZ = 0.1;
	RadioAstronX = -0.3; RadioAstronY = 0.2; RadioAstronZ = 0.5; RadioAstronAngle = 0.0; RadioAstronAngleNaklona = 0.1;
	R_earth = 0.06; R_luna = 1.25*R_earth;
	scale_km = 1.0/200000.0;
	MouseX = 0.5;
	MouseY = 0.5;
	MouseKey = 0;
	MouseX_Earth = 0.5;
	MouseY_Earth = 0.5;
	MouseKey_Earth = 0;
	read_textures();
	read_model();
	read_orbita();
	read_telescopes();
	DRAW_EKL = 1;
	orbita_current = 1;
	Time = orbita[0].time; printf("%f \n", Time);
	dTime = 1000;
	Time_old = dTime;
	//Time = 100;

	for (int i = 0; i < n; ++i)
	{
		mas[i].u = random(2*umax)-umax;
		mas[i].v = random(2*vmax)-vmax;
	};
	m_wnd.SetParameters(mas, n, umax, vmax);
	n_max = UVtime/dTime*1000;
	n = 1;


	int screenW=GetSystemMetrics(SM_CXSCREEN);//Получить ширину экрана
	int  screenH=GetSystemMetrics(SM_CYSCREEN);//Получить высоту  экрана
//	if (!m_wnd.CreateEx(NULL,NULL, "UV-Plane", WS_OVERLAPPEDWINDOW|WS_VISIBLE, CRect(800, 400, 1200, 800), NULL, NULL, NULL))
	if (!m_wnd.CreateEx(NULL,NULL, "UV-Plane", WS_OVERLAPPEDWINDOW|WS_VISIBLE, CRect(screenW*2/3, screenH/2, screenW, screenH), NULL, NULL, NULL))
		AfxMessageBox("---not");


/***************EARTH*********************************/
    HRC_earth = wglCreateContext(dc_earth );
    wglMakeCurrent(dc_earth , HRC_earth );
    // Clear background color to black
    glClearColor(0.0, 0.0, 0.0, 0.0);
    // Clear the depth buffer
    glClearDepth(2.0);
    // Type of depth test
    glDepthFunc(GL_LESS);
    // Smooth color shading
    glShadeModel(GL_SMOOTH);
    // Depth test
    glEnable(GL_DEPTH_TEST);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
//    glOrtho(-10,10,-10,10,-10,10);
			gluPerspective(90,1.0,0.05,25);
//	gluPerspective(90,(double)WndWidth(hWnd)/WndHeight(hWnd),0.05,25);
	glViewport(0, 0, WndWidth(wnd), WndHeight(wnd));
    glLoadIdentity();  
    glMatrixMode(GL_MODELVIEW);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
	read_textures();
	read_telescope_model();
/*****************************************************/
	DrawAntenna(0,0,0);
	return 0;
}

int CalculateUVPlane()
{
//	for (int i = 0; i < n; ++i)
//	{
//		mas[i].u = random(2*umax)-umax;
//		mas[i].v = random(2*vmax)-vmax;
//	};
	double alpha,x2,y2,x3,y3,x2j,y2j,x3j,y3j;
	alpha = -(Time + T0)/24.0/60.0/60.0;
	alpha = alpha - (int)(alpha);
	alpha = alpha*2.0*Pi;
	double cs = cos(alpha);
	double sn = sin(alpha);
	double H, sigma;
//	double H_source, sigma_source;
	double x_source,y_source,z_source,cos_source;


	//n = 0;
	for (int j = 0; j < NumTelescops; ++j){
	  for (int i = 0; i < NumTelescops; ++i){

		VECTOR uv;
		if (i != j){

//	for (int i = 1; i < NumTelescops; ++i){
			CoordTelescops[i].Y = -CoordTelescops[i].Y; CoordTelescops[j].Y = -CoordTelescops[j].Y;
			x2 = CoordTelescops[i].X*cs-CoordTelescops[i].Y*sn;
			y2 = CoordTelescops[i].X*sn+CoordTelescops[i].Y*cs;
//			y2 = -y2;
			

			x3 = CoordTelescops[i].X; y3 = CoordTelescops[i].Y;
			CoordTelescops[i].X = x2; CoordTelescops[i].Y = y2;
			x2j = CoordTelescops[j].X*cs-CoordTelescops[j].Y*sn;
			y2j = CoordTelescops[j].X*sn+CoordTelescops[j].Y*cs;
//			y2j = -y2j;
			x3j = CoordTelescops[j].X; y3j = CoordTelescops[j].Y;
			CoordTelescops[j].X = x2j; CoordTelescops[j].Y = y2j;
//	}
			CoordTelescops[0].X =  RadioAstronX / scale_km*1000.0;	
			CoordTelescops[0].Y =  RadioAstronY / scale_km*1000.0;
			CoordTelescops[0].Z =  RadioAstronZ / scale_km*1000.0;

			if (abs(y2)> 0.000000001){H = acos(x2/sqrt(x2*x2+y2*y2))*y2/abs(y2);}else{H = 0;}
			if ( (x2*x2+y2*y2+CoordTelescops[i].Z*CoordTelescops[i].Z) > 0.000001){
				sigma = asin(CoordTelescops[i].Z/sqrt(x2*x2+y2*y2+CoordTelescops[i].Z*CoordTelescops[i].Z));
			}else{
				sigma = 0;
				//printf("***********************888\n");
			}
///			H_source = 1.5;
//			sigma_source = 0.7;
			double delH,delS;
//			delH = abs(H_source + H); if (delH > 2*Pi){delH = delH - 2*Pi;}if (delH < -2*Pi){delH = delH + 2*Pi;}
//			delS = abs(sigma_source - sigma); if (delS > 2*Pi){delS = delS - 2*Pi;}if (delS < -2*Pi){delS = delS + 2*Pi;}
			delH = (H_source - H ); delH = delH - 2.0*Pi*(int)(delH/Pi/2.0);  if (delH > Pi){delH = delH - 2*Pi;}if (delH < -Pi){delH = delH + 2*Pi;}
//			delH = (H_source - H ); if (delH > 2*Pi){delH = delH - 2*Pi;}if (delH < -2*Pi){delH = delH + 2*Pi;}
			delS = (sigma_source - sigma); delS = delS - 2.0*Pi*(int)(delS/Pi/2.0); if (delS > 2*Pi){delS = delS - 2*Pi;}if (delS < -2*Pi){delS = delS + 2*Pi;}
//			printf("%f %f %f %f\n",H_source,H,delH,delS);
//			printf("%f %f %f %f\n",sigma_source,sigma,delH,delS);


			if ((abs(delH) < Pi/2)&&(abs(delS) < Pi/2)){
				CoordTelescops[i].flag = 1;
			}else{
				CoordTelescops[i].flag = 0;
			}
//				CoordTelescops[i].flag = 1;
			x_source = cos(H_source)*cos(sigma_source);			
			y_source = sin(H_source)*cos(sigma_source);			
			z_source = sin(sigma_source);			

/*			double x_ant_u,y_ant_u,z_ant_u,x_ant,y_ant,z_ant,x_ant_v,y_ant_v,z_ant_v,cos_u,cos_v;
			x_ant_u = -sin(H);			
			y_ant_u = cos(H);			
			z_ant_u = 0;			
			x_ant_v = -cos(H)*sin(sigma);			
			y_ant_v = -sin(H)*sin(sigma);			
			z_ant_v = cos(sigma);			
			cos_u = (x_ant_u*x_source+y_ant_u*y_source+z_ant_u*z_source)/
				sqrt(x_ant_u*x_ant_u+y_ant_u*y_ant_u+z_ant_u*z_ant_u);
			cos_v = (x_ant_v*x_source+y_ant_v*y_source+z_ant_v*z_source)/
				sqrt(x_ant_v*x_ant_v+y_ant_v*y_ant_v+z_ant_v*z_ant_v); */
			cos_source = (CoordTelescops[i].X*x_source+CoordTelescops[i].Y*y_source+CoordTelescops[i].Z*z_source)/
				sqrt(CoordTelescops[i].X*CoordTelescops[i].X+CoordTelescops[i].Y*CoordTelescops[i].Y+CoordTelescops[i].Z*CoordTelescops[i].Z);


/*                                                                                    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

			CoordTelescops[i].flag = 1;
		if (i != 0){
//			if (acos(cos_source) > 0.8*Pi/2.0){
//			if (acos(cos_source) > Pi/2.0){
			if (acos(cos_source) > Pi/2.0){
				CoordTelescops[i].flag = 0;
			}//else{
//			if (acos(cos_source) < CoordTelescops[i].Zmin*Pi/180.0){
//			//	CoordTelescops[i].flag = 0;
//			}
		}
*/
                                                                                    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

			CoordTelescops[i].flag = 1;
		if (i != 0){
//			if (acos(cos_source) > 0.8*Pi/2.0){
//			if (acos(cos_source) > Pi/2.0){
			if (acos(cos_source) > CoordTelescops[i].Zmax*Pi/180.0){
				CoordTelescops[i].flag = 0;
			}//else{
			if (acos(cos_source) < CoordTelescops[i].Zmin*Pi/180.0){
			//	CoordTelescops[i].flag = 0;
			}
		}


			//	CoordTelescops[i].flag = 1;
/*			double koeff_A,koeff_B,koeff_C; // плоскость
			koeff_A = 0; koeff_B = 0; koeff_C = 1;
			CoordTelescops[i].flag = 0;
			if ((x_source*koeff_A+y_source*koeff_B+z_source*koeff_C) > 0 ){
				CoordTelescops[i].flag = 0;
			} */
//				CoordTelescops[i].flag = 0;
		//	}
//			if (i == 2){
//				printf("c %f \n",acos(cos_source)/Pi*180);}

//			printf("%f %f %f %f %f\n",H_source,H,delH,delS,acos(cos_source));

			//if (i == 3){
			//printf("%i %f %f %f %f\n",i,CoordTelescops[i].Zmin,CoordTelescops[i].Zmax,CoordTelescops[i].AZmin,CoordTelescops[i].AZmax);
//				printf("%f %f %f %f %f\n",H_source/Pi,H/Pi,delH/Pi,delS/Pi,acos(cos_source)/Pi);
			///	printf("%f %f %f %f %f\n",acos(cos_u)/Pi,acos(cos_v)/Pi,delH/Pi,delS/Pi,acos(cos_source)/Pi);
//				CoordTelescops[i].flag = 1;
			//}else{
//				CoordTelescops[i].flag = 0;
			//}
//			if (( acos(cos_u)/Pi < 0.5)&&(acos(cos_v)/Pi < 0.5)){
//				CoordTelescops[i].flag = 1;
//			}else{
//				CoordTelescops[i].flag = 0;
//			}


/***			double koeff_a,koeff_b,koeff_c; // конус
			double koeff_A,koeff_B,koeff_C; // плоскость
			koeff_a = 1;
			koeff_b = 1;
			koeff_c = 0.001;
//			x_source = cos(delH)*cos(delS);			
//			y_source = sin(delH)*cos(delS);			
//			z_source = sin(delS);			

			double x_source_t, y_source_t, z_source_t;
	//		x_source = CoordTelescops[i].X;
	//		y_source = CoordTelescops[i].Y;
	//		z_source = CoordTelescops[i].Z;


			x_source_t = x_source*cos(-H)  -  y_source*sin(-H);
			y_source_t = x_source*sin(-H)  +  y_source*cos(-H);
			x_source = x_source_t; y_source = y_source_t;

			x_source_t = x_source*cos(Pi/2-sigma)  -  z_source*sin(Pi/2-sigma);
			z_source_t = x_source*sin(Pi/2-sigma)  +  z_source*cos(Pi/2-sigma);

			x_source = x_source_t; z_source = z_source_t;

			if (i == 3){
			printf("%f %f %f \n",x_source, y_source, z_source);
			}
			koeff_A = 0; koeff_B = 0; koeff_C = 1;
			CoordTelescops[i].flag = 0;
			if ( (sqr(x_source/koeff_a)+sqr(y_source/koeff_b)-sqr(z_source/koeff_c) ) < 0){
				if ((x_source*koeff_A+y_source*koeff_B+z_source*koeff_C) > 0 ){
					CoordTelescops[i].flag = 1;
				}
			}
//			}else{
//				CoordTelescops[i].flag = 0;
//			}
*/

			uv = CalculateUVPoint(j, i, H_source, sigma_source, 400000000);
			CoordTelescops[i].X = x3; CoordTelescops[i].Y = y3;
			CoordTelescops[j].X = x3j; CoordTelescops[j].Y = y3j;
			CoordTelescops[i].Y = -CoordTelescops[i].Y; CoordTelescops[j].Y = -CoordTelescops[j].Y;

			//CoordTelescops[i].flag = i - i / 2 * 2; 
			//printf("%i\n",CoordTelescops[i].flag);
			
//UVtime
//			continue;

		if ( CoordTelescops[i].flag == 1)
		{
			if (n == 1){
				for (int k = 1; k <= n_max; ++k){
					mas[k].u = 0;
					mas[k].v = 0;
					mas[k].t = 0;
				}		
				n++;
			}
			if (mas[1].t > Time-UVtime){
				mas[n].u = uv.x;
				mas[n].v = uv.y;
				mas[n].t = Time;
				n++;
				//printf("!!!!!!!! %i !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n",n);
				
			}
			else
			{
				//printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n");
				for (int k = 1; k <= n_max; ++k){
					mas[k-1].u = mas[k].u;
					mas[k-1].v = mas[k].v;
					mas[k-1].t = mas[k].t;
				}
				mas[n].u = uv.x;
				mas[n].v = uv.y;
				mas[n].t = Time;
			}
		}

/*		if ( CoordTelescops[i].flag == 1)
		{
			if (n < n_max){
				mas[n].u = uv.x;
				mas[n].v = uv.y;
				mas[n].t = Time;
				n++;
				//printf("!!!!!!!! %i !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n",n);
				
			}
			else
			{
				//printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n");
				for (int k = 1; k <= n_max; ++k){
					mas[k-1].u = mas[k].u;
					mas[k-1].v = mas[k].v;
					mas[k-1].t = mas[k].t;
				}
				mas[n].u = uv.x;
				mas[n].v = uv.y;
				mas[n].t = Time;
			}
		}

/***/

		}

	  }
	}

/*	printf("*******************************************begin\n");
	for (int i = 1; i < NumTelescops; i++){

		if (CoordTelescops[i].flag == 1){
			printf("** N %i %s %s %f- %f %f - %f ***\n", i, CoordTelescops[i].Name2,CoordTelescops[i].Name,CoordTelescops[i].AZmin,CoordTelescops[i].AZmax,CoordTelescops[i].Zmin,CoordTelescops[i].Zmax);
		}
	}
	printf("*******************************************end\n"); */
	m_wnd.SetParameters(mas, n, umax, vmax);

	return 0;
}

int CalculatePositions()
{
	EarthAngle = - 2.0*Pi*(Time+T0)/(24.0*60.0*60.0);
	EarthAngle = EarthAngle - ( (int)(EarthAngle / (2.0*Pi) ) * (2.0*Pi) );
//	printf("&&& %f",EarthAngle);
	double Lf;
	double Moon_F0 = -12*24.0*60.0*60.0;
	Lf =  2.0*Pi*(Time + Moon_F0)/(27.321582*24.0*60.0*60.0);
	Lf = Lf - ( (int)(Lf / (2.0*Pi) ) * (2.0*Pi) );
	LunaX = EarthX + 384000*scale_km*cos(Lf);
	LunaY = EarthY + 384000*scale_km*sin(Lf);
	LunaZ = EarthZ ;//+ 384000*scale_km*sin(Time);
	double Sf;
	double Sun_F0 = -27.0*24.0*60.0*60.0;
	Sf = 2.0*Pi*(Time + Sun_F0)/(365.2422*24.0*60.0*60.0);
	Sf = Sf - ( (int)(Sf / (2.0*Pi) ) * (2.0*Pi) );
	SunX = EarthX - 1.7*cos(Sf);
	SunY = EarthY + 1.7*sin(Sf);
	SunZ = EarthZ;// + 10.7*sin(5*Time);
	LunaAngle =  Lf +Pi/2+Pi;
	for (int i = 0; i < orbita_last;i++){
		if (Time < orbita[i].time){
			orbita_current = i;
//			printf("%f %f \n", Time,orbita[0].time);
			break;
		}
	}
///	if(Time > orbita[orbita_last].time){
//		orbita_current = orbita_last;
//	}
				//orbita_current++;
	double dx,dy,dz;
	dx = (orbita[orbita_current].r.x-orbita[orbita_current-1].r.x)/(orbita[orbita_current].time-orbita[orbita_current-1].time)*(Time-orbita[orbita_current-1].time);
	dy = (orbita[orbita_current].r.y-orbita[orbita_current-1].r.y)/(orbita[orbita_current].time-orbita[orbita_current-1].time)*(Time-orbita[orbita_current-1].time);
	dz = (orbita[orbita_current].r.z-orbita[orbita_current-1].r.z)/(orbita[orbita_current].time-orbita[orbita_current-1].time)*(Time-orbita[orbita_current-1].time);
	RadioAstronX = (orbita[orbita_current-1].r.x+dx)*scale_km;
	RadioAstronY = (orbita[orbita_current-1].r.y+dy)*scale_km;
	RadioAstronZ = (orbita[orbita_current-1].r.z+dz)*scale_km;
	RadioAstronAngleNaklona = sigma_source;
	RadioAstronAngle = H_source-Pi/2;

	VECTOR uv;
	CoordTelescops[0].X =  RadioAstronX / scale_km*1000.0;
	CoordTelescops[0].Y =  RadioAstronY / scale_km*1000.0;
	CoordTelescops[0].Z =  RadioAstronZ / scale_km*1000.0;

	
	uv = CalculateUVPoint(0, 1, 0.3, 0.2, 400000000);
//	printf("U %f V %f \n",uv.x,uv.y);
	CalculateUVPlane();

//	if ( (Time*100) > 11000){Time = 10;}
	if ( Time > orbita[orbita_last-100].time){Time = orbita[0].time; orbita_current = 1;}
	//if ( orbita_current >= orbita_last-100){Time = orbita[0].time; orbita_current = 1;}
	return 0;
}

int wm_paint_func(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam, HDC dc, HGLRC HRC, double CamX, double CamY, double CamZ, double CamAngleNaklona, double CamAngle,double MouseX, double MouseY, int MouseKey, int DRAW_EKL)
{
PAINTSTRUCT Ps;
        BeginPaint (hWnd, &Ps);
        wglMakeCurrent(dc, HRC);
		glViewport(0, 0, WndWidth(hWnd), WndHeight(hWnd));
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);                   // Clear color & depth buffers
		char sTitle[256];// = (char*)malloc(256);
		sTitle[0] = '\0';
		char sTime[20];

        NewCount = GetTickCount();
        FrameCount += 1;
		if ((NewCount - LastCount) > 1000){
	        itoa(FrameCount * 1000 / (NewCount - LastCount),szString,10);
			//szString = szString + "123fs";
//			strcat(szString,"asd");sgf
		    LastCount = NewCount;
			FrameCount = 0;
		}
		int tttt_day, tttt_month, tttt_year;
		double tttt_t;
		if(hWnd==wnd){
			strcat(sTitle,"Fps: ");
			strcat(sTitle,szString);
			strcat(sTitle," Radioastron orbita");
			strcat(sTitle," Date: ");
			strcat(sTitle,JulDateToDate(Time/24.0/60.0/60.0, tttt_day, tttt_month, tttt_year));
			strcat(sTitle," Time: ");
			tttt_t = (Time/24.0/60.0/60.0-0.5 -(int)(Time/24.0/60.0/60.0 -0.5))*24;
			itoa(tttt_t,sTime,10);
			strcat(sTitle,sTime);
			strcat(sTitle,":");
			tttt_t = (tttt_t - (int)(tttt_t))*60;
			itoa(tttt_t,sTime,10);
			strcat(sTitle,sTime);
			strcat(sTitle,":");
			tttt_t = (tttt_t - (int)(tttt_t))*60;
			itoa(tttt_t,sTime,10);
			strcat(sTitle,sTime);
		}else{
//			strcat(sTitle,"Fps: ");
//			strcat(sTitle,szString);
			strcat(sTitle,"Telescopes: ");
			for (int i = 1; i < NumTelescops; i++){

				if (CoordTelescops[i].flag == 1){
					strcat(sTitle,CoordTelescops[i].Name2);
					strcat(sTitle," ");
//					printf("** N %i %s %s %f- %f %f - %f ***\n", i, CoordTelescops[i].Name2,CoordTelescops[i].Name,CoordTelescops[i].AZmin,CoordTelescops[i].AZmax,CoordTelescops[i].Zmin,CoordTelescops[i].Zmax);
				}
			}
		}
		SendMessage(hWnd, WM_SETTEXT, 0, LPARAM(&sTitle));
//       CString tttt_s;// = "rewqrwerfwerfwewerwerfwe";
//		double jdate;
//		tttt_day = 25;
//		tttt_month = 3;
//		tttt_year = 2000;
//        printf("*");
//        tttt_s = JulDateToDate(Time/24.0/60.0/60.0, tttt_day, tttt_month, tttt_year);
  //      printf("\n %i %i %i \n",tttt_day, tttt_month, tttt_year);

//		printf(tttt_s);
		GLUquadric * quadObj;
		
        glLoadIdentity();
		glRotated(-90,0.0,0.0,1.0);
		glRotated(90,0.0,1.0,0.0);
		glRotated(CamAngleNaklona*180/Pi,0.0,1.0,0.0);
		glRotated(CamAngle*180/Pi,0.0,0.0,1.0);
		glRotated(EarthAngleNaklona*180/Pi,0.0,1.0,0.0);

		//glScalef(1,1,-1);
		//glScalef(1,1,-1);
		double Alpha_gal;
		Alpha_gal = -280+90;
		glRotated(Alpha_gal,0.0,0.0,1.0);  //поворот так чтобы ось Х была направлена на точку весеннего равноденствия+280
		glRotated(57,1.0,0.0,0.0);


        glDisable(GL_LIGHTING);



		glDisable(GL_DEPTH_TEST);
		glEnable(GL_LIGHT_MODEL_TWO_SIDE);//gl_light_model_two_side
//		glCullFace(GL_BACK);


//		glCullFace(GL_FRONT);

		glCullFace(GL_BACK);

		quadObj = gluNewQuadric();
		gluQuadricTexture(quadObj,TRUE);
		gluQuadricOrientation(quadObj,GLU_OUTSIDE);
        glCallList(lists[0]);
        glEnable(GL_TEXTURE_2D);

				gluSphere(quadObj,14.5,250,250);                                                  /******/

        glDisable(GL_TEXTURE_2D);
		gluDeleteQuadric(quadObj);

		glLineStipple(2,0xF0F0);  //galaxy
		if(DRAW_EKL == 1){
			glColor3f(0,0,1);
	        glBegin(GL_LINES);
		    glVertex3f(-15,-15,0);
		    glVertex3f(15,-15,0);
		    glVertex3f(15,-15,0);
		    glVertex3f(15,15,0);
		    glVertex3f(15,15,0);
		    glVertex3f(-15,15,0);
		    glVertex3f(-15,15,0);
		    glVertex3f(-15,-15,0);
	        glEnd();
			glColor3f(1,1,1);
		}
		glLineStipple(2,0xFFFF); 

		if(DRAW_EKL == 1){		
			glScalef(1,1,-1);
			glRotatef(160,0,0,1);
			glTranslatef(14,0,0.5);

			glCallList(lists[6]); // надпись млечный путь
	        glEnable(GL_TEXTURE_2D);
		    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
			glEnable(GL_BLEND);
	        glBegin(GL_POLYGON);
		    glNormal3f(0.0, 0.0, -1.0);
			glTexCoord2f(0,0); glVertex3f(-1.0,2.0, 0.2);
	        glTexCoord2f(1,0); glVertex3f(-1.0,2.0, -0.2);
		    glTexCoord2f(1,-1); glVertex3f(-1.0,-2.0, -0.2);
			glTexCoord2f(0,-1); glVertex3f(-1.0,-2.0, 0.2);
	        glEnd();
			glDisable(GL_BLEND);
		    glDisable(GL_TEXTURE_2D); 
			glTranslatef(-14,0,-0.5);
			glRotatef(-160,0,0,1);
			glScalef(1,1,-1);
		}
		glRotated(-57,1.0,0.0,0.0);
		glRotated(-Alpha_gal,0.0,0.0,1.0);  //поворот так чтобы ось Х была направлена на точку весеннего равноденствия

		if(DRAW_EKL == 1){		
		glScalef(1,1,-1);
			glTranslatef(14,0,0.5);

		    glCallList(lists[4]); // надпись эклиптика
			glEnable(GL_TEXTURE_2D);
	        glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
			glEnable(GL_BLEND);
			glBegin(GL_POLYGON);
	        glNormal3f(0.0, 0.0, -1.0);
		    glTexCoord2f(0,0); glVertex3f(-1.0,2.0, 0.2);
			glTexCoord2f(1,0); glVertex3f(-1.0,2.0, -0.2);
	        glTexCoord2f(1,-1); glVertex3f(-1.0,-2.0, -0.2);
		    glTexCoord2f(0,-1); glVertex3f(-1.0,-2.0, 0.2);
			glEnd();
			glDisable(GL_BLEND);
		    glDisable(GL_TEXTURE_2D); 
		}
        glLoadIdentity();
        glLoadIdentity();

		glRotated(-90,0.0,0.0,1.0);
		glRotated(90,0.0,1.0,0.0);
		glRotated(CamAngleNaklona*180/Pi,0.0,1.0,0.0);
		glRotated(CamAngle*180/Pi,0.0,0.0,1.0);
		glRotated(EarthAngleNaklona*180/Pi,0.0,1.0,0.0);
		glScalef(1,1,-1);

		if(DRAW_EKL == 1){					//эклиптика
		    glDisable(GL_LIGHTING);
   			glDisable(GL_DEPTH_TEST);
			glColor3f(1,0,0);
			glLineStipple(2,0xF0F0);
			glEnable(GL_LINE_STIPPLE);
	        glBegin(GL_LINE_LOOP);
		    glVertex3f(-15,-15,0);
		    glVertex3f(15,-15,0);
		    glVertex3f(15,15,0);
		    glVertex3f(-15,15,0);
	        glEnd();
			glColor3f(1,1,1);
		} 
        glLoadIdentity();
		glRotated(-90,0.0,0.0,1.0);
		glRotated(90,0.0,1.0,0.0);
		glRotated(CamAngleNaklona*180/Pi,0.0,1.0,0.0);
		glRotated(CamAngle*180/Pi,0.0,0.0,1.0);
		glRotated(EarthAngleNaklona*180/Pi,0.0,1.0,0.0);

 		glTranslatef(CamX,CamY,CamZ);
		glScalef(1,1,-1);
		glRotated(EarthAngleNaklona*180/Pi,0.0,1.0,0.0);
		glRotated(90,0.0,0.0,1.0);  //поворот так чтобы ось Х была направлена на точку весеннего равноденствия

        glLoadIdentity();
		glRotated(-90,0.0,0.0,1.0);
		glRotated(90,0.0,1.0,0.0);
		glRotated(CamAngleNaklona*180/Pi,0.0,1.0,0.0);
		glRotated(CamAngle*180/Pi,0.0,0.0,1.0);
		glRotated(EarthAngleNaklona*180/Pi,0.0,1.0,0.0);

		glScalef(1,1,-1);
		glRotated(EarthAngleNaklona*180/Pi,0.0,1.0,0.0);
		glRotated(90,0.0,0.0,1.0);  //поворот так чтобы ось Х была направлена на точку весеннего равноденствия


		if(DRAW_EKL == 1){					//небесный экватор
			glColor3f(0,1,1);
	        glBegin(GL_LINES);
		    glVertex3f(-15,-15,0);
		    glVertex3f(15,-15,0);
		    glVertex3f(15,-15,0);
		    glVertex3f(15,15,0);
		    glVertex3f(15,15,0);
		    glVertex3f(-15,15,0);
		    glVertex3f(-15,15,0);
		    glVertex3f(-15,-15,0);
	        glEnd();
			glColor3f(1,1,1);
		}
		glLineStipple(2,0xFFFF);

		if(DRAW_EKL == 1){		
			glRotatef(-40,0,0,1);
			glTranslatef(14,0,0.5);
			glCallList(lists[5]); // надпись небесный экватор
			glEnable(GL_TEXTURE_2D);
	        glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
			glEnable(GL_BLEND);
			glBegin(GL_POLYGON);
	        glNormal3f(0.0, 0.0, -1.0);
		    glTexCoord2f(0,0); glVertex3f(-1.0,2.0, 0.2);
			glTexCoord2f(1,0); glVertex3f(-1.0,2.0, -0.2);
	        glTexCoord2f(1,-1); glVertex3f(-1.0,-2.0, -0.2);
		    glTexCoord2f(0,-1); glVertex3f(-1.0,-2.0, 0.2);
			glEnd();
			glDisable(GL_BLEND);
		    glDisable(GL_TEXTURE_2D); 
			glTranslatef(-14,0,-0.5);
			glRotatef(40,0,0,1);
		}


		if(DRAW_EKL == 1){		
			glTranslatef(14,0,0.5);
		    glCallList(lists[7]); // надпись гамма
	       glEnable(GL_TEXTURE_2D);
			glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
			glEnable(GL_BLEND);
	        glBegin(GL_POLYGON);
			glNormal3f(0.0, 0.0, -1.0);
		    glTexCoord2f(0,0); glVertex3f(-1.0,0.2, 0.3);
	        glTexCoord2f(1,0); glVertex3f(-1.0,0.2, -0.3);
			glTexCoord2f(1,-1); glVertex3f(-1.0,-0.2, -0.3);
		    glTexCoord2f(0,-1); glVertex3f(-1.0,-0.2, 0.3);
	        glEnd();
			glDisable(GL_BLEND);
		    glDisable(GL_TEXTURE_2D); 
			glTranslatef(-14,0,-0.5);
		}
	   glColor3f(1,1,0);
		glColor3f(1,1,1);


		
		
		double SunRnd;

        glLoadIdentity();     // Солнце 
		glRotated(-90,0.0,0.0,1.0);
		glRotated(90,0.0,1.0,0.0);
		glRotated(CamAngleNaklona*180/Pi,0.0,1.0,0.0);
		glRotated(CamAngle*180/Pi,0.0,0.0,1.0);
		glRotated(EarthAngleNaklona*180/Pi,0.0,1.0,0.0);

		glScalef(1,1,-1);
		glScalef(1,1,-1);
		Alpha_gal = -280+90;

		global_ambient[0] = 2;
		global_ambient[1] = 2;
		global_ambient[2] = 2;
		global_ambient[3] = 1;
		glLightModelfv(GL_LIGHT_MODEL_AMBIENT,global_ambient);

		LightPos[0] = -SunX*1000;
		LightPos[1] = -SunY*1000;
		LightPos[2] = -SunZ*1000;
		LightPos[3] = 0.0;
		glLightfv(GL_LIGHT0,GL_POSITION,LightPos);
		glColor3f(1,1,0);

		glTranslatef(SunX,SunY,SunZ);
		glDisable(GL_CULL_FACE);
		glCullFace(GL_FRONT);
		quadObj = gluNewQuadric();
		gluQuadricOrientation(quadObj,GLU_OUTSIDE);
		gluDeleteQuadric(quadObj);
		SunRnd = random(0.1);
		glColor3f(SunRnd+0.9,SunRnd+0.9,0);
			glCallList(lists[3]); // солнце
	        glEnable(GL_TEXTURE_2D);
		    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
			glEnable(GL_BLEND);
	        glBegin(GL_POLYGON);
			double CosSunH, SinSunH, SinSunSigma, CosSunSigma;
			SinSunH = -SunY / sqrt(SunX*SunX+SunY*SunY);
			CosSunH = SunX / sqrt(SunX*SunX+SunY*SunY);
			SinSunSigma = SunZ / sqrt(SunX*SunX+SunY*SunY+SunZ*SunZ);
			CosSunSigma = sqrt(SunX*SunX+SunY*SunY) / sqrt(SunX*SunX+SunY*SunY+SunZ*SunZ);
			glTexCoord2f(0,0); glVertex3f(0.2*SinSunH + 0.2*(-CosSunH*SinSunSigma),0.2*CosSunH + 0.2*(SinSunH*SinSunSigma),0.2*CosSunSigma);
	        glTexCoord2f(1,0); glVertex3f(0.2*SinSunH - 0.2*(-CosSunH*SinSunSigma),0.2*CosSunH + -0.2*(SinSunH*SinSunSigma),-0.2*CosSunSigma);
		    glTexCoord2f(1,-1); glVertex3f(-0.2*SinSunH - 0.2*(-CosSunH*SinSunSigma),-0.2*CosSunH + -0.2*(SinSunH*SinSunSigma),-0.2*CosSunSigma);
			glTexCoord2f(0,-1); glVertex3f(-0.2*SinSunH + 0.2*(-CosSunH*SinSunSigma),-0.2*CosSunH + 0.2*(SinSunH*SinSunSigma),0.2*CosSunSigma);
	        glEnd();
	        glEnd();
			glDisable(GL_BLEND);
		    glDisable(GL_TEXTURE_2D); 
		glColor3f(1,1,1);

		glTranslatef(-SunX,-SunY,-SunZ);
		
		
		
		glLoadIdentity();
		glRotated(-90,0.0,0.0,1.0);
		glRotated(90,0.0,1.0,0.0);
		glRotated(CamAngleNaklona*180/Pi,0.0,1.0,0.0);
		glRotated(CamAngle*180/Pi,0.0,0.0,1.0);







		glTranslatef(CamX,CamY,CamZ);
		glScalef(1,1,-1);

       glEnable(GL_LIGHTING);
		glCullFace(GL_FRONT);


        glDisable(GL_LIGHTING);
        glEnable(GL_LIGHTING);
		glEnable(GL_DEPTH_TEST);

        glDisable(GL_LIGHTING);
   		glEnable(GL_DEPTH_TEST);

		glTranslatef(EarthX,EarthY,EarthZ);
//		glRotated(EarthAngleNaklona*180/Pi,0.0,1.0,0.0);
		glRotated(90,0.0,0.0,1.0);  //поворот так чтобы ось Х была направлена на точку весеннего равноденствия

		glColor3f(0.70,0.70,0.0);
		scale_km = 1.0/100000.0;
		int startDi;// = 0;
		startDi = 0;
		if ( orbita_current > 300){
			startDi = orbita_current - 300;
		}
		for (int i = startDi;i< orbita_current-1; i++){			//орбита
	        glBegin(GL_LINES);
			glVertex3f(orbita[i].r.x*scale_km,orbita[i].r.y*scale_km,orbita[i].r.z*scale_km);
			glVertex3f(orbita[i+1].r.x*scale_km,orbita[i+1].r.y*scale_km,orbita[i+1].r.z*scale_km);
	        glEnd();
		}
        glBegin(GL_LINES);
		glVertex3f(orbita[orbita_current-1].r.x*scale_km,orbita[orbita_current-1].r.y*scale_km,orbita[orbita_current-1].r.z*scale_km);
		glVertex3f(RadioAstronX,RadioAstronY,RadioAstronZ);
        glEnd();

		glTranslatef(15*cos(H_source)*cos(sigma_source),15*sin(H_source)*cos(sigma_source),15*sin(sigma_source));

			glCallList(lists[3]); // источник
//			glCallList(lists[8]); // источник
	        glEnable(GL_TEXTURE_2D);
		    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
			glEnable(GL_BLEND);
	        glBegin(GL_POLYGON);
		    glNormal3f(0.0, 0.0, -1.0);

			SunRnd = random(0.4);
			glColor3f(SunRnd+0.6,SunRnd+0.6,SunRnd+0.6);
	//		glColor3f(1,1,1);
			H_source = - H_source;
			glTexCoord2f(0,0); glVertex3f(0.3*sin(H_source) + 0.3*(-cos(H_source)*sin(sigma_source)),0.3*cos(H_source) + 0.3*(sin(H_source)*sin(sigma_source)),0.3*cos(sigma_source));
	        glTexCoord2f(1,0); glVertex3f(0.3*sin(H_source) - 0.3*(-cos(H_source)*sin(sigma_source)),0.3*cos(H_source) + -0.3*(sin(H_source)*sin(sigma_source)),-0.3*cos(sigma_source));
		    glTexCoord2f(1,-1); glVertex3f(-0.3*sin(H_source) - 0.3*(-cos(H_source)*sin(sigma_source)),-0.3*cos(H_source) + -0.3*(sin(H_source)*sin(sigma_source)),-0.3*cos(sigma_source));
			glTexCoord2f(0,-1); glVertex3f(-0.3*sin(H_source) + 0.3*(-cos(H_source)*sin(sigma_source)),-0.3*cos(H_source) + 0.3*(sin(H_source)*sin(sigma_source)),0.3*cos(sigma_source));
	        glEnd();
			H_source = - H_source;
			glDisable(GL_BLEND);
		    glDisable(GL_TEXTURE_2D); 
			glTranslatef(-15*cos(H_source)*cos(sigma_source),-15*sin(H_source)*cos(sigma_source),-15*sin(sigma_source));
		//}

		double s_x,s_y,s_z;
		s_x = cos(H_source)*cos(sigma_source);
		s_y = sin(H_source)*cos(sigma_source);
		s_z = sin(sigma_source);			glColor3f(1,1,1);			//Source
		glLineStipple(1,0x9999);
		glColor3f(0,1,0);
		glBegin(GL_LINES);
		glVertex3f(0,0,0);
		glVertex3f(15*s_x,15*s_y,15*s_z);
		glVertex3f(RadioAstronX,RadioAstronY,RadioAstronZ);
		glVertex3f(15*s_x,15*s_y,15*s_z);
		glEnd();
		glColor3f(1,1,1);
		glLineStipple(2,0xFFFF);

		glColor3f(1,1,1);
    		glEnable(GL_DEPTH_TEST);



		glEnable(GL_CULL_FACE);
//		glDisable(GL_LIGHT_MODEL_TWO_SIDE);//gl_light_model_two_side//		
//		glCullFace(GL_FRONT);

//		glCullFace(GL_BACK);
    	glDisable(GL_DEPTH_TEST);
		glEnable(GL_BLEND);
//		glColor4f(1.0,1.0,1,0.10);
		glColor4f(0.5,0.5,1,0.25);
		quadObj = gluNewQuadric();
		gluQuadricOrientation(quadObj,GLU_INSIDE);
		gluSphere(quadObj,R_earth+0.0015,150,150);  // Отрисовка Атмосферы
		gluDeleteQuadric(quadObj);
    	glEnable(GL_DEPTH_TEST);
		glDisable(GL_BLEND);
		glColor4f(1,1,1,1);
		glDisable(GL_CULL_FACE);
/**/
//		glRotated(EarthAngleNaklona*180/Pi,0.0,1.0,0.0);
//		glRotated(90,0.0,0.0,1.0);
		glRotated(EarthAngle*180/Pi,0.0,0.0,1.0);
//		glRotated(-150,0.0,0.0,1.0);  //подгон текстуры земли
//		glRotated(-77,0.0,0.0,1.0);  //подгон текстуры земли
//		glRotated(-65,0.0,0.0,1.0);  //подгон текстуры земли
		glRotated(90,0.0,0.0,1.0);  //подгон текстуры земли
		quadObj = gluNewQuadric();
		gluQuadricTexture(quadObj,TRUE);
		gluQuadricOrientation(quadObj,GLU_INSIDE);
        glCallList(lists[1]);
        glEnable(GL_TEXTURE_2D);
		gluSphere(quadObj,R_earth,150,150);  // Отрисовка Земли
        glDisable(GL_TEXTURE_2D);
		gluDeleteQuadric(quadObj);
		glRotated(-90,0.0,0.0,1.0); //подгон текстуры земли

		double temp_R,temp_x,temp_y,temp_z, temp_sigma, temp_H;
		for (int i = 1; i < NumTelescops; i++){
//		  if (i == 6) //arecibo
//		  if (i == 7) //bd
//		  if ((i == 10)||(i == 7))
		  {
			temp_R = sqrt(CoordTelescops[i].X*CoordTelescops[i].X+CoordTelescops[i].Y*CoordTelescops[i].Y+CoordTelescops[i].Z*CoordTelescops[i].Z);
			temp_x = CoordTelescops[i].X/temp_R;
			temp_y = -CoordTelescops[i].Y/temp_R;
			temp_z = CoordTelescops[i].Z/temp_R;


/*			temp_z = temp_z/sqrt(temp_x*temp_x+temp_y*temp_y)/5;
			temp_R = sqrt(temp_x*temp_x+temp_y*temp_y+temp_z*temp_z);
			temp_x = temp_x/temp_R;
			temp_y = temp_y/temp_R;
			temp_z = temp_z/temp_R;*/

			//if (abs(temp_y) > 0.00001){
//			glTranslatef(1.03*R_earth*temp_x,1.03*R_earth*temp_y,1.03*R_earth*temp_z);
				temp_H = acos(temp_x)*temp_y/abs(temp_y);
			//}else{
			//	temp_H = acos(temp_x);
			//}
			temp_sigma = asin(temp_z);//+Pi/2;
			//temp_sigma = tan(temp_sigma);//random(2);
//			temp_H = random(6.28);
			//temp_sigma = 0;//random(2);
			temp_H =  temp_H-Pi/2.0;

//			temp_H = H_source;
//			temp_sigma = sigma_source;
//			printf("%f\n",sqrt(temp_x*temp_x+temp_y*temp_y+temp_z*temp_z)/R_earth);
			glPushMatrix();
//			glRotated(90,0.0,0.0,1.0);

			glTranslatef(1.01*R_earth*temp_x,1.01*R_earth*temp_y,1.01*R_earth*temp_z);   //+++++
//			glTranslatef(1.03*R_earth*temp_x,1.03*R_earth*temp_y,1.03*R_earth*temp_z);
			glRotated(temp_sigma*180/Pi,cos(temp_H),sin(temp_H),0);
			glRotated(temp_H*180/Pi,0.0,0.0,1.0);
			glRotated(-90,1.0,0.0,0.0);   //+++++++

			glScalef(0.2,0.2,0.2);  //++++
//			glScalef(0.4,0.4,0.4);  
			if(CoordTelescops[i].flag == 1){
				DrawAntenna(0,1,0);
			}else{
				DrawAntenna(1,0,0);
			}
			glColor3f(1,1,1);
			glPopMatrix();
		  }
		}

		glRotated(-EarthAngle*180/Pi,0.0,0.0,1.0);
		glRotated(-90,0.0,0.0,1.0);
		glRotated(-EarthAngleNaklona*180/Pi,0.0,1.0,0.0);
		glTranslatef(-EarthX,-EarthY,-EarthZ);

		glTranslatef(LunaX,LunaY,LunaZ);
		glRotated(LunaAngleNaklona*180/Pi,0.0,1.0,0.0);
		glRotated(LunaAngle*180/Pi,0.0,0.0,1.0);

		quadObj = gluNewQuadric();
		gluQuadricTexture(quadObj,TRUE);
		gluQuadricOrientation(quadObj,GLU_INSIDE);
        glCallList(lists[2]);
        glEnable(GL_TEXTURE_2D);
		gluSphere(quadObj,R_luna,15,15);
        glDisable(GL_TEXTURE_2D);
		gluDeleteQuadric(quadObj);
		glRotated(-LunaAngle*180/Pi,0.0,0.0,1.0);
		glRotated(-LunaAngleNaklona*180/Pi,0.0,1.0,0.0);
		glTranslatef(-LunaX,-LunaY,-LunaZ);
		glRotated(EarthAngleNaklona*180/Pi,0.0,1.0,0.0);
		glRotated(90,0.0,0.0,1.0);  //поворот так чтобы ось Х была направлена на точку весеннего равноденствия
		glTranslatef(RadioAstronX,RadioAstronY,RadioAstronZ);
		glRotated(RadioAstronAngleNaklona*180/Pi,cos(RadioAstronAngle),sin(RadioAstronAngle),0.0);

		glRotated(RadioAstronAngle*180/Pi,0.0,0.0,1.0);

//		printf("asd %i \n",Num_Squares*2+Num_Triangles);
				VECTOR Normal;
        glEnable(GL_LIGHTING);
		for(int i = 0; i<Num_Squares;i++){
			Normal = CalcNormal(Vertexs[Squares[i].v1],Vertexs[Squares[i].v2],Vertexs[Squares[i].v3]);
			glBegin(GL_POLYGON);
	        glNormal3f(Normal.x, Normal.y, Normal.z);
			glVertex3f(Vertexs[Squares[i].v1].x,Vertexs[Squares[i].v1].y,Vertexs[Squares[i].v1].z);
			glVertex3f(Vertexs[Squares[i].v2].x,Vertexs[Squares[i].v2].y,Vertexs[Squares[i].v2].z);
			glVertex3f(Vertexs[Squares[i].v3].x,Vertexs[Squares[i].v3].y,Vertexs[Squares[i].v3].z);
			glVertex3f(Vertexs[Squares[i].v4].x,Vertexs[Squares[i].v4].y,Vertexs[Squares[i].v4].z);
			glEnd();
		}
		for(int i = 0; i<Num_Triangles;i++){
			Normal = CalcNormal(Vertexs[Triangles[i].v1],Vertexs[Triangles[i].v2],Vertexs[Triangles[i].v3]);
			glBegin(GL_TRIANGLES);
	        glNormal3f(Normal.x, Normal.y, Normal.z);
			glVertex3f(Vertexs[Triangles[i].v1].x,Vertexs[Triangles[i].v1].y,Vertexs[Triangles[i].v1].z);
			glVertex3f(Vertexs[Triangles[i].v2].x,Vertexs[Triangles[i].v2].y,Vertexs[Triangles[i].v2].z);
			glVertex3f(Vertexs[Triangles[i].v3].x,Vertexs[Triangles[i].v3].y,Vertexs[Triangles[i].v3].z);
			glEnd();
		}
		for(int i = 0; i<Num_Lines-20;i++){
			glBegin(GL_LINES);
			glVertex3f(Vertexs[Lines[i].v1].x,Vertexs[Lines[i].v1].y,Vertexs[Lines[i].v1].z);
			glVertex3f(Vertexs[Lines[i].v2].x,Vertexs[Lines[i].v2].y,Vertexs[Lines[i].v2].z);
			glEnd();
		}

		glRotated(-RadioAstronAngle*180/Pi,0.0,0.0,1.0);
		glRotated(-RadioAstronAngleNaklona*180/Pi,cos(RadioAstronAngle),sin(RadioAstronAngle),0.0);
		glTranslatef(-RadioAstronX,-RadioAstronY,-RadioAstronZ);
		glRotated(-90,0.0,0.0,1.0);  
		glRotated(-EarthAngleNaklona*180/Pi,0.0,1.0,0.0);

        SwapBuffers(dc) ;
		wglMakeCurrent(0,0);
//		ReleaseDC(hWnd, dc);
		EndPaint(hWnd, &Ps);  

		return 0;
}

LRESULT CALLBACK WndProc_3d_view_earth(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam )
{
double	DeltaMouseX,DeltaMouseY;
//PAINTSTRUCT Ps;
switch (message)
	{

    case WM_MOUSEMOVE:

		DeltaMouseX = MouseX_Earth - (double)LOWORD(lParam)/(double)WndWidth(hWnd);
		DeltaMouseY = MouseY_Earth - (double)HIWORD(lParam)/(double)WndHeight(hWnd);
		DeltaMouseX = DeltaMouseX * 3.0;
		DeltaMouseY = DeltaMouseY * 3.0;
		MouseX_Earth = (double)LOWORD(lParam)/(double)WndWidth(hWnd);
		MouseY_Earth = (double)HIWORD(lParam)/(double)WndHeight(hWnd);
		double rrr,rrf,rr,rrt;
		if(MouseKey_Earth == 1){
		rrr = sqrt(CamX_Earth*CamX_Earth+CamY_Earth*CamY_Earth+CamZ_Earth*CamZ_Earth);
		rr = sqrt(CamX_Earth*CamX_Earth+CamY_Earth*CamY_Earth);
		rrf = acos(CamY_Earth/rr)*asin(CamX_Earth/rr)/abs(asin(CamX_Earth/rr));
		CamAngle_Earth =  CamAngle_Earth+DeltaMouseX;///180*3.14;
		CamAngleNaklona_Earth = CamAngleNaklona_Earth-DeltaMouseY;
		CamX_Earth = rrr*sin(CamAngle_Earth-3.14/2*3)*(cos(CamAngleNaklona_Earth));
		CamY_Earth = rrr*cos(CamAngle_Earth-3.14/2*3)*(cos(CamAngleNaklona_Earth));
		CamZ_Earth = rrr*sin(CamAngleNaklona_Earth);
		}

//		printf("*******************   %f %f ***************\n",MouseX,rrr);
  //      mousey := HIWORD(lP)/DeviceMode_real.dmPelsHeight*100;
	    break;	
    case WM_LBUTTONDOWN: 
		MouseKey_Earth = 1;
//        mousekey := wp;
  
	    break;	
	case WM_MOUSEWHEEL:
		if ( (HIWORD(wParam)>>15)  > 0 ){
			CamX_Earth = CamX_Earth -0.1*cos(CamAngle_Earth)*cos(CamAngleNaklona_Earth);
			CamY_Earth = CamY_Earth +0.1*sin(CamAngle_Earth)*cos(CamAngleNaklona_Earth);
			CamZ_Earth = CamZ_Earth -0.1*sin(CamAngleNaklona_Earth);
		}else{
			CamX_Earth = CamX_Earth +0.1*cos(CamAngle_Earth)*cos(CamAngleNaklona_Earth);
			CamY_Earth = CamY_Earth -0.1*sin(CamAngle_Earth)*cos(CamAngleNaklona_Earth);
			CamZ_Earth = CamZ_Earth +0.1*sin(CamAngleNaklona_Earth);
		}
  
	    break;	
    case WM_LBUTTONUP: 
		MouseKey_Earth = 0;
       // mousekey := 0;
       
	    break;	
	case WM_PAINT:
		wm_paint_func(hWnd, message, wParam, lParam, dc_earth, HRC_earth, CamX_Earth, CamY_Earth, CamZ_Earth ,CamAngleNaklona_Earth, CamAngle_Earth, MouseX_Earth, MouseY_Earth, MouseKey_Earth, 0);

	    break;	
	case WM_DESTROY:
	    PostQuitMessage(0);
	    break;	
	case WM_SIZE:
        wglMakeCurrent(dc_earth, HRC_earth);
		glViewport(0, 0, WndWidth(hWnd), WndHeight(hWnd));
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
		gluPerspective(90,(double)WndWidth(hWnd)/WndHeight(hWnd),0.05,25);
        glMatrixMode(GL_MODELVIEW);
		InvalidateRect(hWnd, NULL, 0); 
		break;	
	default:
		return DefWindowProc(hWnd, message, wParam, lParam);
	}
	return 0;
}

LRESULT CALLBACK WndProc_3d_view(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
//    printf("a");
double	DeltaMouseX,DeltaMouseY;
//PAINTSTRUCT Ps;
switch (message)
	{

    case WM_MOUSEMOVE:
        //int mousex = LOWORD(lParam);//DeviceMode_real.dmPelsWidth*100; WndWidth(hWnd)

		DeltaMouseX = MouseX - (double)LOWORD(lParam)/(double)WndWidth(hWnd);
		DeltaMouseY = MouseY - (double)HIWORD(lParam)/(double)WndHeight(hWnd);
		DeltaMouseX = DeltaMouseX * 3.0;
		DeltaMouseY = DeltaMouseY * 3.0;
		MouseX = (double)LOWORD(lParam)/(double)WndWidth(hWnd);
		MouseY = (double)HIWORD(lParam)/(double)WndHeight(hWnd);
		double rrr,rrf,rr,rrt;
		if(MouseKey == 1){
//			CamAngle = CamAngle+DeltaMouseX;
		rrr = sqrt(CamX*CamX+CamY*CamY+CamZ*CamZ);
		rr = sqrt(CamX*CamX+CamY*CamY);
		rrf = acos(CamY/rr)*asin(CamX/rr)/abs(asin(CamX/rr));
		//rrt = asin(CamZ/rrr);//*acos(rr/rrr)/abs(acos(rr/rrr));
		//CamAngle =  3.14/2-rrf;///180*3.14;
		CamAngle =  CamAngle+DeltaMouseX;///180*3.14;
//		CamAngle =  rrf+DeltaMouseX+3.14/2*3;///180*3.14;
//		CamAngleNaklona = rrt-DeltaMouseY;
		CamAngleNaklona = CamAngleNaklona-DeltaMouseY;
//		CamX = rrr*sin(rrf+DeltaMouseX)*cos(rrt-DeltaMouseY);
//		CamY = rrr*cos(rrf+DeltaMouseX)*cos(rrt-DeltaMouseY);
//		CamZ = rrr*sin(rrt-DeltaMouseY);
		CamX = rrr*sin(CamAngle-3.14/2*3)*(cos(CamAngleNaklona));
		CamY = rrr*cos(CamAngle-3.14/2*3)*(cos(CamAngleNaklona));
		CamZ = rrr*sin(CamAngleNaklona);
		}
		if(MouseKey == 3){
//			CamAngle = CamAngle+DeltaMouseX;
		rrr = sqrt(CamX*CamX+CamY*CamY+CamZ*CamZ);
		rr = sqrt(CamX*CamX+CamY*CamY);
		rrf = acos(CamY/rr)*asin(CamX/rr)/abs(asin(CamX/rr));
		//rrt = asin(CamZ/rrr);//*acos(rr/rrr)/abs(acos(rr/rrr));
		//CamAngle =  3.14/2-rrf;///180*3.14;
		CamAngle =  CamAngle+DeltaMouseX;///180*3.14;
//		CamAngle =  rrf+DeltaMouseX+3.14/2*3;///180*3.14;
//		CamAngleNaklona = rrt-DeltaMouseY;
		CamAngleNaklona = CamAngleNaklona-DeltaMouseY;
//		CamX = rrr*sin(rrf+DeltaMouseX)*cos(rrt-DeltaMouseY);
//		CamY = rrr*cos(rrf+DeltaMouseX)*cos(rrt-DeltaMouseY);
//		CamZ = rrr*sin(rrt-DeltaMouseY);
//		CamX = rrr*sin(CamAngle-3.14/2*3)*(cos(CamAngleNaklona));
//		CamY = rrr*cos(CamAngle-3.14/2*3)*(cos(CamAngleNaklona));
//		CamZ = rrr*sin(CamAngleNaklona);
		}

//		printf("*******************   %f %f ***************\n",MouseX,rrr);
  //      mousey := HIWORD(lP)/DeviceMode_real.dmPelsHeight*100;
	    break;	
    case WM_LBUTTONDOWN: 
		MouseKey = 1;
//        mousekey := wp;
  
	    break;	
    case WM_RBUTTONDOWN: 
		MouseKey = 3;
//        mousekey := wp;
  
	    break;	
	case WM_MOUSEWHEEL:
		if ( (HIWORD(wParam)>>15)  > 0 ){
			CamX = CamX -0.1*cos(CamAngle)*cos(CamAngleNaklona);
			CamY = CamY +0.1*sin(CamAngle)*cos(CamAngleNaklona);
			CamZ = CamZ -0.1*sin(CamAngleNaklona);
		}else{
			CamX = CamX +0.1*cos(CamAngle)*cos(CamAngleNaklona);
			CamY = CamY -0.1*sin(CamAngle)*cos(CamAngleNaklona);
			CamZ = CamZ +0.1*sin(CamAngleNaklona);
		}
	//		printf("@@@@@@@@@@@@@@@@@@@@@@%i %i %i %i@@@@@@@@@@@@@@@@\n",HIWORD(wParam)>>15,LOWORD(wParam),HIWORD(lParam),LOWORD(lParam));

//		MouseKey = 1;
//        mousekey := wp;
  
	    break;	
    case WM_LBUTTONUP: 
		MouseKey = 0;
       // mousekey := 0;
       
	    break;	
    case WM_RBUTTONUP: 
		MouseKey = 0;
       // mousekey := 0;
       
	    break;	


	case WM_TIMER:
		//printf("asd");
//		Time = Time + 0.01;
		Time = Time + dTime;
		CalculatePositions();
	   // PostQuitMessage(0);
	    break;	
	case WM_DESTROY:
		glDeleteLists(GLF_START_LIST,256);
	    PostQuitMessage(0);
	    break;	

	case WM_CREATE:
//		InitWndGraphSettings();
	    break;	
	case WM_CHAR:
		if ((char)wParam == 'x'){
			//Angle_x = Angle_x + 0.01;
			//CamAlpha = CamAlpha+0.01;
		}
		if ((char)wParam == 'y'){
			//Angle_y = Angle_y + 0.05;
			//CamTeta = CamAlpha+0.01;
		}
		if ((char)wParam == '1'){
			dTime = dTime*2;
			Time_old = dTime;
		}
		if ((char)wParam == '2'){
			dTime = dTime/2;
			Time_old = dTime;
		}
		if ((char)wParam == '3'){
			if (dTime < 0.001){
				dTime = Time_old;
			}else{
				dTime = 0;
			}
//			KillTimer(wnd,100);
		}
//		if ((char)wParam == '4'){
//			SetTimer (wnd, 100, 1, NULL);
//		}
		if ((char)wParam == '4'){
			dTime = -dTime;
		}
		if ((char)wParam == 'd'){
			CamAngle = CamAngle-0.05;
//			Angle_z = Angle_z + 0.01;
			//CamAlpha = CamAlpha+0.01;
		}
		if ((char)wParam == 'a'){
			CamAngle = CamAngle+0.05;
//			Angle_z = Angle_z - 0.01;
			//CamAlpha = CamAlpha+0.01;
		}
		if ((char)wParam == 'w'){
			CamAngleNaklona = CamAngleNaklona-0.05;
		}
		if ((char)wParam == 's'){
			CamAngleNaklona = CamAngleNaklona+0.05;
		}
		if ((char)wParam == 'q'){
			//CamX = CamX - 0.01;
			CamX = CamX -0.1*cos(CamAngle)*cos(CamAngleNaklona);
			CamY = CamY +0.1*sin(CamAngle)*cos(CamAngleNaklona);
			CamZ = CamZ -0.1*sin(CamAngleNaklona);
		}
		if ((char)wParam == 'z'){
			CamX = CamX +0.1*cos(CamAngle)*cos(CamAngleNaklona);
			CamY = CamY -0.1*sin(CamAngle)*cos(CamAngleNaklona);
			CamZ = CamZ +0.1*sin(CamAngleNaklona);
//			DRAW_EKL = 0;
		}

		if ((char)wParam == 'm'){
			if (DRAW_EKL == 1){DRAW_EKL = 0;}else{DRAW_EKL = 1;}
//			DRAW_EKL = 0;
		}
	    break;	
	case WM_SIZE:
        wglMakeCurrent(dc, HRC);
		glViewport(0, 0, WndWidth(hWnd), WndHeight(hWnd));
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
		gluPerspective(90,(double)WndWidth(hWnd)/WndHeight(hWnd),0.05,25);
        glMatrixMode(GL_MODELVIEW);
		InvalidateRect(hWnd, NULL, 0); 
		break;	
	case WM_PAINT:
		wm_paint_func(hWnd, message, wParam, lParam, dc, HRC,CamX, CamY, CamZ ,CamAngleNaklona, CamAngle, MouseX, MouseY, MouseKey,DRAW_EKL);
	    break;	
	default:
		return DefWindowProc(hWnd, message, wParam, lParam);
	}
	return 0;
}

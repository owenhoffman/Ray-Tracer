#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define HIGHT 640
#define WIDTH 640
#define numSpheres 4

struct Sphere
{
int red;
int green;
int blue;
double X;
double Y;
double Z;
double R;
double dTheta;
double dPhi;
};

struct Ray
{
double mag;
double dx;
double dy;
double dz;
};

struct Point
{
double X;
double Y;
double Z;
};

struct Image
{
int*** pixels;
int rows;
int cols;
};

void moveSphere(struct Sphere* s, double dx, double dy, double dz)
{
s->X = s->X + dx;
s->Y = s->Y + dy;
s->Z = s->Z + dz;
}

void rotateSphere(struct Sphere* s, double dT, double dP)
{
s-> dTheta += dT;
s->dPhi += dP;
if(s->dTheta < 0)
s->dTheta += 2*3.1415926535;
if(s->dPhi < 0)
s->dPhi += 2*3.1415926535;
}

struct Sphere* newSphere(int r, int g, int b, double x, double y, double z, double rad)
{
struct Sphere* s;
s->red = r;
s->green = g;
s->blue = b;
s->X = x;
s->Y = y;
s->Z = z;
s->R = rad;
s->dTheta = 0;
s->dPhi = 0;
return s;
}


//takes two points returns magnitude and unit vectors in each direction
struct Ray rayTrace(struct Point p0, struct Point p1)
{
        double rMag = sqrt(pow(p1.X-p0.X, 2) + pow(p1.Y-p0.Y, 2) + pow(p1.Z-p0.Z, 2)) ;
        double rX = (p1.X-p0.X)/rMag;
        double rY = (p1.Y-p0.Y)/rMag;
        double rZ = (p1.Z-p0.Z)/rMag;
		
        struct Ray r = {rMag, rX, rY, rZ};
		return r;
}

struct Ray sphereNormal(struct Sphere s, struct Point p)
{
	struct Point c = {s.X, s.Y, s.Z};
	struct Ray ray = rayTrace(c, p);
	return ray;
}


//struct Ray reflectedRay(

//returns the distance to the nearest sphere (t) and coordinates of intersection
double* intersect(struct Point p, struct Ray r, struct Sphere* spheres)
{
int i;
double tArr[numSpheres];

        for(i = 0; i < numSpheres; i++)
        {
        double a = 1.0;
        double b = 2*(r.dx*(p.X - spheres[i].X) + r.dy*(p.Y - spheres[i].Y) + r.dz*(p.Z-spheres[i].Z));
        double c = pow(spheres[i].X,2) + pow(spheres[i].Y,2) + pow(spheres[i].Z,2) + pow(p.X,2) + pow(p.Y,2) + pow(p.Z,2) - 2*(spheres[i].X*p.X + spheres[i].Y*p.Y + spheres[i].Z*p.Z) - pow(spheres[i].R,2);

        double t = -1;
        if(b*b - 4*a*c >= 0)
        {
                t = (-1*b - sqrt(b*b - 4*a*c))/(2*a);
        }
        tArr[i] = t;
		}
		int minIndex = 0;
        double minT = 999;
        for(i = 0; i < numSpheres; i++)
        {
                if(tArr[i] < minT && tArr[i] >= 0)
                {
                minIndex = i;
                minT = tArr[i];
                }
        }
		if(minT == 999)
		minT = -1;
		
	//double* arr = malloc(5);
	static double arr[5];
	
	double xi = p.X + tArr[minIndex] * r.dx;
	double yi = p.Y + tArr[minIndex] * r.dy;
	double zi = p.Z + tArr[minIndex] * r.dz;
	arr[0] = tArr[minIndex];
	arr[1] = minIndex;
	arr[2] = xi;
	arr[3] = yi;
	arr[4] = zi;
return arr;
}

struct Image* newImage(char* filename)
{
FILE *imgRd;
imgRd=fopen(filename,"r");
char format[10];
int imgCols, imgRows, imgMax;
fscanf(imgRd, "%s ", format);
fscanf(imgRd, "%d %d %d ", &imgCols, &imgRows, &imgMax);

struct Image* pointer = malloc(sizeof(struct Image));
int*** px = (int***)malloc(imgRows*sizeof(int**));
int x,y;
for(x = 0; x < imgRows; x++)
{
px[x] = (int**)malloc(imgCols*sizeof(int*));
for(y = 0; y < imgCols; y++)
px[x][y] = (int*)malloc(3*sizeof(int));
}
pointer->pixels = px;
pointer->rows = imgRows;
pointer->cols = imgCols;
int imgy, imgx;
for(imgy = 0; imgy < imgRows; imgy++)
{
for(imgx = 0; imgx < imgCols; imgx++)
{
int imgr, imgg, imgb;
fscanf(imgRd, "%d %d %d", &imgr, &imgg, &imgb);
pointer->pixels[imgy][imgx][0] = imgr;
pointer->pixels[imgy][imgx][1] = imgg;
pointer->pixels[imgy][imgx][2] = imgb;
}
}
return pointer;
}

double max(double a, double b)
{
if(a > b)
return a;
return b;
}

double min(double a, double b)
{
if(a < b)
return a;
return b;
}

double dot(struct Ray r0, struct Ray r1)
{
return r0.dx*r1.dx + r0.dy*r1.dy + r0.dz*r1.dz;
}

int* floorColor(struct Point p)
{
int* c = malloc(3);

c[0] = 150;
c[1] = 255;
c[2] = 50;

if((fmod(fabs(p.X), .2) > .1 && fmod(fabs(p.Z), .2) < .1) || (fmod(fabs(p.X), .2) < .1 && fmod(fabs(p.Z), .2) > .1) )
		{
		c[0] = 0; 
		}

return c;

}

void writeImage(char* filename, int*** rgb)
{
FILE* fout;
fout = fopen( filename , "w" ) ;
   //
   fprintf( fout , "P3\n" ) ;
   fprintf( fout , "%d %d\n" , WIDTH , HIGHT ) ;
   fprintf( fout , "255\n" ) ;
	int y, x;
   for( y = 0 ; y < HIGHT ; y++ )
   {
      for( x = 0 ; x < WIDTH ; x++)
      {
         fprintf( fout , "%d %d %d\n" ,
          rgb[y][x][0] , rgb[y][x][1] , rgb[y][x][2] );
      }
   }
   close( fout ) ;


}





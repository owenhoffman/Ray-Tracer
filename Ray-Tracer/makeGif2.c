#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "owengl.h"

//
#define WIDTH 640
#define HIGHT 640
#define numSpheres 4
//


int main(void)
{
struct Sphere blue = {125, 238, 205, .5, .5, 1.0/6, 1.0/6};
struct Sphere green = {82, 24, 250, 5.0/6, .5, .5, 1.0/6};
struct Sphere red = {152, 26, 227, 1.0/3, 2.0/3, 2.0/3, 1.0/3};
struct Sphere tiny = {255, 200, 255, .1, .5, .1, .1};

struct Sphere* balls = malloc(sizeof(struct Sphere) * numSpheres);
balls[0] = red;
balls[1] = green;
balls[2] = blue;
balls[3] = tiny;

double eyeX = .5;
double eyeY = .5;
double eyeZ = -1.0;

double lightX = 0.0;
double lightY = 1.0;
double lightZ = -0.5;

double screenZ = 0.0;
double floorY = 1.0/3;

double ambient = .2;

double bulb = .8;

struct Image img0 = *newImage("owen.ppm");



double angleOffset = 0;
for(angleOffset = 0; angleOffset <= 360; angleOffset+=20)
{
int*** rgb = (int***)malloc(HIGHT*sizeof(int**));
int g,h;
for(g = 0; g < HIGHT; g++)
{
rgb[g] = (int**)malloc(WIDTH*sizeof(int*));
for(h = 0; h < WIDTH; h++)
rgb[g][h] = (int*)malloc(3*sizeof(int));
}
   //
   int y , x ;
   //
   for( y = 0 ; y < HIGHT ; y++ )
   {
      for( x = 0 ; x < WIDTH ; x++)
      {
	  
        double xScale = ((double)x)/WIDTH;
        double yScale = 1.0-((double)y)/HIGHT;
		//point projected on the screen
		struct Point pixelPoint = {xScale, yScale, screenZ};
		//eye's location
		struct Point eyePoint = {eyeX, eyeY, eyeZ};
		//ray from eye to screen
        struct Ray r0 = rayTrace(eyePoint, pixelPoint);

		double* i0 = intersect(eyePoint, r0, balls);
		int minIndex = i0[1];
		double minT = i0[0];
		int inShadow = 0;
	
		
                //if intersects sphere
                if(minT >= 0)
                {
		struct Point p2 = {lightX, lightY, lightZ};
		struct Point lightPoint = {lightX, lightY, lightZ};
		struct Point p3 = {i0[2], i0[3], i0[4]};
		struct Ray r1 = rayTrace(p3, p2);
		double* i1 = intersect(p3, r1, balls);
		if(i1[0] > 0.01)
		inShadow = 1;
		struct Point surfacePoint = {eyeX + minT*r0.dx, eyeY + minT*r0.dy, eyeZ + minT*r0.dz};
		struct Ray normal = sphereNormal(balls[minIndex], surfacePoint);
          
		struct Ray surfaceToLight = rayTrace(surfacePoint, lightPoint);
		double nDotL = dot(normal, surfaceToLight);
                //double nDotL = (nX*lX + nY*lY + nZ*lZ)/(nMag*lMag);
				
				double vX = 2*nDotL*normal.dx - surfaceToLight.dx;
				double vY = 2*nDotL*normal.dy - surfaceToLight.dy;
				double vZ = 2*nDotL*normal.dz - surfaceToLight.dz;
				double vMag = sqrt(vX*vX+vY*vY+vZ*vZ);
			
			double spec = .5;
			double vDotE = -1*(r0.dx*vX/vMag + r0.dy*vY/vMag + r0.dz*vZ/vMag);
			double p = 15;
			double phong = spec*pow(max(0, vDotE),p)*bulb*255;


                rgb[y][x][0] = min(255,balls[minIndex].red*(ambient + bulb*max(0,nDotL)) + phong);
                rgb[y][x][1] = min(255,balls[minIndex].green*(ambient + bulb*max(0,nDotL)) + phong);
                rgb[y][x][2] = min(255,balls[minIndex].blue*(ambient + bulb*max(0,nDotL)) + phong);
				
				if(minIndex == 0)
				{
				struct Ray VP = {normal.mag, -1*normal.dx, -1*normal.dy, -1*normal.dz};
				
				struct Point center = {balls[0].X, balls[0].Y, balls[0].Z};
				struct Point p10 = {eyeX + minT*r0.dx, eyeY + minT*r0.dy, eyeZ + minT*r0.dz};
				struct Ray kp = rayTrace(center, p10);
				double pi = 3.14159265;	
				double phi = acos(kp.dy);
				double theta = atan2(kp.dz, kp.dx);
				double baseLongitude = 0;
				double tempX = (fmod(angleOffset*pi/360.0 + theta + pi, 2*pi))/(2*pi);
				double tempY = (fmod(angleOffset*pi/360.0 + phi, pi))/(pi);
				
				rgb[y][x][0] = min(255,img0.pixels[(int)(tempY*img0.rows-0.01)][(int)(tempX*img0.cols-.01)][0] * ((ambient + bulb*max(0, nDotL)) + phong/50));
				rgb[y][x][1] = min(255,img0.pixels[(int)(tempY*img0.rows-0.01)][(int)(tempX*img0.cols-.01)][1] * ((ambient + bulb*max(0, nDotL)) + phong/50));
				rgb[y][x][2] = min(255,img0.pixels[(int)(tempY*img0.rows-0.01)][(int)(tempX*img0.cols-.01)][2] *((ambient + bulb*max(0, nDotL)) + phong/50));

				
				}
              if(inShadow == 1)
                {
                        rgb[y][x][0] = rgb[y][x][0]/2;
                        rgb[y][x][1] = rgb[y][x][1]/2;
                        rgb[y][x][2] = rgb[y][x][2]/2;
                }

                }
                else
                {
                //hitting the floor
                if(r0.dy < 0)
                {
                //floor normal
                double fnY = 1;
                //floor to light
				double ft = (floorY - eyeY)/r0.dy;
                double flX = (eyeX + ft*r0.dx);
                double flY = (eyeY + ft*r0.dy);
                double flZ = (eyeZ + ft*r0.dz);
                double flMag = sqrt(flX*flX + flY*flY + flZ*flZ);
				
					struct Point p2 = {lightX, lightY, lightZ};
		struct Point p3 = {flX, flY, flZ};
		struct Ray r1 = rayTrace(p3, p2);
		double* i1 = intersect(p3, r1, balls);
		if(i1[0] > 0)
				inShadow = 1;
                double floorDist = sqrt((r0.dz*r0.dx + r0.dy*r0.dy + r0.dz*r0.dz)*r0.mag*r0.mag);
                double fnDotfl = flY/flMag;

		double floorRed = 150;
		double floorGreen = 255;
		double floorBlue = 50;
		if((fmod(fabs(p3.X), .2) > .1 && fmod(fabs(p3.Z), .2) < .1) || (fmod(fabs(p3.X), .2) < .1 && fmod(fabs(p3.Z), .2) > .1) )
		{
		floorRed = 0;
		}
                rgb[y][x][0] = floorRed *(ambient + bulb*(fnDotfl));
                rgb[y][x][1] = floorGreen *(ambient + bulb*(fnDotfl));
                rgb[y][x][2] = floorBlue *(ambient + bulb*(fnDotfl));

               if(inShadow == 1)
                {
                        rgb[y][x][0] = ambient*floorRed;
                        rgb[y][x][1] = ambient*floorGreen;
                        rgb[y][x][2] = ambient*floorBlue;
                }
                }
                else
                {
                rgb[y][x][0] = 0;
                rgb[y][x][1] = 0;
                rgb[y][x][2] = 0;
                }
                }
      }
   }

char filename1[3];
char filename2[3];
int f0 =  (int)fmod(angleOffset,10);
int f1 = (int)fmod(angleOffset/10, 10);
int f2 = (int)fmod(angleOffset/100, 10);
sprintf(filename1, "%d%d%d", f2, f1, f0);
sprintf(filename2, "%d%d%d.ppm", f2, f1, f0);
  writeImage(filename2, rgb);
char command[50];
sprintf(command, "convert %s.ppm %s.png", filename1, filename1);
system(command);
moveSphere(&balls[0], .03, 0, 0);
}

//system("convert *0.ppm 0.png");
system("convert -delay 10 *.png out.gif");
   return 0 ;
}
//
// end of file
//

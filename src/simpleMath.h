#ifndef SIMPLEMATH_H
#define SIMPLEMATH_H

#include <math.h>
#include <string>

//Vectors
struct Vec3
{
    double x, y, z;
    Vec3() {x=0; y=0; z=0;}
    Vec3(double a, double b, double c):x(a),y(b),z(c){}
	std::string toString() { return ("(" + std::to_string(x) + "," + std::to_string(y) + "," + std::to_string(z) + ")"); }

	Vec3 operator+(Vec3 other);
};

struct Vec4
{
    double x, y, z, w;
    Vec4() {x=0; y=0; z=0; w=0;}
    Vec4(double a, double b, double c, double d):x(a),y(b),z(c),w(d){}
};

struct Vec2
{
    double x, y;
    Vec2() {x=0; y=0;}
    Vec2(double a, double b):x(a),y(b){}
};

struct Mat3
{
    double values[9];
    Mat3()
    {
        for (unsigned i=0; i<9; i++)
            values[i]=0;
    }
};

struct Mat4
{
    double values[16];
    Mat4(double values[16])
    {
		for (unsigned i = 0; i < 16; i++)
		{
			this->values[i] = values[i];
		}
    }
};


//Vector operations
double dot2(Vec2 a, Vec2 b);

double dot3(Vec3 a, Vec3 b);

double dot4(Vec4 a, Vec4 b);

void trans(float* m, float* v, float* r);

void axisToMat(float* a, float ang, float* m);

void mult(float* m1, float* m2, float* r);
void mult(double* m1, float* m2, double* r);
void multv(float *m, float *v, float *r, float w=1.0);

void scalar(float s,float* m, float *r);

void add(float *s, float*m, float *r);

float* interpolate(float a[16], float b[16], float t);

#endif

#pragma once
#include "simpleMath.h"

class Quaternion
{
public:
	Quaternion(double x, double y, double z, double w = 0);
	Quaternion(double angle, Vec3 axis);
	~Quaternion();

	// Operator overloads
	Quaternion operator*(Quaternion other);
	Quaternion operator*(double scalar);
	Quaternion operator/(double scalar);

	// Custom operations
	static double dot(Quaternion a, Quaternion b);
	Quaternion conjugate();
	double norm();
	Quaternion normalize();

	// Conversions
	Mat4 toMatrix();
	Vec3 toEulerAngles();

private:
	double w;
	double x;
	double y;
	double z;
};


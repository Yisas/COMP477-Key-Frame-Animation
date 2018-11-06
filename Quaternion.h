#pragma once
#include "simpleMath.h"

class Quaternion
{
public:
	Quaternion(double x, double y, double z, double w = 0);
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

private:
	double w;
	double x;
	double y;
	double z;
};


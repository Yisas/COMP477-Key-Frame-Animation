#pragma once
#include "simpleMath.h"
#include <iostream>
#include <string>

class Quaternion
{
public:
	Quaternion();
	Quaternion(double x, double y, double z, double w);
	Quaternion(double angle, Vec3 axis);
	Quaternion(Vec3 point);
	Quaternion(float eulerX, float eulerY, float eulerZ);
	Quaternion(float rotationMatrix[16]);
	Quaternion(std::string quaternion);
	~Quaternion();
	std::string toString() { return ("(" + std::to_string(w) + "," + std::to_string(x) + "i," + std::to_string(y) + "j," + std::to_string(z) + "k)"); };

	// Attribute accesors
	double getW() { return w; }
	double getX() { return x; }
	double getY() { return y; }
	double getZ() { return z; }

	// Operator overloads
	Quaternion operator*(Quaternion other);
	Quaternion operator+(Quaternion other);
	Quaternion operator-(Quaternion other);
	Quaternion operator*(float scalar);
	Quaternion operator/(float scalar);

	// Custom operations
	static double dot(Quaternion a, Quaternion b);
	Quaternion conjugate();
	double norm();
	Quaternion normalize();
	Vec3 rotatePoint(Vec3 point);

	// Conversions
	Mat4 toMatrix();
	float* toFloatMatrix();
	Vec3 toEulerAngles();
	static Quaternion rotationMatrixToQuaternion(float* rotationMatrix);
	static Quaternion interpolateLineraly(Quaternion a, Quaternion b, float t);
	static Quaternion SLERP(Quaternion a, Quaternion b, float t);

	//Misc
	bool isUnit() { return (norm() == 1); };

private:
	double w;
	double x;
	double y;
	double z;
};


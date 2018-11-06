#include "Quaternion.h"
#define _USE_MATH_DEFINES
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

Quaternion::Quaternion(double x, double y, double z, double w)
{
	this->w = w;
	this->x = x;
	this->y = y;
	this->z = z;
}

Quaternion::Quaternion(double angle, Vec3 axis)
{
	this->x = axis.x * sin(angle / 2);
	this->y = axis.y * sin(angle / 2);
	this->z = axis.z * sin(angle / 2);
	this->w = cos(angle / 2);
}

// Will produce a point form quaternion with w = 0
Quaternion::Quaternion(Vec3 point)
{
	this->w = 0;
	this->x = point.x;
	this->y = point.y;
	this->z = point.z;
}

Quaternion::~Quaternion()
{
}

Quaternion Quaternion::operator*(Quaternion other)
{
	return Quaternion(
		w*other.x + x * other.w + y * other.z - z * other.y,	//x
		w*other.y + y * other.w + z * other.x - x * other.z,	//y
		w*other.z + z * other.w + x * other.y - y * other.x,	//z
		w*other.w - x * other.x - y * other.y - z * other.z		//w
	);
}

Quaternion Quaternion::operator*(double scalar)
{
	return Quaternion(x * scalar, y * scalar, z * scalar, w * scalar);
}

Quaternion Quaternion::operator/(double scalar)
{
	return Quaternion(x / scalar, y / scalar, z / scalar, w / scalar);
}

double Quaternion::dot(Quaternion a, Quaternion b)
{
	return (
		a.x * b.x +
		a.y * b.y +
		a.z * b.z +
		a.w * b.w
		);
}

Quaternion Quaternion::conjugate()
{
	return Quaternion(-x, -y, -z, w);
}

double Quaternion::norm()
{
	return sqrt(w*w + x * x + y * y + z * z);
}

Quaternion Quaternion::normalize()
{
	return Quaternion(*this / norm());
}

Vec3 Quaternion::rotatePoint(Vec3 point)
{
	if (!isUnit())
		std::cout << "Warning! Rotating point " + point.toString() + " with a non-unit quaternion " + this->toString() + "!" << std::endl;

	// Convert point to a quaternion with w = 0
	Quaternion pointInQuatForm = Quaternion(point);
	Quaternion qv = *this * pointInQuatForm;				// q x v
	Quaternion rotatedPoint = qv * this->conjugate();		// (q x v) x q^-1

	// TODO: check that w was indeed = 0 after operations?

	return Vec3(rotatedPoint.x, rotatedPoint.y, rotatedPoint.z);
}

Mat4 Quaternion::toMatrix()
{
	double values[16] = {
		1 - (2 * y*y) - (2 * z*z), 2 * x*y - 2 * z*w, 2 * x*z + 2 * y*w, 0,
		2 * x*y + 2 * z*w, 1 - 2 * x*x - 2 * z*z, 2 * y*z - 2 * x*w, 0,
		2 * x*z - 2 * y*w, 2 * y*z + 2 * x*w, 1 - 2 * x*x - 2 * y*y, 0,
		0, 0, 0, 1
	};

	return Mat4(values);
}

Vec3 Quaternion::toEulerAngles()
{
	// Atan2 implementation is necessary because normal atan only returns  -pi/2<value<pi/2
	Vec3 result;

	// roll (x-axis rotation)
	result.x = atan2(2 * (w * x + y * z), 1 - 2 * (x * x + y * y));

	// pitch (y-axis rotation)
	double sinp = 2 * (w * y - z * x);
	if (fabs(sinp) >= 1)
		result.y = copysign(M_PI / 2, sinp); // use 90 degrees if out of range
	else
		result.y = asin(sinp);

	// yaw (z-axis rotation)
	result.z = atan2(2 * (w * z + x * y), 1 - 2 * (y * y + z * z));

	return result;
}
#include "Quaternion.h"
#include <cmath>

Quaternion::Quaternion(double x, double y, double z, double w = 0)
{
	this->w = w;
	this->x = x;
	this->y = y;
	this->z = z;
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
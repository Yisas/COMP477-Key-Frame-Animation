#include "Quaternion.h"
#define _USE_MATH_DEFINES
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

Quaternion::Quaternion()
{
	this->w = 1.0;
	this->x = 0.0;
	this->y = 0.0;
	this->z = 0.0;
}

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

/**
Will produce a point form quaternion with w = 0
**/
Quaternion::Quaternion(Vec3 point)
{
	this->w = 0;
	this->x = point.x;
	this->y = point.y;
	this->z = point.z;
}

Quaternion::Quaternion(float rotationMatrix[16])
{
	this->rotationMatrixToQuaternion(rotationMatrix);
}

/**
String format should coincide with toString method
**/
Quaternion::Quaternion(std::string quaternion)
{
	std::size_t current, previous = 0;
	current = quaternion.find(",");
	w = stod(quaternion.substr(1, current - 1));					// Exclude "(" at beginning 
	previous = current + 1;
	current = quaternion.find(",", previous);
	x = stod(quaternion.substr(previous, current - previous - 1));	// Exclude "(" at beginning and symbol at end
	previous = current + 1;
	current = quaternion.find(",", previous);
	y = stod(quaternion.substr(previous, current - previous - 1));	// Exclude "(" at beginning and symbol at end
	previous = current + 1;
	current = quaternion.find(",", previous);
	z = stod(quaternion.substr(previous, current - previous - 1));	// Exclude "(" at beginning and symbol at end
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

float * Quaternion::toFloatMatrix()
{
	float values[16] = {
		1 - (2 * y*y) - (2 * z*z), 2 * x*y - 2 * z*w, 2 * x*z + 2 * y*w, 0,
		2 * x*y + 2 * z*w, 1 - 2 * x*x - 2 * z*z, 2 * y*z - 2 * x*w, 0,
		2 * x*z - 2 * y*w, 2 * y*z + 2 * x*w, 1 - 2 * x*x - 2 * y*y, 0,
		0, 0, 0, 1
	};

	return values;
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

Quaternion Quaternion::rotationMatrixToQuaternion(float* rotationMatrix)
{
	// Calculate the trace
	float t = rotationMatrix[0] + rotationMatrix[5] + rotationMatrix[10];

	if (t > 0) {
		float S = sqrt(t + 1.0) * 2;						// S=4*qw 
		return Quaternion(
			(rotationMatrix[9] - rotationMatrix[6]) / S,	// X = (m21 - m12) / S
			(rotationMatrix[2] - rotationMatrix[8]) / S,	// Y = (m02 - m20) / S
			(rotationMatrix[4] - rotationMatrix[1]) / S,	// Z = (m10 - m01) / S
			0.25 / S
		);
	}
	else if ((rotationMatrix[0] > rotationMatrix[5]) && (rotationMatrix[0] > rotationMatrix[10])) {	// (m00 > m11)&(m00 > m22)
		float S = sqrt(1.0 + rotationMatrix[0] - rotationMatrix[5] - rotationMatrix[10]) * 2;		// S=4*qx 
		return Quaternion(
			0.25 * S,
			(rotationMatrix[1] + rotationMatrix[5]) / S,	// Y = (m01 + m10) / S
			(rotationMatrix[2] + rotationMatrix[8]) / S,	// Z = (m02 + m20) / S
			(rotationMatrix[9] - rotationMatrix[6]) / S		// w = (m21 - m12) / S
		);
	}
	else if (rotationMatrix[5] > rotationMatrix[10]) {		// m11 > m22
		float S = sqrt(1.0 + rotationMatrix[5] - rotationMatrix[0] - rotationMatrix[10]) * 2;	// S=4*qy
		return Quaternion(
			(rotationMatrix[1] + rotationMatrix[4]) / S,	// X = (m01 + m10) / S
			0.25 * S,										// Y = 0.25 * S;
			(rotationMatrix[6] + rotationMatrix[9]) / S,	// Z = (m12 + m21) / S;
			(rotationMatrix[2] - rotationMatrix[8]) / S		// W = (m02 - m20) / S
		);
	}
	else {
		float S = sqrt(1.0 + rotationMatrix[10] - rotationMatrix[0] - rotationMatrix[5]) * 2; // S=4*qz
		return Quaternion(
			(rotationMatrix[2] + rotationMatrix[8]) / S,	// X = (m02 + m20) / S
			(rotationMatrix[6] + rotationMatrix[9]) / S,	// Y = (m12 + m21) / S
			0.25 * S,										// Z = 0.25 * S
			(rotationMatrix[4] - rotationMatrix[1]) / S		// W = (m10 - m01) / S;
		);
	}
}

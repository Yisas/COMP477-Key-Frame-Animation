#pragma once
class Quaternion
{
public:
	Quaternion(double x, double y, double z, double w = 0);
	~Quaternion();

	Quaternion operator*(Quaternion other);
	Quaternion operator*(double scalar);
	Quaternion operator/(double scalar);

	static double dot(Quaternion a, Quaternion b);
	Quaternion conjugate();
	double norm();
	Quaternion normalize();

private:
	double w;
	double x;
	double y;
	double z;
};


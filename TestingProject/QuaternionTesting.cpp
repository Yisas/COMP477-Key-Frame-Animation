#include "stdafx.h"
#include "CppUnitTest.h"
#include "../src/Quaternion.h"
#include "../src/Quaternion.cpp"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace TestingProject
{
	TEST_CLASS(TestingQuaternion)
	{
	public:

		TEST_METHOD(TestIsUnitQuaternion)
		{
			Quaternion q1 = Quaternion(0.0, 0.0, 0.0, 1.0);
			Quaternion q2 = Quaternion(0.0, 0.0, 1.0, 0.0);
			Quaternion q3 = Quaternion(0.0, 1.0, 0.0, 0.0);
			Quaternion q4 = Quaternion(1.0, 0.0, 0.0, 0.0);
			Quaternion q5 = Quaternion(0.0, 0.0, 0.0, -1.0);
			Quaternion q6 = Quaternion(sqrt(2.0) / 3.0, sqrt(3.0) / 3.0, (1.0 / 3.0), sqrt(3.0) / 3.0);
			Quaternion q7 = Quaternion(sqrt(3.0) / 3.0, sqrt(3.0) / 3.0, 0, sqrt(3.0) / 3.0);
			Quaternion q8 = Quaternion(sqrt(2.0) / 3.0, sqrt(3.0) / 3.0, 0, sqrt(3.0) / 3.0);

			Assert::IsTrue(q1.isUnit());
			Assert::IsTrue(q2.isUnit());
			Assert::IsTrue(q3.isUnit());
			Assert::IsTrue(q4.isUnit());
			Assert::IsTrue(q5.isUnit());
			Assert::IsTrue(q6.isUnit());
			Assert::IsTrue(q7.isUnit());
			Assert::IsFalse(q8.isUnit());
		}

		TEST_METHOD(TestRotationMatrixToQuaternion) {
			float testMatrix[16] = { 0.0760397, 0.906544, -0.415206, 0, -0.994975, 0.0417863,-0.0909827,0,-0.0651299,0.420038,0.905166,0,0,0,0,1 };

			Quaternion quat = Quaternion(testMatrix);
			Assert::AreEqual(quat.getW(), 0.71115964651107788);
			Assert::AreEqual(quat.getX(), 0.17964346706867218);
			Assert::AreEqual(quat.getY(), -0.12306522578001022);
			Assert::AreEqual(quat.getZ(), -0.66845715045928955);
		}
	};
}
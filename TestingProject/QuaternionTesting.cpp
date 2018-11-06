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
			Quaternion q6 = Quaternion(sqrt(2.0)/3.0, sqrt(3.0)/3.0, (1.0/3.0), sqrt(3.0)/3.0);
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

	};
}
#include "G_Z.h"
#include <iostream>
#include <cmath>
#include <vector>
#define N 3 // ����������� ������� ���������
using namespace std;
namespace G_Z
{
	vector<double> grad(vector<double>(3, 0));
	vector<vector<double>> points(3, vector<double>(3, 0));
	vector<double> e(vector<double>(3, 0));
	vector<double> point(vector<double>(3, 0));//�������� �����������
	vector<double> point_i(vector<double>(3, 0));//�������� �����������
	int k;//����� �������� ������ �����
	int j;//����� ����� ����������
	int n;//������������ ����� ��������
	const double ac = 0.00001;
	double d12;
	double d13;
	double d32;
	void step1();
	void step2();
	void step3();
	void step4();
	void step5();
	void find_grad_in_point();
	double funkx();//��� �� �
	double funky();//��� �� y
	double funkz();//��� �� z
	void ras();
	
	void All(double t1, double t2, double t3, double x1, double x2, double x3, double y1, double y2, double y3, double z1, double z2, double z3) {

		points[0][0] = x1; points[0][1] = y1; points[0][2] = z1;
		points[1][0] = x2; points[1][1] = y2; points[1][2] = z2;
		points[2][0] = x3; points[2][1] = y3; points[2][2] = z3;
		d12 = t1;
		d13 = t2;
		d32 = t3;

		n = 1000;//������������ ����� ��������
		j = 0;
		//��� 1
		step1();		//��� 2
				//��� 3
		//��� 4
		//��� 5
		//��� 6
		//��� 7
	}


	void step1()
	{
		k = 0;//����� �������� � ������;
		step2();
	}
	
	void step2()
	{
		if (k < n)
		{
			step3();
		}
		else
		{
			j++;
			step1();
		}
	}
	
	void step3()
	{
		find_grad_in_point();
		step4();
	}

	void step4()
	{
		if (sqrt((pow(grad[0], 2) + pow(grad[1], 2) + pow(grad[2], 2))) < ac)
		{
			//�����??
		}
		else
		{
			step5();
		}
	}
	void step5()
	{
		ras();
		k++;
		step2();
	}

	void find_grad_in_point()
	{
		
			grad[0] = funkx();
			grad[1] = funky();
			grad[2] = funkz();
		////��� 2 ���������
		//grad[1][0] = funkx(points[1]);
		//grad[1][1] = funky(points[1]);
		//grad[1][2] = funkz(points[1]);
		////��� 3 ���������
		//grad[2][0] = funkx(points[2]);
		//grad[2][1] = funky(points[2]);
		//grad[2][2] = funkz(points[2]);

	}
	double funkx()//��� �� �
	{
		return ((2 * (point[0] - points[0][0]) / (sqrt(pow(point[0] - points[0][0], 2) + pow(point[1] - points[0][1], 2) + pow(point[2] - points[0][2], 2)))
			- 2 * (point[0] - points[1][0]) / (sqrt(pow(point[0] - points[1][0], 2) + pow(point[1] - points[1][1], 2) + pow(point[2] - points[1][2], 2))))
			* (sqrt(pow(point[0] - points[0][0], 2) + pow(point[1] - points[0][1], 2) + pow(point[2] - points[0][2], 2))
				- sqrt(pow(point[0] - points[1][0], 2) + pow(point[1] - points[1][1], 2) + pow(point[2] - points[1][2], 2))
				- d12)
			+
			(2 * (point[0] - points[0][0]) / (sqrt(pow(point[0] - points[0][0], 2) + pow(point[1] - points[0][1], 2) + pow(point[2] - points[0][2], 2)))
				- 2 * (point[0] - points[2][0]) / (sqrt(pow(point[0] - points[2][0], 2) + pow(point[1] - points[2][1], 2) + pow(point[2] - points[2][2], 2))))
			* (sqrt(pow(point[0] - points[0][0], 2) + pow(point[1] - points[0][1], 2) + pow(point[2] - points[0][2], 2))
				- sqrt(pow(point[0] - points[2][0], 2) + pow(point[1] - points[2][1], 2) + pow(point[2] - points[2][2], 2))
				- d13) 
			+
			(2 * (point[0] - points[2][0]) / (sqrt(pow(point[0] - points[2][0], 2) + pow(point[1] - points[2][1], 2) + pow(point[2] - points[2][2], 2)))
				- 2 * (point[0] - points[1][0]) / (sqrt(pow(point[0] - points[1][0], 2) + pow(point[1] - points[1][1], 2) + pow(point[2] - points[1][2], 2))))
			* (sqrt(pow(point[0] - points[2][0], 2) + pow(point[1] - points[2][1], 2) + pow(point[2] - points[2][2], 2))
				- sqrt(pow(point[0] - points[1][0], 2) + pow(point[1] - points[1][1], 2) + pow(point[2] - points[1][2], 2))
				- d32));
	}

	double funky()//��� �� �
	{
		return ((2 * (point[1] - points[0][1]) / (sqrt(pow(point[0] - points[0][0], 2) + pow(point[1] - points[0][1], 2) + pow(point[2] - points[0][2], 2)))
			- 2 * (point[1] - points[1][1]) / (sqrt(pow(point[0] - points[1][0], 2) + pow(point[1] - points[1][1], 2) + pow(point[2] - points[1][2], 2))))
			* (sqrt(pow(point[0] - points[0][0], 2) + pow(point[1] - points[0][1], 2) + pow(point[2] - points[0][2], 2))
				- sqrt(pow(point[0] - points[1][0], 2) + pow(point[1] - points[1][1], 2) + pow(point[2] - points[1][2], 2))
				- d12)
			+
			(2 * (point[1] - points[0][1]) / (sqrt(pow(point[0] - points[0][0], 2) + pow(point[1] - points[0][1], 2) + pow(point[2] - points[0][2], 2)))
				- 2 * (point[1] - points[2][1]) / (sqrt(pow(point[0] - points[2][0], 2) + pow(point[1] - points[2][1], 2) + pow(point[2] - points[2][2], 2))))
			* (sqrt(pow(point[0] - points[0][0], 2) + pow(point[1] - points[0][1], 2) + pow(point[2] - points[0][2], 2))
				- sqrt(pow(point[0] - points[2][0], 2) + pow(point[1] - points[2][1], 2) + pow(point[2] - points[2][2], 2))
				- d13)
			+
			(2 * (point[1] - points[2][1]) / (sqrt(pow(point[0] - points[2][0], 2) + pow(point[1] - points[2][1], 2) + pow(point[2] - points[2][2], 2)))
				- 2 * (point[1] - points[1][1]) / (sqrt(pow(point[0] - points[1][0], 2) + pow(point[1] - points[1][1], 2) + pow(point[2] - points[1][2], 2))))
			* (sqrt(pow(point[0] - points[2][0], 2) + pow(point[1] - points[2][1], 2) + pow(point[2] - points[2][2], 2))
				- sqrt(pow(point[0] - points[1][0], 2) + pow(point[1] - points[1][1], 2) + pow(point[2] - points[1][2], 2))
				- d32));
	}

	double funkz()//��� �� �
	{
		return ((2 * (point[2] - points[0][2]) / (sqrt(pow(point[0] - points[0][0], 2) + pow(point[1] - points[0][1], 2) + pow(point[2] - points[0][2], 2)))
			- 2 * (point[2] - points[1][2]) / (sqrt(pow(point[0] - points[1][0], 2) + pow(point[1] - points[1][1], 2) + pow(point[2] - points[1][2], 2))))
			* (sqrt(pow(point[0] - points[0][0], 2) + pow(point[1] - points[0][1], 2) + pow(point[2] - points[0][2], 2))
				- sqrt(pow(point[0] - points[1][0], 2) + pow(point[1] - points[1][1], 2) + pow(point[2] - points[1][2], 2))
				- d12)
			+
			(2 * (point[2] - points[0][2]) / (sqrt(pow(point[0] - points[0][0], 2) + pow(point[1] - points[0][1], 2) + pow(point[2] - points[0][2], 2)))
				- 2 * (point[2] - points[2][2]) / (sqrt(pow(point[0] - points[2][0], 2) + pow(point[1] - points[2][1], 2) + pow(point[2] - points[2][2], 2))))
			* (sqrt(pow(point[0] - points[0][0], 2) + pow(point[1] - points[0][1], 2) + pow(point[2] - points[0][2], 2))
				- sqrt(pow(point[0] - points[2][0], 2) + pow(point[1] - points[2][1], 2) + pow(point[2] - points[2][2], 2))
				- d13)
			+
			(2 * (point[2] - points[2][2]) / (sqrt(pow(point[0] - points[2][0], 2) + pow(point[1] - points[2][1], 2) + pow(point[2] - points[2][2], 2)))
				- 2 * (point[2] - points[1][2]) / (sqrt(pow(point[0] - points[1][0], 2) + pow(point[1] - points[1][1], 2) + pow(point[2] - points[1][2], 2))))
			* (sqrt(pow(point[0] - points[2][0], 2) + pow(point[1] - points[2][1], 2) + pow(point[2] - points[2][2], 2))
				- sqrt(pow(point[0] - points[1][0], 2) + pow(point[1] - points[1][1], 2) + pow(point[2] - points[1][2], 2))
				- d32));
	}
	
	void set_e()
	{		
			for (int i = 0; i < point.size(); ++i)
			{
				if (i == j) e[i] = 1;
				else e[i] = 0;
			}
	}
	double getH(vector<double> p)
	{
		return (sqrt(pow(p[0] - points[0][0], 2) + pow(p[1] - points[0][1], 2) + pow(p[2] - points[0][2], 2)) +
			sqrt(pow(p[0] - points[1][0], 2) + pow(p[1] - points[1][1], 2) + pow(p[2] - points[1][2], 2)) +
			sqrt(pow(p[0] - points[2][0], 2) + pow(p[1] - points[2][1], 2) + pow(p[2] - points[2][2], 2)) - d12 - d32 - d13);
	}
	double LA = 10;
	void useLA()
	{
		point_i = point;
		point_i[j] = point[k] - LA * grad[j];
	}
	void ras()
	{
		useLA();
		if (getH(point) > getH(point_i))
		{
			point = point_i;
		}
		else
		{
			LA = LA / 1.05;
			useLA();
			if (getH(point) > getH(point_i))
			{
				point = point_i;
			}
			else
			{
				LA = LA /- 1.05;
				useLA();
				if (getH(point) > getH(point_i))
				{
					point = point_i;
				}
			}
		}
	}




//{
//    // ������� ��� ���������� �������� ���������� xi
//    double calculateXi(double A[N][N], double x[N], double b[N], int i) {
//        double sum = 0.0;
//        for (int j = 0; j < N; ++j) {
//            if (j != i) {
//                sum += A[i][j] * x[j];
//            }
//        }
//        return (b[i] - sum) / A[i][i];
//    }
//
//    // ����� ��������������� ������ ������ � �������
//    void gaussSeidel(double A[N][N], double x[N], double b[N], int maxIterations, double tolerance) {
//        double xNew[N]; // ������ ��� �������� ����� �������� ����������
//        int iteration = 0;
//        double error = tolerance + 1;
//
//        while (iteration < maxIterations && error > tolerance) {
//            error = 0.0;
//            for (int i = 0; i < N; ++i) {
//                double xiNew = calculateXi(A, xNew, b, i);
//                error = std::max(error, std::abs(xiNew - x[i]));
//                x[i] = xiNew;
//            }
//            ++iteration;
//        }
//
//        if (error <= tolerance) {
//            std::cout << "����� ������� � �������:" << std::endl;
//            for (int i = 0; i < N; ++i) {
//                std::cout << "x[" << i << "] = " << x[i] << std::endl;
//            }
//        }
//        else {
//            std::cout << "����� �� ������� � ������� ����� " << maxIterations << " ��������." << std::endl;
//        }
//    }
//
//    void All(double t1, double t2, double t3, double x1, double x2, double x3, double y1, double y2, double y3, double z1, double z2, double z3)
//    {
//        double A[N][N] = { {x1, y1, z1}, {x2, y2, z2}, {x3, y3, z3} }; // ������� �������������
//        double b[N] = { t1, t2, t3 }; // ������ ������ �����
//        double x[N] = { 0 }; // ��������� ����������� �������
//        int maxIterations = 1000; // ������������ ���������� ��������
//        double tolerance = 1e-4; // �����������
//
//        gaussSeidel(A, x, b, maxIterations, tolerance);
//
//    }
};

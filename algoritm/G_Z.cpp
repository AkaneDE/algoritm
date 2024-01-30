#include "G_Z.h"
#include <iostream>
#include <cmath>
#include <vector>
#define N 3 // Размерность системы уравнений
using namespace std;
namespace G_Z
{
	vector<double> grad(vector<double>(3, 0));
	vector<vector<double>> points(3, vector<double>(3, 0));
	vector<double> point;
	int k;//номер итерации внутри цикла
	int j;//номер цикла вычислений
	int n;//максимальное число итераций
	const double ac = 0.00001;
	void All(double t1, double t2, double t3, double x1, double x2, double x3, double y1, double y2, double y3, double z1, double z2, double z3) {

		points[0][0] = x1; points[0][1] = y1; points[0][2] = z1;
		points[1][0] = x2; points[1][1] = y2; points[1][2] = z2;
		points[2][0] = x3; points[2][1] = y3; points[2][2] = z3;
		n = 1000;//максимальное число итераций
		j = 0;
		//шаг 1
		k = 0;//число итераций в начале 1;
		//шаг 2
				//шаг 3
		//шаг 4
		//шаг 5
		//шаг 6
		//шаг 7
	}


	void step1()
	{
		k = 0;//число итераций в начале;
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
	}

	void step4()
	{
		if (sqrt((pow(grad[0], 2) + pow(grad[1], 2) + pow(grad[2], 2))) < ac)
		{
			//конец??
		}
		else
		{
			step5();
		}
	}

	void step5()
	{

	}

	void find_grad_in_point()
	{
		
		if (j < 2)
		{
			grad[0] = funkx(points[k],points[k+1]);
			grad[1] = funky(points[k],points[k+1]);
			grad[2] = funkz(points[k],points[k+1]);
		}
		else
		{
			grad[0] = funkx(points[k], points[0]);
			grad[1] = funky(points[k], points[0]);
			grad[2] = funkz(points[k], points[0]);
		}
		////для 2 микрофона
		//grad[1][0] = funkx(points[1]);
		//grad[1][1] = funky(points[1]);
		//grad[1][2] = funkz(points[1]);
		////для 3 микрофона
		//grad[2][0] = funkx(points[2]);
		//grad[2][1] = funky(points[2]);
		//grad[2][2] = funkz(points[2]);

	}
	double funkx(vector<double> c1, vector<double> c2)//диф по х
	{
		return (point[0] - c1[0]) / (sqrt(pow(point[0] - c1[0], 2) + pow(point[1] - c1[1], 2) + pow(point[2] - c1[2], 2)))
			 + (point[0] - c2[0]) / (sqrt(pow(point[0] - c2[0], 2) + pow(point[1] - c2[1], 2) + pow(point[2] - c2[2], 2)));
	}
	double funky(vector<double> c1, vector<double> c2)//диф по у
	{
		return (point[1] - c1[1]) / (sqrt(pow(point[0] - c1[0], 2) + pow(point[1] - c1[1], 2) + pow(point[2] - c1[2], 2)))
			 + (point[1] - c2[1]) / (sqrt(pow(point[0] - c2[0], 2) + pow(point[1] - c2[1], 2) + pow(point[2] - c2[2], 2)));
	}
	double funkz(vector<double> c1, vector<double> c2)//диф по z
	{
		return (point[2] - c1[2]) / (sqrt(pow(point[0] - c1[0], 2) + pow(point[1] - c1[1], 2) + pow(point[2] - c1[2], 2)))
			+ (point[2] - c2[2]) / (sqrt(pow(point[0] - c2[0], 2) + pow(point[1] - c2[1], 2) + pow(point[2] - c2[2], 2)));
	}





//{
//    // Функция для нахождения значения переменной xi
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
//    // Метод покоординатного спуска Гаусса — Зейделя
//    void gaussSeidel(double A[N][N], double x[N], double b[N], int maxIterations, double tolerance) {
//        double xNew[N]; // Массив для хранения новых значений переменных
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
//            std::cout << "Метод сошелся к решению:" << std::endl;
//            for (int i = 0; i < N; ++i) {
//                std::cout << "x[" << i << "] = " << x[i] << std::endl;
//            }
//        }
//        else {
//            std::cout << "Метод не сошелся к решению после " << maxIterations << " итераций." << std::endl;
//        }
//    }
//
//    void All(double t1, double t2, double t3, double x1, double x2, double x3, double y1, double y2, double y3, double z1, double z2, double z3)
//    {
//        double A[N][N] = { {x1, y1, z1}, {x2, y2, z2}, {x3, y3, z3} }; // Матрица коэффициентов
//        double b[N] = { t1, t2, t3 }; // Вектор правой части
//        double x[N] = { 0 }; // Начальное приближение решения
//        int maxIterations = 1000; // Максимальное количество итераций
//        double tolerance = 1e-4; // Погрешность
//
//        gaussSeidel(A, x, b, maxIterations, tolerance);
//
//    }
};

#include "G_Z.h"
#include <iostream>
#include <cmath>
#include <vector>
#define N 3 // Ðàçìåðíîñòü ñèñòåìû óðàâíåíèé
using namespace std;
namespace G_Z
{
	struct Point
	{
	public:
		Point(double x1, double y1, double z1)  // конструктор объекта
		{
			x = x1; y = y1; z = z1;
		}
		double x;
		double y;
		double z;
	};
	double sko, sko_last;
	double SKO();
	vector<double> grad(vector<double>(3, 0));
	vector<Point> points;
	Point point(0, 0, 0);//à÷àëüíîå ïðèáëåæåíèå
	Point point_i(0, 0, 0);//à÷àëüíîå ïðèáëåæåíèå
	int k;//íîìåð èòåðàöèè âíóòðè öèêëà
	int j;//íîìåð öèêëà âû÷èñëåíèé
	int n;//ìàêñèìàëüíîå ÷èñëî èòåðàöèé
	const double ac = 0.000000001;
	double d12;
	double d13;
	double d32;
	double getH(Point p);
	/*void step1();
	void step2();
	void step3();
	void step4();
	void step5();
	void find_grad_in_point();
	double funkx();
	double funky();
	double funkz();*/
	void ras();
	void WhyNot();


	void All(double t1, double t2, double t3, double x1, double x2, double x3, double y1, double y2, double y3, double z1, double z2, double z3) {

		points.push_back(Point(x1, y1, z1));
		points.push_back(Point(x2, y2, z2));
		points.push_back(Point(x3, y3, z3));

		d12 = t1;
		d13 = t2;
		d32 = t3;
		point.x = 1000;
		point.y = 1000;
		point.z = 1000;
		point_i = point;
		sko = SKO();
		sko_last = sko;
		n = 1;
		j = 0;
		k = 0;
		WhyNot();
		//step1();		
	}
	double LA = 4;

	void WhyNot()
	{
		do
		{
			
			
			if (j < 3 - 1)
			{
				j++;
			}
			else j = 0;
			

			ras();//основной расчёт
			//вывод промежуточных значений
			cout << point.x << " " << point.y << " " << point.z << endl;
			cout << LA << endl;
			cout << getH(point) << endl;
		} while (getH(point) > ac);//Пока мы не найдём нужную точку
		cout << point.x << " " << point.y << " " << point.z << endl;
		cout << getH(point) << endl;
		cout << "Ну типа";

	}

	double getH(Point p)
	{
		return abs(sqrt(pow(p.x - points[0].x, 2) + pow(p.y - points[0].y, 2) + pow(p.z - points[0].z, 2)) - sqrt(pow(p.x - points[1].x, 2) + pow(p.y - points[1].y, 2) + pow(p.z - points[1].z, 2)) - d12
			+ sqrt(pow(p.x - points[0].x, 2) + pow(p.y - points[0].y, 2) + pow(p.z - points[0].z, 2)) - sqrt(pow(p.x - points[2].x, 2) + pow(p.y - points[2].y, 2) + pow(p.z - points[2].z, 2)) - d13
			+ sqrt(pow(p.x - points[2].x, 2) + pow(p.y - points[2].y, 2) + pow(p.z - points[2].z, 2)) - sqrt(pow(p.x - points[1].x, 2) + pow(p.y - points[1].y, 2) + pow(p.z - points[1].z, 2)) - d32);
	}
	
	double Fd12()
	{
		d12 - pow(sqrt(pow(point.x - points[0].x, 2) + pow(point.y - points[0].y, 2) + pow(point.z - points[0].z, 2)) - sqrt(pow(point.x - points[1].x, 2) + pow(point.y - points[1].y, 2) + pow(point.z - points[1].z, 2)), 2);
	}
	double Fd13()
	{
		d13 - pow(sqrt(pow(point.x - points[0].x, 2) + pow(point.y - points[0].y, 2) + pow(point.z - points[0].z, 2)) - sqrt(pow(point.x - points[2].x, 2) + pow(point.y - points[2].y, 2) + pow(point.z - points[2].z, 2)), 2);
	}
	double Fd32()
	{
		d32 - pow(sqrt(pow(point.x - points[2].x, 2) + pow(point.y - points[2].y, 2) + pow(point.z - points[2].z,2)) - sqrt(pow(point.x - points[1].x, 2) + pow(point.y - points[1].y, 2) + pow(point.z - points[1].z, 2)), 2);
	}
	double SKO()
	{
		return Fd12() + Fd13() + Fd32();
	}
	void useLA()
	{
		point_i = point;
		if (j == 0)
			point_i.x = point.x - LA;
		else if (j == 1)
			point_i.y = point.y - LA;
		else
			point_i.z = point.z - LA;


	}
	void ras()
	{
		useLA();//Расчёт координаты с заданным коэф шага
		sko = SKO();
		if (getH(point) > getH(point_i))//Перебор коэффициента и поск поддходящего
		{
			point = point_i;
		}
		else
		{
			LA = LA * -1;
			useLA();
			if (getH(point) > getH(point_i))
			{
				point = point_i;

			}
			else {
				LA = LA / 2;
				useLA();
				if (getH(point) > getH(point_i))
				{
					point = point_i;

				}
				else
				{
					LA = LA / -1;
					useLA();
					if (getH(point) > getH(point_i))
					{
						point = point_i;

					}
					else
					{
						LA *= 2;
						useLA();
						if (getH(point) > getH(point_i))
						{
							point = point_i;

						}
						else
						{
							LA /= 2;

						}
					}
				}
			}
		}
	}

}


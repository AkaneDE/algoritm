#pragma once
#include <vector>
#include "utils.h"
#include "filters.h"
#include "TDoA.h"
using namespace std;
using namespace filter;

class all
{
private:
	vector<float> v1;
	vector<float> v2;
	vector<float> v3;
	int d1 = 0;
	int d2 = 0;
	int d3 = 0;
	double x1; double y1; double z1;
	double x2; double y2; double z2;
	double x3; double y3; double z3;




	size_t n = 3;
	const float del_t = 1 / 40000;
	void get_Del()
	{
		float R = -1;
		int lag_max = max(v1.size(), max(v2.size(), v3.size()));
		
		for (int lag = 0; lag < lag_max; lag++) {
			if (corr(v1, v2, lag) > R)
			{
				R = corr(v1, v2, lag);
				d1 = lag;
			}
			if (corr(v2, v1, lag) > R)
			{
				R = corr(v2, v1, lag);
				d1 = -lag;
			}
		}

		R = -1;

		for (int lag = 0; lag < lag_max; lag++) {
			
			if (corr(v1, v3, lag) > R)
			{
				R = corr(v1, v3, lag);
				d2 = lag;
			}
			if (corr(v3, v1, lag) > R)
			{
				R = corr(v3, v1, lag);
				d2 = -lag;
			}
		}

		R = -1;
		
		for (int lag = 0; lag < lag_max; lag++) {
			if (corr(v2, v3, lag) > R)
			{
				R = corr(v2, v3, lag);
				d3 = lag;
			}
			if (corr(v3, v2, lag) > R)
			{
				R = corr(v3, v2, lag);
				d3 = -lag;
			}
		}

	}

	vector<float> use_filter(vector<float> v)
	{
		MedianFilter m = MedianFilter(3);
		for (int j = 0; j < v.size(); ++j)
		{
			v[j] = (m.filter(v[j]));

		}
		return v;
	}
public:
	all(vector<float> v1, vector<float> v2, vector<float> v3, double x1, double x2, double x3, double y1, double y2, double y3, double z1, double z2, double z3)
	{
		this->v1 = v1;
		this->v2 = v2;
		this->v3 = v3;
		this->x1 = x1; this->y1 = y1; this->z1 = z1;
		this->x2 = x2; this->y2 = y2; this->z2 = z2;
		this->x3 = x3; this->y3 = y3; this->z3 = z3;
	}
	void TDoA()
	{
		get_Del();
		v1 = use_filter(v1);
		v2 = use_filter(v2);
		v3 = use_filter(v3);
		TDoA_Chan t = TDoA_Chan(d1 * del_t, d2 * del_t, d3 * del_t, x1, x2, x3, y1, y2, y3, z1, z2, z3);
		t.DoTDoA();
	}


};


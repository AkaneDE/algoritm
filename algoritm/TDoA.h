#pragma once
#include <cmath>
#include <Eigen/Dense>
#include <vector>
#include <iostream>

using namespace Eigen;
using namespace std;

class TDoA_Chan
{
private:
    const double c = 3120.77;
    double t1_2, t2_3, t3_1, x1, x2, x3, y1, y2, y3, z1, z2, z3, d1, d2, d3;

public:
    TDoA_Chan(double t1, double t2, double t3, double x1, double x2, double x3, double y1, double y2, double y3, double z1, double z2, double z3)
    {
        this->x1 = x1;
        this->x2 = x2;
        this->x3 = x3;

        this->y1 = y1;
        this->y2 = y2;
        this->y3 = y3;

        this->z1 = z1;
        this->z2 = z2;
        this->z3 = z3;

        this->t1_2 = t1;
        this->t2_3 = t2;
        this->t3_1 = t3;




    }
    MatrixXd DoTDoA()
    {
        MatrixXd xyz(3, 1);
        xyz << 50, 50, 50;
        for (int i = 0; i < xyz.rows(); ++i) {
            for (int j = 0; j < xyz.cols(); ++j) {
                std::cout << "Element at (" << i << "," << j << "): " << xyz(i, j) << std::endl;
            }
        }
        do
        {
            MatrixXd Del;
            d1 = sqrt(pow(xyz(0, 0) - x1, 2) + pow(xyz(1, 0) - y1, 2) + pow(xyz(2, 0) - z1, 2));
            d2 = sqrt(pow(xyz(0, 0) - x2, 2) + pow(xyz(1, 0) - y2, 2) + pow(xyz(2, 0) - z2, 2));
            d3 = sqrt(pow(xyz(0, 0) - x3, 2) + pow(xyz(1, 0) - y3, 2) + pow(xyz(2, 0) - z3, 2));

            MatrixXd r(3, 1);
            r << (d3 - d2) - t2_3 * c,
                (d2 - d1) - t1_2 * c,
                (d1 - d3) - t3_1 * c;
            for (int i = 0; i < r.rows(); ++i) {
                for (int j = 0; j < r.cols(); ++j) {
                    std::cout << "Element at (" << i << "," << j << "): " << r(i, j) << std::endl;
                }
            }
            MatrixXd J(3, 3);
            J << (xyz(0, 0) - x3) / d3 - (xyz(0, 0) - x2) / d2, (xyz(1, 0) - y3) / d3 - (xyz(1, 0) - y2) / d2, (xyz(2, 0) - z3) / d3 - (xyz(2, 0) - z2) / d2,
                (xyz(0, 0) - x2) / d2 - (xyz(0, 0) - x1) / d1, (xyz(1, 0) - y2) / d2 - (xyz(1, 0) - y1) / d1, (xyz(2, 0) - z2) / d2 - (xyz(2, 0) - z1) / d1,
                (xyz(0, 0) - x1) / d1 - (xyz(0, 0) - x3) / d3, (xyz(1, 0) - y1) / d1 - (xyz(1, 0) - y3) / d3, (xyz(2, 0) - z1) / d1 - (xyz(2, 0) - z3) / d3;
            //cout << (xyz(0, 0) - x1) / d1 <<" "<< (xyz(0, 0) - x3) / d3 << endl;
            Del = -(J.transpose() * J).inverse() * J.transpose() * r;
            /*for (int i = 0; i < J.rows(); ++i) {
                for (int j = 0; j < J.cols(); ++j) {
                    std::cout << "Element at JJJJJJ(" << i << "," << j << "): " << J(i, j) << std::endl;
                }
            }*/

            /*for (int i = 0; i < Del.rows(); ++i) {
                for (int j = 0; j < Del.cols(); ++j) {
                    std::cout << "Element at DEL(" << i << "," << j << "): " << Del(i, j) << std::endl;
                }
            }*/
            xyz = xyz + Del;
             for (int i = 0; i < xyz.rows(); ++i) {
                 for (int j = 0; j < xyz.cols(); ++j) {
                     std::cout << "Element at (" << i << "," << j << "): " << xyz(i, j) << std::endl;
                 }
             }
            if (abs(Del(0, 0)) <= 0.001 && abs(Del(1, 0)) <= 0.001 && abs(Del(2, 0)) <= 0.001)
                break;
            std::cout << "!!!" << std::endl;

        } while (true);

        return xyz;

    }
};
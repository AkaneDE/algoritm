


#include <iostream>
#include <vector>
#include "all.h"  
#include <cstdlib>
#include <ctime>
using namespace std;
vector<float> v1;
vector<float> v2;
vector<float> v3;
//Частота дискр 20 кгц - частота дискр. 20 тыс отсч в сек.
//взять реальные значения на графике, выбрать точку, нарисовать отрезки, посчитать количество секунд.
//
void create_vector() {
    srand(time(0));
    int vectorLength = rand() % 10000 + 1000;

    for (int i = 0; i < vectorLength; ++i) {
        v1.push_back(1000-i);
    }
    //vectorLength = rand() % 500 + 10;
    for (int i = 0; i < vectorLength; ++i) {
        v2.push_back(v1[i]);
        v3.push_back(v1[i]);

    }
    
    for (int i = 0; i < 2; ++i) {
        v2.erase(v2.begin());
    }
    vectorLength = rand() % 500 + 10;
    for (int i = 0; i < 7; ++i) {
        v3.insert(v3.begin(), 1000-i);
        
    }

}



int main()
{
    create_vector();
    all a = all(v1, v2, v3, 0, 0, 0, -14, 0, 0, 15, 0, 0);
    a.TDoA();


}


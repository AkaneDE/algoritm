


#include <iostream>
#include <vector>
#include "all.h"  
#include <cstdlib>
#include <ctime>
using namespace std;
vector<float> v1;
vector<float> v2;
vector<float> v3;
void create_vector() {
    srand(time(0));
    int vectorLength = rand() % 100000 + 1000;

    for (int i = 0; i < vectorLength; ++i) {
        v1.push_back(static_cast<float>(rand()) / static_cast<float>(RAND_MAX));
    }
    vectorLength = rand() % 500 + 10;
    
    v2 = v1;
    v3 = v1;
    for (int i = 0; i < vectorLength; ++i) {
        v2.push_back(0);
    }
    vectorLength = rand() % 500 + 10;
    for (int i = 0; i < vectorLength; ++i) {
        v3.insert(v3.begin(), 0);
    }

}



int main()
{
    create_vector();
    all a = all(v1, v2, v3, 100, 0, 200, 100, 0, 100, 100, 0, 100);
    a.TDoA();


}


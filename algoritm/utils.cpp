//
// Created by Horton on 14.11.2023.
//

#include "utils.h"
float corr(std::vector<float>& v, std::vector<float>& f, int lag) {
    float R = 0;

    int N = std::max(f.size(), v.size());

    if (v.size() > f.size()) {
        int n = v.size() - f.size();

        for (int i = 0; i < n; i++) {
            f.push_back(0);
        }
    }
    else if (v.size() < f.size()) {
        int n = f.size() - v.size();

        for (int i = 0; i < n; i++) {
            v.push_back(0);
        }
    }

    float K = 1 / (float)(N - lag);

    for (int k = 0; k < N - lag; k++) {
        R += f[k] * v[k + lag];
    }

    R = R * K;
    return R;
}
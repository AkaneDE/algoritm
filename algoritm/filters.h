//
// Created by justnik on 28.08.2021.
//

#ifndef FILTERS_FILTERS_HPP
#define FILTERS_FILTERS_HPP

#include <vector>

namespace filter {

    class IFilter {
    public:
        virtual double filter(double new_val) = 0;
    };

    class MedianFilter : public IFilter {
    public:
        explicit MedianFilter(size_t size);

        double filter(double new_val) final;

    private:
        std::vector<double> data;
        std::size_t size;
        std::size_t ind = 0;
    };

    class RunningAvg : public IFilter {
    public:
        explicit RunningAvg(size_t size);

        double filter(double new_val) final;

    private:
        std::vector<double> data;
        std::size_t size;
        std::size_t ind = 0;
    };

    class ExpSmooth : public IFilter {
    public:
        explicit ExpSmooth(double coef);

        double filter(double new_val) final;

    private:
        double coef = 1;
        double last_val = std::numeric_limits<double>::quiet_NaN();
    };

    class ABFilter : public IFilter {
    public:
        ABFilter(double dt, double sigma_process, double sigma_noise);

        double filter(double new_val) final;

    private:
        double dt;
        double sigma_process;
        double sigma_noise;
    };

    class KalmanFilter : public IFilter {
    public:
        KalmanFilter(double dt, double errMeasure);
        // ABFilter(double dt, double sigma_process, double sigma_noise);

        double filter(double new_val) final;

    private:
        double dt;
        double err_measure;
    };

}
#endif //FILTERS_FILTERS_HPP

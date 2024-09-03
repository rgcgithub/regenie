#ifndef SURDATA_H
#define SURDATA_H

#include <iostream>
#include "Eigen/Dense"
#include "Eigen/Core"
#include "Eigen/Sparse"
#include <algorithm>
#include <cmath>
#include <map>
#include <vector>

class survival_data {
    public:
        Eigen::VectorXd status;
        ArrayXb keep_sample_order;

        // number of samples and variants
        unsigned int n, neff;
        // number of event
        unsigned int n_events;
        unsigned int n_unique_time;
        // time order
        SpMat permute_mtx;
        // ordered survival data;
        Eigen::VectorXd status_order;
        Eigen::VectorXd time_order_event;
        // tie
        std::map<double, std::vector<int>> unique_time_indices;
        Eigen::VectorXd time_first_index;
        std::map<double, Eigen::VectorXi> ties_index;
        // weights
        Eigen::VectorXd w, w_orig;
        Eigen::VectorXd ww, ww_k;
        Eigen::VectorXd dd;
        // Eigen::SparseMatrix<int, Eigen::RowMajor> R;
        SpMat R; // n by K
        // Eigen::MatrixXd R;
        Eigen::VectorXd rskcount;
        
        survival_data();
        void setup(const Eigen::VectorXd& event_time, const Eigen::VectorXd& event_status, const ArrayXb& mask, const bool& norm_weights = false);
    
    private:
        void _getOrder(const Eigen::VectorXd&, const Eigen::VectorXd&);
        void _findTies(const Eigen::VectorXd&, const Eigen::VectorXi&);
};

Eigen::VectorXd cumulativeSum(const Eigen::VectorXd& x);
Eigen::VectorXd cumulativeSum_reverse2(const Eigen::VectorXd& x);

#endif
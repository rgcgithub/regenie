#include "Regenie.hpp"
#include "survival_data.hpp"

using namespace Eigen;
using namespace std;

survival_data::survival_data(){};

void survival_data::setup(const Eigen::VectorXd& event_time, const Eigen::VectorXd& event_status, const ArrayXb& mask, const bool& norm_weights) {
    if (event_time.size() != event_status.size()) { throw std::runtime_error("event_time and event_status should have same length."); }
    status = mask.select(event_status, -999).matrix();
    n = event_time.size();
    neff = (mask == true).count();

    // order event times
    _getOrder(event_time, status);
    Eigen::VectorXd time_order = permute_mtx * mask.select(event_time, -999).matrix(); // missing or masked samples are in the front
    status_order = permute_mtx * status;
    Eigen::VectorXd keep_sample_double = mask.cast<double>();
    keep_sample_order = ((permute_mtx * keep_sample_double).array() != 0).cast<bool>();

    w = Eigen::VectorXd::Ones(n);
    if (norm_weights) w /= neff;
    w_orig = mask.select(w, 0);
    w = keep_sample_order.select(w, 0);

    // find ties
    // time with event
    n_events = (status_order.array() == 1).count();

    time_order_event.resize(n_events);
    Eigen::VectorXi time_order_event_index(n_events);
    int idx = 0;
    for (unsigned int i = 0; i < n; ++i) {
        if (status_order(i) == 1) {
            time_order_event(idx) = time_order(i);
            time_order_event_index(idx) = i;
            ++idx;
        }
    }
    _findTies(time_order_event, time_order_event_index);

    dd = keep_sample_order.select(status_order, 0);
    ww = w;
    double wsum;
    for (const auto& entry: ties_index) {
        const Eigen::VectorXi& ties = entry.second;
        if (ties.size() > 1) {
            for (int i = 0; i < ties.size(); i++) {
                dd(ties(i)) = 0;
                ww(ties(i)) = 0;
            }
            dd(ties(0)) = 1;

            if (norm_weights) {
                wsum = static_cast<double>(ties.size())/neff;
            } else {
                wsum = static_cast<double>(ties.size());
            }
            ww(ties(0)) = wsum;
        }
    }
    rskcount = cumulativeSum(dd);

	// R matrix
	n_unique_time = time_first_index.size();
	R.resize(n, n_unique_time);
    // R = Eigen::MatrixXd::Zero(n_unique_time, n);
	for (unsigned int k = 0; k < n_unique_time; ++k) {
		if (k < n_unique_time - 1) {
			for (unsigned int i = time_first_index[k]; i < time_first_index[k+1]; ++i){
				R.insert(i, k) = 1;
                // R(k,i) = 1;
			}
		} else {
			for (unsigned int i = time_first_index[k]; i < n; ++i){
				R.insert(i, k) = 1;
                // R(k,i) = 1;
			}
		}
	}
    R.makeCompressed();
	
	ww_k.resize(n_unique_time);
	int idx_t = 0;
	for(unsigned int i = 0; i < n; i++) {
    	if(dd(i) == 1) {
			ww_k(idx_t) = ww(i);
			idx_t += 1;
    	}
	}
}

void survival_data::_getOrder(const Eigen::VectorXd& time, const Eigen::VectorXd& status) {
    // order = Eigen::seq(0, n - 1);
    Eigen::VectorXi order(n);
    for (unsigned int i = 0; i < n; ++i) {
        order(i) = i;
    }
    
    std::sort(order.data(), order.data() + n, [&](int i, int j) {
        if (time(i) < time(j)) {
            return true;
        } else if (time(i) > time(j)) {
            return false;
        } else {
            return status(i) > status(j);
        }
    });

    permute_mtx.resize(n, n);
    permute_mtx.reserve(VectorXi::Constant(n,1));
    for (unsigned int i = 0; i < n; ++i) {
        permute_mtx.insert(i, order(i)) = 1.0;
    }
    permute_mtx.makeCompressed();
}

void survival_data::_findTies(const Eigen::VectorXd& x_sorted, const Eigen::VectorXi& index) {
    int n_times = x_sorted.size();
    // indices for each unique time
    for (int i = 0; i < n_times; ++i) {
        unique_time_indices[x_sorted(i)].push_back(index(i));
    }

    for (const auto& entry: unique_time_indices) {
        double time = entry.first;
        const vector<int>& indices = entry.second;

        time_first_index.conservativeResize(time_first_index.size() + 1);
        time_first_index(time_first_index.size() - 1) = indices[0];

        if (indices.size() > 1) {
            Eigen::VectorXi ties(indices.size());
            for (unsigned int i = 0; i < indices.size(); ++i) {
                ties(i) = indices[i];
            }
            ties_index[time] = ties;
        }
    }
}

Eigen::VectorXd cumulativeSum(const Eigen::VectorXd& x) {
    int n = x.size();
    Eigen::VectorXd x_cumsum(n);
    x_cumsum(0) = x(0);
    for (int i = 1; i < n; ++i) {
        x_cumsum(i) = x_cumsum(i - 1) + x(i);
    }
    return x_cumsum;
}

Eigen::VectorXd cumulativeSum_reverse2(const Eigen::VectorXd& x) {
    int n = x.size();
    Eigen::VectorXd x_cumsum(n);
    x_cumsum(n-1) = x(n-1);
    for (int i = n-2; i >= 0; --i) {
        x_cumsum(i) = x_cumsum(i + 1) + x(i);
    }
    return x_cumsum;
}
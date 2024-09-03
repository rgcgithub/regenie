#include "Regenie.hpp"
#include "survival_data.hpp"
#include "cox_ridge.hpp"
#include "cox_score.hpp"

using namespace Eigen;
using namespace std;
using namespace boost;

cox_mle::cox_mle(){}

void cox_mle::setup(const survival_data& survivalData, const Eigen::MatrixXd& Xmat, const Eigen::VectorXd& offset_val, const ArrayXb& mask, const int& max_iter, const int& max_inner_iter, const double& tolerance, const bool& verbose_obj, const Eigen::VectorXd& beta_init, const Eigen::VectorXd& eta_init) {
    converge = false;
    p = Xmat.cols();

	_niter = max_iter;
    _mxitnr = max_inner_iter;
    _tol = tolerance;
    _verbose = verbose_obj;

    if (beta_init.size() > 0) {
        beta = beta_init;
        eta = eta_init;
    } else {
        beta = Eigen::VectorXd::Zero(p);
        eta = mask.select(offset_val, 0).matrix();
    }
    eta_order = mask.select(survivalData.permute_mtx * eta, 0).matrix();
    lambda0.resize(survivalData.n);
    mu.resize(survivalData.n);
    Y.resize(survivalData.n);
    residual.resize(survivalData.n);
    loglike.resize(_niter + 1);
}

void cox_mle::fit(const survival_data& survivalData, const Eigen::MatrixXd& Xmat, const Eigen::VectorXd& offset_val, const ArrayXb& mask) {
    Eigen::VectorXd beta_old;
    Eigen::VectorXd sqrtWY, XtWY;
    int ii;
	compute_loglike(survivalData);
	loglike(0) = loglik_val;
    // std::cout << "start MLE fitting:\n";
    // std::cout << "beta: " << beta << "\n";
    // std::cout << "loglik_val: " << loglik_val << "\n";
    iter = 0;
    for (int t = 0; t < _niter; ++t) { 
        ii = 0;
        beta_old = beta;
        mu = survivalData.w_orig.array() * lambda0.array() * eta.array().exp();
        residual = survivalData.w_orig.array() * survivalData.status.array() - mu.array();
        Y = mask.select(eta - offset_val, 0).matrix() +
            (mu.array() != 0).select(residual.array()/mu.array(), 0).matrix();
      	
        // update beta
        if (p == 0) {
            sqrtWX.resize(survivalData.n, 0);
            XtWX.resize(0,0);
            converge = true;
            break;
        }
        
        ++iter;
        sqrtWY = Y.array() * (mu.array().sqrt());
		sqrtWX = Xmat.array().colwise() * (mu.array().sqrt());
        XtWX = sqrtWX.transpose() * sqrtWX;
		XtWY = sqrtWX.transpose() * sqrtWY;
		beta = XtWX.colPivHouseholderQr().solve(XtWY);
		eta = mask.select(Xmat * beta + offset_val, 0).matrix();
        eta_order = survivalData.permute_mtx * eta;
        compute_loglike(survivalData);
        // std::cout << "iter: " << iter << "\n";
        // std::cout << "beta: " << beta << "\n";
        // std::cout << "loglik_val: " << loglik_val << "\n";
        // std::cout << "diff loglik_val: " << loglik_val - loglike(iter - 1) << "\n";

        if ((loglike(iter - 1) - loglik_val) > _tol) { // step-halving
			// std::cout << "\nLoglikelihood decreases at iteration " << iter << ", start step-halving.\n";
            while ((loglike(iter - 1) - loglik_val) > _tol) {
                ++ii;
                if (ii > _mxitnr) {
                    std::cout << "Convergence issue, inner loop: cannot correct step size\n";
                    return;
                    // throw std::runtime_error("inner loop: cannot correct step size");
                }
                beta = (beta_old + beta)/2;
                eta = mask.select(Xmat * beta + offset_val, 0).matrix();
                eta_order = survivalData.permute_mtx * eta;
                compute_loglike(survivalData);
                if (_verbose) {
                    std::cout << "beta: " << beta << "\n";
                    std::cout << "Iteration " << iter << " Halved, Objective: " << loglik_val << "\n";
                }
            }
		}
		loglike(iter) = loglik_val;

        // std::cout << "beta: " << beta << "\n";
        // std::cout << "loglik_val: " << loglik_val << "\n";
        // std::cout << "first_der_1: " << first_der_1 << "\n";
        // std::cout << "second_der_1: " << second_der_1 << "\n";
        // std::cout << "loglike val change: " << loglike(iter) - loglike(iter - 1) << "\n";
        // std::cout << "beta change max: " << ((beta.array() - beta_old.array()).abs()/(beta.array().abs() + beta_old.array().abs() + _tol)).maxCoeff() << "\n";

        if( loglike(iter) - loglike(iter - 1) < _tol || (ii <= 1 && ((beta.array() - beta_old.array()).abs()/(beta.array().abs() + beta_old.array().abs() + _tol)).maxCoeff() < _tol )) {
			mu = survivalData.w_orig.array() * lambda0.array() * eta.array().exp();
            residual = survivalData.w_orig.array() * survivalData.status.array() - mu.array();
        	Y = mask.select(eta - offset_val, 0).matrix() +
            	(mu.array() != 0).select(residual.array()/mu.array(), 0).matrix();
			sqrtWX = Xmat.array().colwise() * mu.array().sqrt();
            XtWX = sqrtWX.transpose() * sqrtWX;
            converge = true;
            break;
      	
        }
    }
    // std::cout << "finish fitting\n";
}

void cox_mle::compute_loglike(const survival_data& survivalData) {
    Eigen::VectorXd w_exp_eta = survivalData.w.array() * eta_order.array().exp();

    Eigen::VectorXd S0 = cumulativeSum_reverse2(survivalData.R.transpose() * w_exp_eta); // length K, risk set sum at each unique failure time
    double log_terms_sum = (survivalData.ww_k.array() * S0.array().log()).sum();

    loglik_val = (survivalData.w.array() * eta_order.array() * (survivalData.status_order.array() == 1).cast<double>()).sum() - log_terms_sum;

    Eigen::VectorXd ww_rsk = cumulativeSum(survivalData.ww_k.array() / S0.array()); // length K
    for (unsigned int i = 0; i < survivalData.n; ++i) {
        if (survivalData.rskcount(i) == 0) {
            lambda0(i) = 0;
        } else {
            lambda0(i) = ww_rsk(int(survivalData.rskcount(i)) - 1);
        }
    }
    lambda0 = survivalData.permute_mtx.transpose() * lambda0;
}

void cox_mle::cox_test_prep(const survival_data& survivalData, const Eigen::MatrixXd& Xmat, const Eigen::VectorXd& offset_val, const ArrayXb& mask) {
    double eta_mean = eta_order.array().sum()/eta_order.size();
    w_exp_eta = survivalData.w.array() * (eta_order.array() - eta_mean).exp();

    Eigen::VectorXd rskden = cumulativeSum_reverse2(survivalData.R.transpose() * w_exp_eta); // length K

    Dhalf = (survivalData.ww_k.array().sqrt()) / rskden.array();

    Eigen::MatrixXd Gamma_X = (survivalData.permute_mtx * Xmat).array().colwise() * w_exp_eta.array(); // n by p
    UhalfX.resize(survivalData.n_unique_time, p);

    Eigen::VectorXd RGammaXr;
    for (int r = 0; r < p; ++r) {
        RGammaXr = cumulativeSum_reverse2(survivalData.R.transpose() * Gamma_X.col(r)); // length K
        UhalfX.col(r) = Dhalf.array() * RGammaXr.array();
    }

    cov_inv.resize(p, p);
    if (p > 0) {
        cov_inv = (XtWX - UhalfX.transpose() * UhalfX).colPivHouseholderQr().inverse();
    }

    Eigen::MatrixXd X1 = MatrixXd::Zero(survivalData.n, p + 1);
    X1 << Eigen::VectorXd::Ones(survivalData.n), Xmat;
    WX1 = X1.array().colwise() * mu.array();
    Eigen::MatrixXd X1tWX1 = X1.transpose() * WX1;
    X1_X1WX1inv = X1 * X1tWX1.colPivHouseholderQr().inverse();
    
    double res_mean = residual.mean();
    // Compute the variance
    res_var = (residual.array() - res_mean).square().sum() / (residual.size() - 1);
}

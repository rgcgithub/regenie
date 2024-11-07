#include "Regenie.hpp"
#include "survival_data.hpp"
#include "cox_firth.hpp"

using namespace Eigen;
using namespace std;

cox_firth::cox_firth(){}

void cox_firth::setup(const survival_data& survivalData, const Eigen::MatrixXd& Xmat, const Eigen::VectorXd& offset_val, const int& cols_incl, const int& max_iter, const int& max_inner_iter, const double& tolerance, const double& stephalf_tol, const double& beta_tol, const double& max_step, const bool& use_firth, const bool& verbose_obj, const Eigen::VectorXd& beta_init) {
	converge = false;
    p = Xmat.cols();
    
	_niter = max_iter;
    _mxitnr = max_inner_iter;
    _tol = tolerance;
    _stephalf_tol = stephalf_tol;
    _betatol = beta_tol;
	_maxstep = max_step;
	_usefirth = use_firth;
    _verbose = verbose_obj;
    _cols_incl = cols_incl;
    loglike.resize(_niter + 1);
    first_der.resize(p);
    mu.resize(survivalData.n);
    residual.resize(survivalData.n);

    beta = Eigen::VectorXd::Zero(p);
    if (beta_init.size() > 0) {
        beta.head(_cols_incl) = beta_init.head(_cols_incl);
        eta = Xmat * beta + offset_val;
    } else {
        eta = offset_val;
    }
    eta_order = survivalData.keep_sample_order.select(survivalData.permute_mtx * eta, 0).matrix();
    if (p == 0) {
        _usefirth = false;
    }
}

void cox_firth::cox_firth_likelihood(const survival_data& survivalData, const Eigen::MatrixXd& Xmat) {
    Eigen::VectorXd w_exp_eta, ww_rsk, S0;
    Eigen::VectorXd lambda0(survivalData.n);
    Eigen::MatrixXd S1, GammaX, XtW;
    Eigen::MatrixXd S2 = Eigen::MatrixXd::Zero(p, p);
    std::vector<Eigen::MatrixXd> firth_der;
    double log_terms_sum;
    second_der = Eigen::MatrixXd::Zero(p, p);

    if (_usefirth) {
        firth_der.resize(p);
        for(int i = 0; i < p; i++) {
            firth_der[i] = Eigen::MatrixXd::Zero(p, p);
        }
    }

    exp_eta = eta_order.array().exp();
    w_exp_eta = survivalData.w.array() * exp_eta.array();

    S0 = cumulativeSum_reverse2(survivalData.R.transpose() * w_exp_eta); // length K, risk set sum at each unique failure time
    log_terms_sum = (survivalData.ww_k.array() * S0.array().log()).sum();

    loglik_val = (survivalData.w.array() * eta_order.array() * (survivalData.status_order.array() == 1).cast<double>()).sum() - log_terms_sum;

    // double mean_eta = (eta.array() * survivalData.w_orig.array()).sum()/survivalData.w_orig.array().sum();
    // Eigen::VectorXd eta_center = eta_order.array() - mean_eta;
    // exp_eta = eta_center.array().exp();
    // w_exp_eta = survivalData.w.array() * exp_eta.array();
    // S0 = cumulativeSum_reverse2(survivalData.R.transpose() * w_exp_eta); // length K, risk 

    ww_rsk = cumulativeSum(survivalData.ww_k.array() / S0.array());
    for (unsigned int i = 0; i < survivalData.n; ++i) {
        if (survivalData.rskcount(i) == 0) {
            lambda0(i) = 0;
        } else {
            lambda0(i) = ww_rsk(int(survivalData.rskcount(i)) - 1);
        }
    }
    mu = lambda0.array() * w_exp_eta.array();

    S1 = survivalData.R.transpose() * ((survivalData.permute_mtx * Xmat).array().colwise() * w_exp_eta.array()).matrix(); // K by p

    GammaX = (survivalData.permute_mtx * Xmat).array().colwise() * w_exp_eta.array().sqrt(); // n by p
    for (int k = survivalData.n_unique_time - 1; k >= 0; --k) {
        if (k < survivalData.n_unique_time - 1) {
            S1.row(k) += S1.row(k+1);
        }

        std::vector<int> k_indices;
        for (SpMat::InnerIterator it(survivalData.R, k); it; ++it) {
            k_indices.push_back(it.index());
        }
        
        // for (int j = 0; j < survivalData.R.cols(); ++j) {
        //     if (survivalData.R(k, j) != 0) {
        //         k_indices.push_back(j);
        //     }
        // }

        S2 += GammaX(k_indices, all).transpose() * GammaX(k_indices, all);

        second_der = second_der + survivalData.ww_k(k) * (S2/S0(k) - S1.row(k).transpose() * S1.row(k)/(std::pow(S0(k), 2)));
        if (_usefirth) {
            for (int t = 0; t < p; ++t) {
                firth_der[t] += survivalData.ww_k(k) * ((-S2 * S1(k,t) - S2.col(t) * S1.row(k) - S2.row(t).transpose() * S1.row(k))/(std::pow(S0(k), 2)) + 2 * S1.row(k).transpose() * S1.row(k) * S1(k,t)/(std::pow(S0(k), 3)));
            }
        }
    }
    if (p > 0) qrsd.compute(second_der);
    residual = survivalData.w.array() * (survivalData.status_order - mu).array();
    if(_cols_incl < p) {
        qrsd_incl.compute(second_der.block(0,0,_cols_incl,_cols_incl)); // p-1 by p-1
        if (_usefirth) {
            loglik_val += 0.5 * qrsd.logAbsDeterminant();
            XtW = ((survivalData.permute_mtx * Xmat.leftCols(_cols_incl)).array().colwise() * mu.array().sqrt()).transpose(); // p-1 by n
            first_der = (survivalData.permute_mtx * Xmat.leftCols(_cols_incl)).transpose() * survivalData.keep_sample_order.select(residual + 0.5 * (qrsd_incl.solve(XtW).array() * XtW.array()).colwise().sum().matrix().transpose(), 0); // qrsd.solve(XtW) is p-1 by n
            for (int t = 0; t < _cols_incl; ++t) {
                first_der(t) = first_der(t) + 0.5 * qrsd_incl.solve(firth_der[t].block(0,0,_cols_incl,_cols_incl)).trace();
            }
        } else {
            first_der = (survivalData.permute_mtx * Xmat.leftCols(_cols_incl)).transpose() * residual;
        }
    } else {
        if (_usefirth) {
            loglik_val += 0.5 * qrsd.logAbsDeterminant();
            XtW = ((survivalData.permute_mtx * Xmat).array().colwise() * mu.array().sqrt()).transpose(); // p by n
            first_der = (survivalData.permute_mtx * Xmat).transpose() * survivalData.keep_sample_order.select(residual + 0.5 * (qrsd.solve(XtW).array() * XtW.array()).colwise().sum().matrix().transpose(), 0); // qrsd.solve(XtW) is p by n
            for (int t = 0; t < p; ++t) {
                first_der(t) = first_der(t) + 0.5 * qrsd.solve(firth_der[t]).trace();
            }
        } else {
            first_der = (survivalData.permute_mtx * Xmat).transpose() * residual;
        }
    }
}

void cox_firth::fit(const survival_data& survivalData, const Eigen::MatrixXd& Xmat, const Eigen::VectorXd& offset_val) {
    Eigen::VectorXd steps, betanew;
    int ii;
    cox_firth_likelihood(survivalData, Xmat);
    loglike(0) = loglik_val;
    // std::cout << "start fitting:\n";
    // std::cout << "beta: " << beta << "\n";
    // std::cout << "loglik_val: " << loglik_val << "\n";
    // std::cout << "first_der: " << first_der << "\n";
    // std::cout << "second_der: " << second_der << "\n";
    iter = 0;
    if (p == 0 || _cols_incl == 0) {
        converge = true;
        residual = survivalData.permute_mtx.transpose() * residual;
        loglike.conservativeResize(iter+1);
        return;
    }
    betanew = beta;
    while (iter++ < _niter) {
        // std::cout << "iter: " << iter << "\n";
        ii = 0;
        if (_cols_incl < p) {
            steps = qrsd_incl.solve(first_der);
        } else{
            steps = qrsd.solve(first_der);
        }
        // std::cout << "steps: " << steps << "\n";
        for (int i = 0; i < steps.size(); ++i) {
            if (abs(steps(i)) >= _maxstep) {
                steps(i) = (steps(i) / fabs(steps(i))) * _maxstep;
            }
        }
        // std::cout << "adjusted steps: " << steps << "\n";
        betanew.head(_cols_incl) = beta.head(_cols_incl) + steps;
        // std::cout << "beta: " << betanew << "\n";
        eta = Xmat * betanew + offset_val;
        eta_order = survivalData.keep_sample_order.select(survivalData.permute_mtx * eta, 0);
        cox_firth_likelihood(survivalData, Xmat);
        // std::cout << "loglik_val: " << loglik_val << "\n";
        // std::cout << "diff loglik_val: " << loglik_val - loglike(iter - 1) << "\n";
        if ((loglike(iter - 1) - loglik_val) > _stephalf_tol) { // step-halving
            // std::cout << "\nLoglikelihood decreases at iteration " << iter << ", start step-halving.\n";
            ii = 0;
            while ((loglike(iter - 1) - loglik_val) > _stephalf_tol) {
                ++ii;
                // std::cout << "inner iteration: " << ii << "\n";
                if (ii > _mxitnr) {
                    // std::cout << "Convergence issue, inner loop: cannot correct step size, add eps\n";
                    steps.array() += 1e-6;
                    betanew.head(_cols_incl) = beta.head(_cols_incl) + steps;
                    eta = Xmat * betanew + offset_val;
                    eta_order = survivalData.keep_sample_order.select(survivalData.permute_mtx * eta, 0);
                    cox_firth_likelihood(survivalData, Xmat);
                    break;
                    // throw std::runtime_error("inner loop: cannot correct step size");
                }
                betanew = (beta + betanew)/2;
                eta = Xmat * betanew + offset_val;
                eta_order = survivalData.keep_sample_order.select(survivalData.permute_mtx * eta, 0);
                cox_firth_likelihood(survivalData, Xmat);
                if (_verbose) {
                    std::cout << "beta: " << betanew << "\n";
                    std::cout << "Iteration " << iter << " Halved, Objective: " << loglik_val << "\n";
                }
            }
        }
        loglike(iter) = loglik_val;
        // std::cout << "beta: " << betanew << "\n";
        // std::cout << "loglik_val: " << loglik_val << "\n";
        // std::cout << "loglik_val change: " << loglik_val - loglike(iter - 1) << "\n";
        // std::cout << "first_der max: " << first_der.array().abs().maxCoeff() << "\n";
        // std::cout << "beta change max: " << (beta - betanew).array().abs().maxCoeff() << "\n";
        if( first_der.array().abs().maxCoeff() < _tol || (ii <= 1 && (beta - betanew).array().abs().maxCoeff() < _betatol) ) {
            beta = betanew;
            converge = true;
            break;
        }
        beta = betanew;
    }
    residual = survivalData.permute_mtx.transpose() * residual;
    loglike.conservativeResize(iter+1);
    // std::cout << "finish fitting\n";
}


void cox_firth::cox_firth_likelihood_1(const survival_data& survivalData, const Eigen::VectorXd& g) {
    Eigen::VectorXd w_exp_eta, ww_rsk;
    Eigen::VectorXd lambda0(survivalData.n);
    Eigen::VectorXd S0, S1, S2, S3;
    double log_terms_sum;

    exp_eta = eta_order.array().exp();
    w_exp_eta = survivalData.w.array() * exp_eta.array();

    S0 = cumulativeSum_reverse2(survivalData.R.transpose() * w_exp_eta); // length K, risk set sum at each unique failure time
    log_terms_sum = (survivalData.ww_k.array() * S0.array().log()).sum();

    loglik_val = (survivalData.w.array() * eta_order.array() * (survivalData.status_order.array() == 1).cast<double>()).sum() - log_terms_sum;

    ww_rsk = cumulativeSum(survivalData.ww_k.array() / S0.array());
    for (unsigned int i = 0; i < survivalData.n; ++i) {
        if (survivalData.rskcount(i) == 0) {
            lambda0(i) = 0;
        } else {
            lambda0(i) = ww_rsk(int(survivalData.rskcount(i)) - 1);
        }
    }
    mu = lambda0.array() * w_exp_eta.array();

    S1 = cumulativeSum_reverse2(survivalData.R.transpose() * ((survivalData.permute_mtx * g).array() * w_exp_eta.array()).matrix()); // K by 1

    S2 = cumulativeSum_reverse2(survivalData.R.transpose() * ((survivalData.permute_mtx * g.array().pow(2).matrix()).array() * w_exp_eta.array()).matrix()); // K by 1
    
    second_der_1 = (survivalData.ww_k.array() * (S2.array()/S0.array() - S1.array().pow(2)/S0.array().pow(2))).sum();

    residual = survivalData.w.array() * (survivalData.status_order - mu).array();

    if (_usefirth) {
        loglik_val += 0.5 * log(fabs(second_der_1));
        
        S3 = cumulativeSum_reverse2(survivalData.R.transpose() * ((survivalData.permute_mtx * g.array().pow(3).matrix()).array() * w_exp_eta.array()).matrix());
        
        first_der_1 = (survivalData.permute_mtx * g).dot(residual) + 0.5 * (survivalData.ww_k.array() * (S3.array()/S0.array() - 3 * S2.array() * S1.array()/S0.array().pow(2) + 2 * S1.array().pow(3)/S0.array().pow(3))).sum()/second_der_1;
    } else {
        first_der_1 = (survivalData.permute_mtx * g).dot(residual);
    }
}

void cox_firth::fit_1(const survival_data& survivalData, const Eigen::VectorXd& g, const Eigen::VectorXd& offset_val) {
    Eigen::VectorXd betanew;
    double steps;
    int ii = 0;
    cox_firth_likelihood_1(survivalData, g);
    // std::cout << "start fitting:\n";
    // std::cout << "beta: " << beta << "\n";
    // std::cout << "loglik_val: " << loglik_val << "\n";
    // std::cout << "first_der_1: " << first_der_1 << "\n";
    // std::cout << "second_der_1: " << second_der_1 << "\n";
    loglike(0) = loglik_val;
    iter = 0;
    while (iter++ < _niter) {
        // std::cout << "iter: " << iter << "\n";
        steps = first_der_1/second_der_1;
        // std::cout << "first der: " << first_der_1 << "\n";
        // std::cout << "second der: " << second_der_1 << "\n";
        // std::cout << "steps: " << steps << "\n";
        if (abs(steps) >= _maxstep) {
            steps = (steps / fabs(steps)) * _maxstep;
        }
        // std::cout << "adjusted steps: " << steps << "\n";
        betanew = beta.array() + steps;
        eta = g * betanew + offset_val;
        eta_order = survivalData.keep_sample_order.select(survivalData.permute_mtx * eta, 0);
        cox_firth_likelihood_1(survivalData, g);
        // std::cout << "beta: " << betanew << "\n";
        // std::cout << "loglik_val: " << loglik_val << "\n";
        // std::cout << "diff loglik_val: " << loglik_val - loglike(iter - 1) << "\n";
        
        if ((loglike(iter - 1) - loglik_val) > _stephalf_tol) { // step-halving
            // std::cout << "\nLoglikelihood decreases at iteration " << iter << ", start step-halving.\n";
            ii = 0;
            while ((loglike(iter - 1) - loglik_val) > _stephalf_tol) {
                ++ii;
                // std::cout << "inner iteration: " << ii << "\n";
                if (ii > _mxitnr) {
                    // std::cout << "Convergence issue, inner loop: cannot correct step size, add eps\n";
                    steps += 1e-6;
                    betanew = beta.array() + steps;
                    eta = g * betanew + offset_val;
                    eta_order = survivalData.keep_sample_order.select(survivalData.permute_mtx * eta, 0);
                    cox_firth_likelihood_1(survivalData, g);
                    break;
                    // throw std::runtime_error("inner loop: cannot correct step size");
                }
                betanew = (beta + betanew)/2;
                eta = g * betanew + offset_val;
                eta_order = survivalData.keep_sample_order.select(survivalData.permute_mtx * eta, 0);
                cox_firth_likelihood_1(survivalData, g);
                if (_verbose) {
                    std::cout << "beta: " << betanew << "\n";
                    std::cout << "Iteration " << iter << " Halved, Objective: " << loglik_val << "\n";
                }
            }
        }
        loglike(iter) = loglik_val;
        // std::cout << "beta: " << betanew << "\n";
        // std::cout << "loglik_val: " << loglik_val << "\n";
        // std::cout << "first_der_1: " << first_der_1 << "\n";
        // std::cout << "second_der_1: " << second_der_1 << "\n";
        // std::cout << "first_der max: " << fabs(first_der_1) << "\n";
        // std::cout << "beta change max: " << (beta - betanew).array().abs().maxCoeff() << "\n";
        if (fabs(first_der_1) < _tol || (ii <= 1 && (beta - betanew).array().abs().maxCoeff() < _betatol)) {
            beta = betanew;
            converge = true;
            break;
        }
        beta = betanew;
    }
    residual = survivalData.permute_mtx.transpose() * residual;
    loglike.conservativeResize(iter+1);
    // std::cout << "finish fitting\n";
}
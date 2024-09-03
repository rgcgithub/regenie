#include "Regenie.hpp"
#include "survival_data.hpp"
#include "cox_ridge.hpp"

using namespace Eigen;
using namespace std;

cox_ridge::cox_ridge(const survival_data& survivalData, const Eigen::MatrixXd& Xmat, const Eigen::VectorXd& offset_val, const ArrayXb& mask, const double& lambda_val, const int& max_iter, const int& max_inner_iter, const double& tolerance, const bool& verbose_obj, const Eigen::VectorXd& beta_init, const double& null_deviance) {
    converge = false;

    if (beta_init.size() > 0) {
        beta = beta_init;
    } else {
        beta = Eigen::VectorXd::Zero(Xmat.cols());
    }
    lambda = lambda_val;
    _niter = max_iter;
    _mxitnr = max_inner_iter;
    _tol = tolerance;
    _verbose = verbose_obj;
    _object.resize(_niter + 1);
    _deviance.resize(_niter + 1);

    eta = mask.select(Xmat * beta + offset_val, 0).matrix();
    eta_order = survivalData.permute_mtx * eta;
    
    if (null_deviance == -999) {
        _deviance(0) = _coxDeviance(survivalData);
    } else {
        _deviance(0) = null_deviance;
    }
    _object(0) = _deviance(0) + lambda * (beta.array().pow(2).sum())/2;
}

void cox_ridge::reset(const survival_data& survivalData, const Eigen::MatrixXd& Xmat, const Eigen::VectorXd& offset_val, const ArrayXb& mask, const double& lambda_val, const Eigen::VectorXd& beta_init, const double& null_deviance) {
    // Reset the object's state based on the provided parameters
    converge = false;

    if (beta_init.size() > 0) {
        beta = beta_init;
    } else {
        beta = Eigen::VectorXd::Zero(Xmat.cols());
    }
    lambda = lambda_val;

    // Calculate eta, eta_order, or other necessary calculations
    eta = mask.select(Xmat * beta + offset_val, 0).matrix();
    eta_order = survivalData.permute_mtx * eta;

    _deviance.resize(_niter + 1);
    _object.resize(_niter + 1);
    if (null_deviance == -999){
        _deviance(0) = _coxDeviance(survivalData);
    } else {
        _deviance(0) = null_deviance;
    }
    _object(0) = _deviance(0) + lambda * (beta.array().pow(2).sum())/2;
}

void cox_ridge::coxGrad(const survival_data& survivalData) {
    double mean_eta = (eta.array() * survivalData.w_orig.array()).sum()/survivalData.w_orig.array().sum();
    Eigen::VectorXd eta_center = eta_order.array() - mean_eta;
    Eigen::VectorXd exp_eta = eta_center.array().exp();
    // Eigen::VectorXd exp_eta = eta_order.array().exp();
    Eigen::VectorXd rskden = cumulativeSum_reverse2(survivalData.w.array() * exp_eta.array());

    Eigen::VectorXd ww_rsk = survivalData.ww.array() / rskden.array();
    Eigen::VectorXd ww_rsk2 = survivalData.ww.array() / (rskden.array().pow(2));

    Eigen::VectorXd rskdeninv_n = cumulativeSum(survivalData.dd.array().cast<bool>().select(ww_rsk, 0));
    Eigen::VectorXd rskdeninv2_n = cumulativeSum(survivalData.dd.array().cast<bool>().select(ww_rsk2, 0));

    Eigen::VectorXd gradient_order = survivalData.w.array() * (survivalData.status_order.array() - exp_eta.array() * rskdeninv_n.array());

    Eigen::VectorXd diag_hessian_order = (survivalData.w.array() * exp_eta.array()).pow(2) * rskdeninv2_n.array() - 
        survivalData.w.array() * exp_eta.array() * rskdeninv_n.array();
    // change to original order
    _gradient = survivalData.permute_mtx.transpose() * gradient_order;
    _diagHessian = survivalData.permute_mtx.transpose() * diag_hessian_order;
}

double cox_ridge::_coxLoglik(const survival_data& survivalData) {
    Eigen::VectorXd rsk = cumulativeSum_reverse2(survivalData.w.array() * eta_order.array().exp());

    // take just the terms related to actual death times
    double log_terms_sum = (survivalData.ww.array() * (survivalData.keep_sample_order.select(rsk.array().log(), 0)) * (survivalData.dd.array() == 1).cast<double>()).sum();
    double loglik_val = (survivalData.w.array() * eta_order.array() * (survivalData.status_order.array() == 1).cast<double>()).sum() - log_terms_sum;
    return loglik_val;
}

double cox_ridge::_coxDeviance(const survival_data& survivalData) {
    Eigen::VectorXd w_sub;
    if (survivalData.unique_time_indices.size() == survivalData.n_events) {
        // no tie
        w_sub = Eigen::VectorXd::Ones(survivalData.n_events);
        w_sub /= survivalData.neff;
    } else {
        // tie
        w_sub.resize(survivalData.unique_time_indices.size());
        int idx = 0;
        for (const auto& entry: survivalData.unique_time_indices) {
            const vector<int>& ties = entry.second;
            w_sub(idx) = static_cast<double>(ties.size())/survivalData.neff;
            ++idx;
        }
    }
    double lsat = -(w_sub.array() * (w_sub.array().log())).sum();
    double loglik_val = _coxLoglik(survivalData);

    double deviance = 2 * (lsat - loglik_val);
    return deviance;
}

void cox_ridge::fit(const survival_data& survivalData, const Eigen::MatrixXd& Xmat, const Eigen::VectorXd& offset_val, const ArrayXb& mask) {
    int p = Xmat.cols();
    int ii = 0;
    Eigen::VectorXd beta_old;
    int break_pt = 1;
    for (int t = 1; t < _niter + 1; ++t) {
        ++break_pt;
        beta_old = beta;
        coxGrad(survivalData);
        Eigen::VectorXd z(survivalData.n);
        z = (_diagHessian.array() != 0).select(_gradient.array()/_diagHessian.array(), 0).matrix();
        z = mask.select(eta - offset_val, 0) - z;
        for (unsigned int k = 0; k < p; ++k) {
            Eigen::VectorXd r = _diagHessian.array() * (z - eta + offset_val).array();
            eta = eta - mask.select(Xmat.col(k) * beta(k), 0).matrix();
            beta(k) = (r.dot(Xmat.col(k)) + beta(k) * (Xmat.col(k).array().pow(2) * _diagHessian.array()).sum()) /
                ((Xmat.col(k).array().pow(2) * _diagHessian.array()).sum() - lambda);
            eta = eta + mask.select(Xmat.col(k) * beta(k), 0).matrix();
        }
        eta_order = survivalData.permute_mtx * eta;
        _deviance(t) = _coxDeviance(survivalData);
        _object(t) = _deviance(t) + lambda * (beta.array().pow(2).sum())/2;
        if (_verbose) {
            std::cout << "Iteration " << t << " objective: " << _object(t) << "; diff: " << _object(t) - _object(t-1) << "; rel diff: " << abs(_object(t) - _object(t - 1)) / (0.1 + abs(_object(t))) << "; score: " << (_gradient.transpose() * Xmat - lambda * beta.transpose()).cwiseAbs().maxCoeff() << "\n";
        }

        if ( (_deviance(t) - _deviance(t-1)) > _tol ) {
            std::cout << "\nDeviance increases at iteration " << t << ".\n";
            ii = 0;
            while ( (_deviance(t) - _deviance(t-1)) > _tol ) {
                ++ii;
                if (ii > _mxitnr) {
                    std::cout << "Convergence issue, inner loop: cannot correct step size\n";
                    return;
                    // throw std::runtime_error("inner loop: cannot correct step size");
                }
                beta = (beta + beta_old)/2;
                eta = mask.select(Xmat * beta + offset_val, 0).matrix();
                eta_order = survivalData.permute_mtx * eta;
                _deviance(t) = _coxDeviance(survivalData);
                _object(t) = _deviance(t) + lambda * (beta.array().pow(2).sum())/2;
                if (_verbose) {
                    std::cout << "Iteration " << t << " Halved, Objective: " << _object(t) << "; diff: " << _object(t) - _object(t-1) << ".\n";
                }
            }
        }

        if (abs(_object(t) - _object(t - 1)) / (0.1 + abs(_object(t))) < _tol || (_gradient.transpose() * Xmat - lambda * beta.transpose()).cwiseAbs().maxCoeff() < _tol ) {
            converge = true;
            break;
        }
    }

    if (break_pt < (_niter + 1)) {
        _deviance.conservativeResize(break_pt);
        _object.conservativeResize(break_pt);
    }
    dev_ratio = 1 - _deviance(_deviance.size() - 1) / _deviance(0);
}

Eigen::VectorXd cox_ridge::get_gradient() {
    return _gradient;
}

double cox_ridge::get_deviance() {
    return _deviance(_deviance.size() - 1);
}

double cox_ridge::get_null_deviance() {
    return _deviance(0);
}

double cox_ridge::get_object() {
    return _object(_object.size() - 1);
}

Eigen::VectorXd cox_ridge::get_object_all() {
    return _object;
}

Eigen::VectorXd cox_ridge::get_deviance_all() {
    return _deviance;
}


cox_ridge_path::cox_ridge_path(const survival_data& survivalData, const Eigen::MatrixXd& Xmat, const Eigen::VectorXd& offset_val, const ArrayXb& mask, const int& nlambda, const double& lambda_min_max_ratio, const Eigen::VectorXd& lambda, const int& max_iter, const int& max_inner_iter, const double& tolerance, const bool& verbose_fit) {
    int p = Xmat.cols();
    // set lambda_vec
    if (lambda.size() > 0) {
        _user_define_lambda = true;
        _lambda_len = lambda.size();
        if (lambda.minCoeff() < 0) { throw std::runtime_error("lambda must >= 0."); }
        lambda_vec = lambda;
        std::sort(lambda_vec.data(), lambda_vec.data() + lambda_vec.size(), std::greater<double>());
    } else {
        double lambda_min_ratio;
        if (lambda_min_max_ratio >= 1) { 
            throw std::runtime_error("lambda_min_max_ratio should be less than 1."); 
        } else if (lambda_min_max_ratio == -1) {
            if (survivalData.neff < p) {
                lambda_min_ratio = 1e-2;
            } else {
                lambda_min_ratio = 1e-4;
            }
        } else {
            lambda_min_ratio = lambda_min_max_ratio;
        }
        _lambda_len = nlambda;
        cox_ridge coxRidge_null_lamb0(survivalData, Xmat, offset_val, mask, 0, max_iter, max_inner_iter, tolerance);
        coxRidge_null_lamb0.coxGrad(survivalData);
        Eigen::VectorXd gradient = coxRidge_null_lamb0.get_gradient();
        double lambda_max = _getCoxLambdaMax(Xmat, gradient);
        // lambda_vec = (Eigen::seq(0, _lambda_len - 1) * log(lambda_min_ratio) + log(lambda_max)).exp();
        Eigen::VectorXd index(nlambda);
        for (int i = 0; i < nlambda; ++i) {
            if (i > 0) {
                index(i) = static_cast<double>(i)/(nlambda - 1);
            } else {
                index(i) = i;
            }
        }
        lambda_vec = (index.array() * log(lambda_min_ratio) + log(lambda_max)).exp();
    }
    beta_mx.resize(p, _lambda_len);
    eta_mx.resize(survivalData.n, _lambda_len);
    object_val.resize(_lambda_len);
    dev_ratio.resize(_lambda_len);
    deviance.resize(_lambda_len);
    converge.resize(_lambda_len);
    niter = max_iter;
    mxitnr = max_inner_iter;
    tol = tolerance;
    verbose = verbose_fit;
}

double cox_ridge_path::_getCoxLambdaMax(const Eigen::MatrixXd& Xmat, const Eigen::VectorXd& gradient) {
    Eigen::VectorXd g = (Xmat.transpose() * gradient).array().abs();
    return g.maxCoeff() / 1e-3;
}

void cox_ridge_path::fit(const survival_data& survivalData, const Eigen::MatrixXd& Xmat, const Eigen::VectorXd& offset_val, const ArrayXb& mask) {
    int break_pt = 0;
    double cur_lambda = lambda_vec(0);
    cox_ridge coxRidge(survivalData, Xmat, offset_val, mask, cur_lambda, niter, mxitnr, tol, verbose);
    double nulldev_old = -999;
    Eigen::VectorXd beta_old(Xmat.cols());

    for (int k = 0; k < _lambda_len; ++k) {
        ++break_pt;
        if (k > 0) {
            cur_lambda = lambda_vec(k);
            coxRidge.reset(survivalData, Xmat, offset_val, mask, cur_lambda, beta_old, nulldev_old);
        }
        if (verbose) {
            std::cout << "lambda: " << cur_lambda << "\n";
        }
        coxRidge.fit(survivalData, Xmat, offset_val, mask);
        std::cout << "converge: " << coxRidge.converge << "\n";

        if (coxRidge.converge == false) {
            converge(k) = false;
            std::cout << "Warning: lambda " << cur_lambda << " failed to converge.\n";
        } else {
            converge(k) = true;
        }
        beta_old = coxRidge.beta;
        nulldev_old = coxRidge.get_null_deviance();
        beta_mx.col(k) = coxRidge.beta;
        eta_mx.col(k) = coxRidge.eta;
        deviance(k) = coxRidge.get_deviance();
        dev_ratio(k) = coxRidge.dev_ratio;
        object_val(k) = coxRidge.get_object();

        if (k > 4 && _user_define_lambda == false) {
            if (dev_ratio(k) > 0.99) { break; }
            if (k > 0 && (dev_ratio(k) - dev_ratio(k - 3)) <
                1e-3 * dev_ratio(k)) { break; }
        }
    }

    if (break_pt < _lambda_len) {
        beta_mx.conservativeResize(Xmat.cols(), break_pt);
        eta_mx.conservativeResize(survivalData.n, break_pt);
        deviance.conservativeResize(break_pt);
        dev_ratio.conservativeResize(break_pt);
        lambda_vec.conservativeResize(break_pt);
        object_val.conservativeResize(break_pt);
        converge.conservativeResize(break_pt);
    }
}
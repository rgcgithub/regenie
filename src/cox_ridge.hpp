#ifndef COXL2_H
#define COXL2_H
#include "Regenie.hpp"

class cox_ridge_path {
    public:
        // coefficients
        Eigen::MatrixXd beta_mx;
        Eigen::MatrixXd eta_mx;
        
        Eigen::VectorXd lambda_vec;

        // objective value
        Eigen::VectorXd object_val;
        Eigen::VectorXd deviance;
        Eigen::VectorXd dev_ratio;
        Eigen::Array<bool, Eigen::Dynamic, 1> converge;
        cox_ridge_path(const survival_data& survivalData, const Eigen::MatrixXd& Xmat, const Eigen::VectorXd& offset_val, const ArrayXb& mask, const int& nlambda = 100, const double& lambda_min_max_ratio = -1, const Eigen::VectorXd& lambda = Eigen::VectorXd(), const int& max_iter = 100, const int& max_inner_iter = 30, const double& tolerance = 1e-6, const bool& verbose_fit = false);
        void fit(const survival_data& survivalData, const Eigen::MatrixXd& Xmat, const Eigen::VectorXd& offset_val, const ArrayXb& mask);
        
        // fitting info
        int niter, mxitnr;
        double tol;
        bool verbose;

    private:
        double _getCoxLambdaMax(const Eigen::MatrixXd& Xmat, const Eigen::VectorXd& gradient);
        int _lambda_len;
        bool _user_define_lambda = false;
};

class cox_ridge {
    public:
        // coefficients
        Eigen::VectorXd beta;
        // prediction
        Eigen::VectorXd eta, eta_order;
        double lambda;
        bool converge;
        double dev_ratio;

        cox_ridge(const survival_data& survivalData, const Eigen::MatrixXd& Xmat, const Eigen::VectorXd& offset_val, const ArrayXb& mask, const double& lambda_val, const int& max_iter = 100, const int& max_inner_iter = 30, const double& tolerance = 1e-6, const bool& verbose_obj = false, const Eigen::VectorXd& beta_init = Eigen::VectorXd(), const double& null_deviance = -999);
        void fit(const survival_data& survivalData, const Eigen::MatrixXd& Xmat, const Eigen::VectorXd& offset_val, const ArrayXb& mask);
        void reset(const survival_data& survivalData, const Eigen::MatrixXd& Xmat, const Eigen::VectorXd& offset_val, const ArrayXb& mask, const double& lambda_val, const Eigen::VectorXd& beta_init = Eigen::VectorXd(), const double& null_deviance = -999);
        void coxGrad(const survival_data& survivalData);
        Eigen::VectorXd get_gradient();
        Eigen::VectorXd get_object_all();
        Eigen::VectorXd get_deviance_all();
        double get_object();
        double get_deviance();
        double get_null_deviance();

    private:
        // gradient
        Eigen::VectorXd _gradient, _diagHessian;
        int _niter, _mxitnr;
        double _tol;
        bool _verbose;
        // objective value
        Eigen::VectorXd _object;
        Eigen::VectorXd _deviance;
        
        double _coxDeviance(const survival_data& survivalData);
        double _coxLoglik(const survival_data& survivalData);
};

#endif
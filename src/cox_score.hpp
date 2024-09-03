#ifndef COXSCORE_H
#define COXSCORE_H

class cox_mle {
    public:
        int p;
        // coefficients
        Eigen::VectorXd beta;
        Eigen::VectorXd lambda0;
        Eigen::VectorXd mu;
        Eigen::VectorXd residual;
        Eigen::VectorXd Y;
        Eigen::MatrixXd XtWX;
        Eigen::MatrixXd sqrtWX;
        double loglik_val;
		Eigen::VectorXd loglike;
        // prediction
        Eigen::VectorXd eta, eta_order;
        int iter;
        bool converge;

        // prepare for test
        int n_events_dd;
        Eigen::VectorXd w_exp_eta;
        Eigen::MatrixXd UhalfX;
        Eigen::VectorXd Dhalf;
        Eigen::MatrixXd cov_inv;
        Eigen::MatrixXd X1_X1WX1inv;
        Eigen::MatrixXd WX1;
        double res_var;
        // Eigen::MatrixXd X_XtVXinv;
        // Eigen::MatrixXd X1_X1tX1invhalf;

        cox_mle();
        void setup(const survival_data& survivalData, const Eigen::MatrixXd& Xmat, const Eigen::VectorXd& offset_val, const ArrayXb& mask, const int& max_iter = 100, const int& max_inner_iter = 30, const double& tolerance = 1e-6, const bool& verbose_obj = false, const Eigen::VectorXd& beta_init = Eigen::VectorXd(), const Eigen::VectorXd& eta_init = Eigen::VectorXd());
        void compute_loglike(const survival_data& survivalData);
        void fit(const survival_data& survivalData, const Eigen::MatrixXd& Xmat, const Eigen::VectorXd& offset_val, const ArrayXb& mask);
        void cox_test_prep(const survival_data& survivalData, const Eigen::MatrixXd& Xmat, const Eigen::VectorXd& offset_val, const ArrayXb& mask);
        
    private:
        int _niter, _mxitnr;
        double _tol;
        bool _verbose;
};

#endif
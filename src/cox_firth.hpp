#ifndef COXFIRTH_H
#define COXFIRTH_H

class cox_firth {
    public:
		int p;
        // coefficients
        Eigen::VectorXd beta;
        // prediction
        Eigen::VectorXd eta, eta_order, residual;
        int iter;
        bool converge = false;

        // prepare for test
        Eigen::VectorXd exp_eta;
        Eigen::VectorXd mu;
		Eigen::VectorXd loglike;
		double loglik_val;
		Eigen::VectorXd first_der;
		Eigen::MatrixXd second_der;
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qrsd, qrsd_incl;
        double first_der_1;
		double second_der_1;

        cox_firth();
        void setup(const survival_data& survivalData, const Eigen::MatrixXd& Xmat, const Eigen::VectorXd& offset_val, const int& cols_incl, const int& max_iter = 100, const int& max_inner_iter = 30, const double& tolerance = 1e-6, const double& beta_tol = 1e-6, const double& max_step = 1, const bool& use_firth = true, const bool& verbose_obj = false, const Eigen::VectorXd& beta_init = Eigen::VectorXd());
        void cox_firth_likelihood(const survival_data& survivalData, const Eigen::MatrixXd& Xmat);
		void fit(const survival_data& survivalData, const Eigen::MatrixXd& Xmat, const Eigen::VectorXd& offset_val);
        void cox_firth_likelihood_1(const survival_data& survivalData, const Eigen::VectorXd& g);
        void fit_1(const survival_data& survivalData, const Eigen::VectorXd& g, const Eigen::VectorXd& offset_val);

    private:
        int _niter, _mxitnr, _cols_incl;
        double _tol, _betatol;
		double _maxstep;
        bool _usefirth, _verbose;
};

#endif
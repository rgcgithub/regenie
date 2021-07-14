/*

   This file is part of the regenie software package.

   Copyright (c) 2020-2021 Joelle Mbatchou, Andrey Ziyatdinov & Jonathan Marchini

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in all
   copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

*/

#include "Regenie.hpp"
#include "Files.hpp"
#include "Geno.hpp"
#include "Step1_Models.hpp"
#include "Pheno.hpp"
#include "HLM.hpp"

using namespace std;
using namespace Eigen;
using namespace boost;
using namespace LBFGSpp;
using boost::math::normal;
using boost::math::chi_squared;

HLM::HLM(){
}

HLM::~HLM(){
}

void HLM::prep_run(struct phenodt const* pheno_data, struct param const* params){

  // store Vlin = (1, E) 
  allocate_mat(Vlin, pheno_data->interaction_cov.rows(), params->ncov_interaction + 1);
  Vlin << MatrixXd::Ones(pheno_data->interaction_cov.rows(),1), pheno_data->interaction_cov; 

  if(params->hlm_vquad && params->int_add_esq_term){

    // set V = (1, E, E^2) [ apply QR to U = (E,E^2), center & scale, and set V = (1, U)]
    MatrixXd Vtmp (Vlin.rows(), params->ncov_interaction * 2);
    Vtmp << pheno_data->interaction_cov, pheno_data->interaction_cov.array().square().matrix();
    //cerr << "pre:\n" << Vtmp.topRows(5) << "\n\n";
    apply_QR(Vtmp, params, true);
    allocate_mat(V, Vtmp.rows(), Vtmp.cols() + 1);
    V << MatrixXd::Ones(Vtmp.rows(),1), Vtmp; 
    //cerr << "post:\n" << V.topRows(5) << "\n\n";

  } else {

    // set V = (1, E) and center & scale E
    allocate_mat(V, Vlin.rows(), Vlin.cols());
    V = Vlin;
    rescale_mat(V.rightCols(params->ncov_interaction), params); 

  }

  // set X = (covs, E^2, blup) - covs may include E
  MatrixXd Xtmp (pheno_data->new_cov.rows(), params->ncov + (params->int_add_esq_term ? params->ncov_interaction : 0 ));
  if(params->int_add_esq_term) {
    //cerr << "pre:\n" << Xtmp.topRows(5) << "\n\n";
    Xtmp << pheno_data->new_cov, pheno_data->interaction_cov.array().square().matrix();
    apply_QR(Xtmp, params, false);
  } else
    Xtmp = pheno_data->new_cov;
  //cerr << "post:\n" << Xtmp.topRows(5) << "\n\n";

  allocate_mat(X, Xtmp.rows(), Xtmp.cols() + (params->skip_blups ? 0 : 1) );
  X.leftCols(Xtmp.cols()) = Xtmp; 

  // for projection under null
  Px.resize(params->n_pheno);
  allocate_mat(yres, pheno_data->interaction_cov.rows(), params->n_pheno);

}

// For each phenotype, fit the null HLM 
//  Y = Xa + e, where e ~ N(0, exp(Vb) )
//  Have this be outside of the class so can use LBFGS solver
void HLM_fitNull(HLM& nullHLM, struct ests const& m_ests, struct phenodt const& pheno_data, struct in_files const& files, struct param const& params, mstream& sout){

  // if no blup predictions are given, this should only be ran once
  if(params.skip_blups && !nullHLM.first_fit) 
    return;

  sout << "   -fitting null HLMs for each trait..." << flush;
  auto t1 = std::chrono::high_resolution_clock::now();

  double fx;
  VectorXd beta(nullHLM.V.cols());
  allocate_mat(nullHLM.Dinv_sqrt, pheno_data.phenotypes_raw.rows(), pheno_data.phenotypes_raw.cols());

  LBFGSParam<double> bfgs_param;
  bfgs_param.max_iterations = nullHLM.max_iter;
  bfgs_param.max_linesearch = nullHLM.linesearch_try; // use more lenient number
  LBFGSSolver<double> solver(bfgs_param);

  
  for(int i = 0; i < params.n_pheno; i++){

    nullHLM.n = pheno_data.Neff(i);
    nullHLM.mask = pheno_data.masked_indivs.col(i);
    nullHLM.y = pheno_data.phenotypes_raw.col(i);
    if(!params.skip_blups) // add blup as a covariate
      nullHLM.X.rightCols(1) = m_ests.blups.col(i);
    beta.array() = 0;

    try {
      
      // get starting value for b
      nullHLM.get_alpha(beta);
      nullHLM.get_beta_approx(beta);
      // LBFGS
      solver.minimize(nullHLM, beta, fx);
      nullHLM.store_null_est(i);

    } catch(...){

      // redo with higher number of line search trials
      try {

        if(params.verbose) 
          sout << "Retrying HLM null model fitting for " << files.pheno_names[i] << endl;

        LBFGSParam<double> bfgs_param_retry;
        bfgs_param_retry.max_iterations = nullHLM.max_iter_retry;
        bfgs_param_retry.max_linesearch = nullHLM.linesearch_retry; 
        bfgs_param_retry.max_step = nullHLM.max_step_retry; 
        LBFGSSolver<double> solver_retry(bfgs_param_retry);

        // get starting value for b
        beta.array() = 0;
        nullHLM.get_alpha(beta);
        nullHLM.get_beta_approx(beta);
        // LBFGS
        solver_retry.minimize(nullHLM, beta, fx);
        nullHLM.store_null_est(i);

      } catch(...){
        throw "LFBGS could not fit HLM null model for trait " + files.pheno_names[i];
      }

    }
    //cerr << "\nFinal--\nalpha=\n"<<nullHLM.alpha << "\n\nbeta=\n" << beta <<"\n\nfx=" << fx << "\t" << std::boolalpha << isnan(fx);

  }

  //cerr << "\n\n" << nullHLM.yres.topRows(5)<<"\n\n";

  sout << "done";
  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  sout << " (" << duration.count() << "ms) "<< endl;

  nullHLM.first_fit = false;
}


void HLM::get_alpha(VectorXd const& beta){

  Vb = (V * beta).array();
  Dinv = (-Vb).exp() * mask.cast<double>();
  if( (Dinv == 0).all() ) // will cause Xd = 0
    throw std::underflow_error("D=0 occurred");
  MatrixXd Xd = (X.array().colwise() * Dinv).matrix().transpose();
  alpha = (Xd * X).colPivHouseholderQr().solve( Xd * y );

}

void HLM::get_beta_approx(VectorXd& beta){

  ArrayXd esq = ((y - X * alpha).array() * mask.cast<double>()).square();
  //cerr << "\nE=\n" << esq.head(10) << "\n\n";
  beta = (V.transpose() * esq.matrix().asDiagonal() * V).colPivHouseholderQr().solve( V.transpose() * ((esq - 1) * mask.cast<double>()).matrix() );
  //cerr << "alpha:\n" << alpha << "\n\nbeta:\n" << beta <<"\n\n";

}

// get projection matrix
void HLM::store_null_est(int const& ph){

  Dinv_sqrt.col(ph) = Dinv.sqrt().matrix();
  MatrixXd Xd = (X.array().colwise() * Dinv_sqrt.col(ph).array()).matrix();
  SelfAdjointEigenSolver<MatrixXd> es(Xd.transpose() * Xd);
  VectorXd eigD = es.eigenvalues();

  Px[ph] = ((Xd * es.eigenvectors()).array().rowwise() / eigD.transpose().array().sqrt()).matrix();
  //cerr << "\nP=\n" << Px[ph].topRows(5) << "\n\n";

  residualize(ph, y, yres.col(ph));
  
}

void HLM::residualize(int const& ph, Ref<MatrixXd> mat_orig, Ref<MatrixXd> mat_res){

  MatrixXd m = (mat_orig.array().colwise() * Dinv_sqrt.col(ph).array()).matrix();
  //cerr << "Y" << ph+1 << "\norig:\n" << print_mat_dims(mat_orig) << 
   // "\nm:\n" << print_mat_dims(m) << "\nPx:\n" << print_mat_dims(Px[ph]) << endl;
  mat_res = m -  Px[ph] * (Px[ph].transpose() * m);

}

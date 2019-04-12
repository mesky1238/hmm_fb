#include <limits>
#include <RcppEigen.h>
#include "hmmfb.h"

typedef std::numeric_limits<double> double_lim;

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R

using namespace Rcpp;

typedef Eigen::Array<double, Eigen::Dynamic, 1> ColumnArray;
typedef Eigen::Array<double, 1,Eigen::Dynamic > RowArray;

//[[Rcpp::export]]
Eigen::MatrixXd sweep_add (const Eigen::MatrixXd &A,
			   const Eigen::VectorXd &alpha){
  return(A.rowwise()+alpha.transpose());
}

//[[Rcpp::export]]
Eigen::ArrayXd sweep_add_logexp (const Eigen::MatrixXd &A,
				 const Eigen::VectorXd &alpha){
  Eigen::MatrixXd tA = A.rowwise()+alpha.transpose();
  ColumnArray tmax= A.rowwise().maxCoeff();
  return((tA.array().colwise()-tmax).exp().rowwise().sum().log()+tmax);
}

//[[Rcpp::export]]
double logsumexp(const Eigen::ArrayXd &x){
  auto m=x.minCoeff();
  return(m+std::log((x-m).exp().sum()));
}


Eigen::ArrayXd fbiter(const Eigen::Ref<const Eigen::MatrixXd> A,
		      const Eigen::Ref<const Eigen::ArrayXd> alpha, 
		      const Eigen::Ref<const Eigen::ArrayXd> B){
  auto tA=A.array().rowwise()+alpha.transpose();
  return(sweep_add_logexp(A,alpha)+B);
}

//[[Rcpp::export(name="fbiter")]]
Eigen::ArrayXd efbiter(Eigen::MatrixXd A, Eigen::ArrayXd alpha, 
Eigen::ArrayXd B){
  return(fbiter(A,alpha,B));
}

//[[Rcpp::export]]
Eigen::ArrayXi cpp_cumsum(Eigen::ArrayXi x){
  Eigen::ArrayXi csum(x.size());
  for(int i=0; i<x.size();i++){
    csum(i)=x.head(i+1).sum();
  }
  return(csum);
}

//[[Rcpp::export]]
Eigen::ArrayXi gen_bt(const Eigen::MatrixXd &linit,const Eigen::MatrixXd &lA,
		      const Eigen::MatrixXd &B, Eigen::ArrayXi ntimes){
  size_t lt =ntimes.size();
  size_t ns=linit.cols();
  size_t nt=B.rows();
  
  Eigen::ArrayXi et=cpp_cumsum(ntimes);  
  Eigen::ArrayXi bt(et.size());
  bt(0)=0;
  if(et.size()>1){
    bt.tail(et.size()-1)=et.head(et.size()-1);
  }
  return(bt);
}

//[[Rcpp::export]]
Rcpp::List forward_rcpp (const Eigen::MatrixXd &linit,
			 const Eigen::MatrixXd &lA,
			 const Eigen::MatrixXd &B, 
			 Eigen::ArrayXi ntimes,
			 bool return_all=false) {
  using namespace Rcpp;
  using namespace Eigen;
  
  // # Forward-Backward algorithm (used in Baum-Welch)
  // # Returns alpha, beta, and full data likelihood
  //
  // # NOTE THE CHANGE IN FROM ROW TO COLUMN SUCH THAT TRANSPOSING A IS NOT NECCESSARY ANYMORE
  // # IN COMPUTING ALPHA AND BETA BUT IS NOW NECCESSARY IN COMPUTING XI
  // # A = T*K*K array with transition probabilities, from row to column!!!!!!!
  // # B = T*D*K matrix with elements ab_{ij} = P(y_i|s_j)
  // # init = N*K vector with initial probabilities
  //
  // # T = total number of time points
  // # K = number of states
  // # D = dimension of observations (D>1 is multivariate)
  // # N = number of participants
  
  size_t lt =ntimes.size();
  size_t ns=linit.cols();
  size_t nt=B.rows();
  
  Eigen::ArrayXi et=cpp_cumsum(ntimes);  
  Eigen::ArrayXi bt(et.size());
  bt(0)=0;
  if(et.size()>1){
    bt.tail(et.size()-1)=et.head(et.size()-1);
  }
  
  Eigen::MatrixXd alpha(ns,nt);
  
  alpha.setZero();
  
  Eigen::ArrayXd sca(nt);
  
  for(int i=0; i<lt; i++){
    int col_num =bt(i);
    alpha.col(col_num) = linit.row(i)+B.row(col_num);
    sca(col_num) = -logsumexp(alpha.col(col_num));
    alpha.col(col_num) = alpha.col(col_num).array()+sca(col_num);
    
    if(ntimes(i)>1){
      int ecase=et(i)-1;
      if(col_num<ecase){
        for(int j=col_num;j<ecase;j++){
          alpha.col(j+1) = fbiter(lA,alpha.col(j),B.row(j+1));
          sca(j+1)=-logsumexp(alpha.col(j+1));
          alpha.col(j+1).array()+= sca(j+1);
        }
      }else{
        for(int j=ecase-1;j>=col_num;j--){
          alpha.col(j+1) = fbiter(lA,alpha.col(j),B.row(j+1));
          sca(j+1)=-logsumexp(alpha.col(j+1));
          alpha.col(j+1).array()+= sca(j+1);
        }
      }
    }
  }
  double like = -(sca.sum());
  
  if(return_all){
    return Rcpp::List::create(_["alpha"]=alpha.transpose(),
			      _["sca"]=sca,_["logLike"] = like);
  }else{
    return Rcpp::List::create(_["logLike"]=like);    
  }
}

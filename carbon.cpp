// Separable covariance on lattice with AR1 structure in each direction.
#include <TMB.hpp>


template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_ARRAY(p_c);
  DATA_ARRAY(p_b);
  DATA_VECTOR(c_i);
  DATA_VECTOR(temp);
  
  PARAMETER_VECTOR(beta_c);
  PARAMETER_MATRIX(A_par);
  PARAMETER_VECTOR(alpha);
  PARAMETER_VECTOR(gamma);
  PARAMETER(ln_sig);

  Type ff = 0.;  //joint likelihood
  
  int n_c = p_c.dim[0];
  int n_i = p_c.dim[1];
  int n_b = p_b.dim[0];

  vector<Type> beta(n_c);
  Type beta_sum = 0.;
  beta(0) = 0.;
  for(int cc=1;cc<n_c;cc++){
    beta(cc) = beta_c(cc-1);
    beta_sum += exp(beta(cc));
  }

  matrix<Type> A(n_c,n_b);
  Type A_sum = 0.;
  A.col(0).fill(0.);
  for(int bb=1;bb<n_b;bb++){
    A.col(bb) = A_par.col(bb-1);
  }

  matrix<Type> phat_c(n_c,n_i);
  for(int i=0;i<n_i;i++){
    for(int cc=0;cc<n_c;cc++){
      phat_c(cc,i) = exp(beta(cc))/(1+beta_sum);
      ff -= p_c(cc,i)*log(phat_c(cc,i));
    }
  }

  matrix<Type> phat_b(n_b,n_i);
  vector<Type> tmp_b(n_b);
  Type b_sum;
  for(int i=0;i<n_i;i++){
    b_sum = 0.;
    tmp_b = 1.;
    for(int bb=0;bb<n_b;bb++){
      for(int cc=0;cc<n_c;cc++){
        tmp_b(bb) *= exp(A(cc,bb)*phat_c(cc,i));
      }
      b_sum += tmp_b(bb);
    }
    for(int bb=0;bb<n_b;bb++){
      phat_b(bb,i) = tmp_b(bb)/b_sum;
      ff -= p_b(bb,i)*log(phat_b(bb,i));
    }
  }
  
  vector<Type> chat_i(n_i);
  chat_i.fill(0.);
  for(int i=0;i<n_i;i++){
    b_sum = 0.;
    tmp_b = 1.;
    for(int bb=0;bb<n_b;bb++){
      for(int cc=0;cc<n_c;cc++){
        tmp_b(bb) *= exp(A(cc,bb)*phat_c(cc,i));
      }
      b_sum += tmp_b(bb);
    }
    for(int bb=0;bb<n_b;bb++){
      phat_b(bb,i) = tmp_b(bb)/b_sum;
      chat_i(i) += alpha(bb) + gamma(bb)*temp(i)*phat_b(bb,i);
      ff -= p_b(bb,i)*log(phat_b(bb,i));
    }
    ff -= dnorm(log(c_i(i)),chat_i(i),exp(ln_sig),true);
  }
  
  REPORT(beta);
  REPORT(A);
  REPORT(alpha);
  REPORT(gamma);
  REPORT(phat_b);
  REPORT(phat_c);
  REPORT(chat_i);

  return ff;
}

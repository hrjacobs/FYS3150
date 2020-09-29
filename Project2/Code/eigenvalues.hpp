#ifndef EIGENVALUES_HPP
#define EIGENVALUES_HPP

#include <armadillo>

using namespace arma;

class eigenvalues {
  private:
    int n;
    int N;
    int p, q;

    double m_max;
    double s,c;
    double t,tau;
    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;

    mat R;
    mat A;

    bool running;
    bool col_swap;

    void Initialize(mat A_, int n_){
      // Initializing step size, numbers etc.
      n = n_;
      A = A_;
      R = zeros(n,n);
      for (int i = 0; i < n; i++){
        R(i,i) = 1;
      }
    }
    // overload function
    void Initialize(){
      n = 4;
      A = zeros(n);
      for (int i = 0; i < n; i++){
        if (i != 0){
          A(i,i-1) = -1;
        }
        A(i,i) = 2;
        if (i != n-1){
          A(i,i+1) = -1;
        }
      }
      std::cout << "Running on overload function\n" << std::endl;
      Initialize(A,n);
    }

  public:
    void offdiag();
    void Jacobi_rotate();
    void solve(double tolerance, int maxiter);

    vec get_eigenvectors(int n_);
    float get_eigenvalues(int n_);

    void order_eigenvalues();

    // setting up the overload
    eigenvalues(mat A, int n_){
      Initialize(A,n_);
    }
    eigenvalues(){
      Initialize();
    }
};
#endif

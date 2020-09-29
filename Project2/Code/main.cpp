#include <iostream>
#include <armadillo>
#include "eigenvalues.hpp"

using namespace arma;
using namespace std;

int main(int argc, char const *argv[]){
  int n = 100;
  int N = n+1;
  float h = 1./N;
  float a = -1*1./pow(h,2);
  float d = 2*1./pow(h,2);

  mat A = zeros(n,n);
  for (int i = 0; i < n; i++){
    if (i != 0){
      A(i,i-1) = a;
    }
    A(i,i) = d;
    if (i != n-1){
      A(i,i+1) = a;
    }
  }
  double tolerance = 1.0E-5;
  int maxiter = 100000;

  eigenvalues test(A,n);    // class defined in eigenvalues.hpp
  test.solve(tolerance,maxiter);
  test.order_eigenvalues();
// Printing to file as rows: eigvals, eigvec^T
  ofstream outfile ("eigenpairs.txt");
  outfile << "eigenvalue first element in row, then corresponding eigenvector." << endl;
  for (int i = 0; i < n; i++){
      float eigval = test.get_eigenvalues(i);
      vec eigvec = test.get_eigenvectors(i);
      outfile << eigval << eigvec.t() << endl;
  }
  outfile.close();

}

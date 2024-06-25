#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
double numer_C(NumericVector Y, double gam) {
  double rst = 0;
  for (int i = 0; i < Y.size(); ++i) {
    rst += Y[i] * std::pow(gam, i);
  }
  return rst;
}

// [[Rcpp::export]]
double Cy(NumericVector Y, double gam) {
  return numer_C(Y, gam) / (1 - std::pow(gam, 2 * Y.size())) * (1 - std::pow(gam, 2));
}

// [[Rcpp::export]]
double Dy(NumericVector Y, double gam) {
  double cy = Cy(Y, gam);
  double nc = numer_C(Y, gam);
  double len = Y.size();
  return sum(pow(Y, 2)) / 2 - cy * nc + std::pow(cy, 2) / 2 * (1 - std::pow(gam, 2 * len)) / (1 - std::pow(gam, 2));
}


// [[Rcpp::export]]
List estspike_vanilla(NumericMatrix dat, double gam, double lam, int trial, double power, NumericVector st_gauss) {
  NumericVector Y = dat(trial - 1, _);
  //Rcout << "The value of Y1: " << Y[0] << "\n";
  int n = Y.size();
  //Rcout << "The value of n: " << n << "\n";
  NumericVector Fset(n + 1); // create vector object
  Fset[0] = -lam;
  std::vector< std::set<int> > cp(n + 1);
  cp[0].insert(0);
  
  double pen = lam;
  NumericVector w_t = exp(-pow(st_gauss, power));
  NumericVector lam_t = w_t / sum(w_t) * lam * n;
  double pen1 = (power != 0) ? lam_t[0] : lam;
  
  //Rcout << "Okay!\n";
  //Rcout << Y[seq(0,1)].size() << '\n';
  //Rcout << Y[0] << '\n';
  //Rcout << n << '\n';
  for (int i = 2; i <= (n+1); ++i) { // i = 2, 3, ..., (n+1)
    double Fmin = Fset[0] + Dy(Y[seq(0, i - 2)], gam) + pen1; // index of Y starts from 0!
    //if (i < 10){
    //Rcout << i <<"Okay!\n";
    //}
    int sprime = 1;
    
    if (i > 2) {
      for (int j = 2; j <= (i-1); ++j) {
        if (power != 0) {
          pen = lam_t[j - 1];
        }
        double Fset_temp = Fset[j - 1] + Dy(Y[seq(j - 1, i - 2)], gam) + pen;
        if (Fset_temp <= Fmin) {
          Fmin = Fset_temp;
          sprime = j;
        }
      }
    }
    //Rcout << "Okay!\n";
    Fset[i-1] = Fmin;
    cp[i-1] = cp[sprime - 1];
    cp[i-1].insert(sprime - 1);
  }
  
  std::vector<int> cpset(cp[n].begin(), cp[n].end());
  for (int &c : cpset) {
    c++;
  }
  
  NumericVector ct(n);
  if (cpset.size() <= 1) {
    ct[0] = Cy(Y, gam);
    for (int i = 1; i < n; ++i) {
      ct[i] = ct[0] * std::pow(gam, i);
    }
    cpset.clear();
  } else {
    cpset.erase(cpset.begin());
    ct[0] = Cy(Y[seq(0, cpset[0] - 1)], gam);
    for (int i = 1; i < cpset[0]; ++i) {
      ct[i] = ct[0] * std::pow(gam, i);
    }
    for (size_t i = 0; i < cpset.size(); ++i) {
      if (i == cpset.size() - 1) {
        ct[cpset[i]] = Cy(Y[seq(cpset[i], n - 1)], gam);
        for (int j = cpset[i]; j < n; ++j) {
          ct[j] = ct[cpset[i]] * std::pow(gam, j - cpset[i]);
        }
      } else {
        ct[cpset[i]] = Cy(Y[seq(cpset[i], cpset[i + 1] - 1)], gam);
        for (int j = cpset[i]; j < cpset[i + 1]; ++j) {
          ct[j] = ct[cpset[i]] * std::pow(gam, j - cpset[i]);
        }
      }
    }
  }
  
  std::vector<int> rm_idx;
  for (size_t i = 0; i < cpset.size(); ++i) {
    if (ct[cpset[i]] - ct[cpset[i] - 1] < 0) {
      rm_idx.push_back(i);
    }
  }
  
  //Rcout << rm_idx.size() << '\n';
  for (size_t i = 0; i < rm_idx.size(); ++i) {
    cpset.erase(cpset.begin() + rm_idx[i] - i);
  }
  
  NumericVector st(n);
  for (int i = 1; i < n; ++i) {
    st[i] = ct[i] - gam * ct[i - 1];
  }
  
  return List::create(
    Named("cp") = cpset,
    Named("ct") = ct,
    Named("st") = st,
    Named("lam_t") = lam_t,
    Named("Fset") = Fset,
    Named("cpall") = cp
  );
}

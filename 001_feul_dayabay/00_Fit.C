#include<iostream>
using namespace std;

#include<fstream>
#include<cstdlib>
#include "math.h"

// Minuit2
#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"

// ROOT
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TMatrixD.h"

// ----------------------------------------------------- Global variables

#define Countor

const int sigma_AA_N = 8;
const int iso_N = 4;

double sigmaf_AA_fit[iso_N] = {0};
TMatrixD sigma_AA(sigma_AA_N, 1);
TMatrixD f_AA(sigma_AA_N, iso_N);
TMatrixD sigmaf_AA(iso_N, 1); // sigma_AA - f_AA * sigmaf_AA;
TMatrixD cov_AA(sigma_AA_N, sigma_AA_N);

// -----------------------------------------------------
  double vratio( double a, double ea, double b, double eb )
{
  return a/b;
}

double eratio( double a, double ea, double b, double eb )
{
  double r=a/b;
  double dr=sqrt( r*r*( ea*ea/a/a + eb*eb/b/b )  );
  return dr;
}

// 
double fcn(const double *par);

// ----------------------------------------------------- ----------------------------------------------------- MAIN

int main()
{
  TString roostr = "";

  // -----------------------------------------------------

  ifstream file_AA_flux("./inputs/flux_data.dat");
  for(int idx=0; idx<8; idx++) {
    double lowF, hghF, f239, f235, f238, f241, sigma;    
    file_AA_flux>>lowF>>hghF>> f239>> f235>> f238>> f241>> sigma;

    sigma_AA[idx][0] = sigma;
    f_AA[idx][0] = f239;
    f_AA[idx][1] = f235;
    f_AA[idx][2] = f238;
    f_AA[idx][3] = f241;
  }

  ifstream file_AA_covstat("./inputs/cov_stat.dat");
  ifstream file_AA_covsyst("./inputs/cov_syst.dat");
  for( int i=0; i<8; i++ )
    for( int j=0; j<8; j++) {

      double f239_array[8] = {0.344, 0.332, 0.321, 0.311, 0.299, 0.287, 0.274, 0.252};
      double sigma_bar = 5.90;
      double ds2df = -1.86;
      double Fbar = 0.299;
      double factor_a = sigma_bar + ds2df * (  f239_array[i] - Fbar );
      double factor_b = sigma_bar + ds2df * (  f239_array[j] - Fbar );
      
      double temp1 = 100;
      double temp2 = 100;

      file_AA_covstat >> temp1;
      file_AA_covsyst >> temp2;
      
      //factor_a = sigma_AA[i][0];
      //factor_b = sigma_AA[j][0];

      cout<<factor_a<<"\t"<<sigma_AA[i][0]<<endl;

      cov_AA[i][j] += temp1;
      cov_AA[i][j] += temp2 * factor_a * factor_b;

      //cout<<i<<"\t"<<j<<"\t"<<temp1<<"\t"<<temp2<<endl;
    }
  
  // -----------------------------------------------------

  cout<<"Fit initialization"<<endl;
  
  ROOT::Minuit2::Minuit2Minimizer min( ROOT::Minuit2::kMigrad );
  min.SetPrintLevel(2);
  min.SetStrategy(1); //0- cursory, 1- default, 2- thorough yet no more successful
  min.SetMaxFunctionCalls(500000);
  min.SetMaxIterations(500000);
  min.SetTolerance(2e-6); // tolerance*2e-3 = edm precision
  //min.SetTolerance(2e-8); // tolerance*2e-3 = edm precision
  min.SetPrecision(1e-18); //precision in the target function

  ROOT::Math::Functor Chi2Functor( &fcn, 4 );
  min.SetFunction(Chi2Functor);

  int par_index = 0;
  min.SetVariable( par_index, "sigmaDYB_239", 4, 1e-4 ); par_index++;
  min.SetVariable( par_index, "sigmaDYB_235", 6, 1e-4 ); par_index++;
  min.SetVariable( par_index, "sigmaDYB_238", 11, 1e-4 ); par_index++;
  min.SetVariable( par_index, "sigmaDYB_241", 8, 1e-4 ); par_index++;
 
  min.Minimize();
  const double *par_val = min.X();
  const double *par_err = min.Errors();
  
  for( int idx=0; idx<par_index; idx++ ) {
    ofstream out_AA_best("out_AA_best.dat", ios::out | ios::app);
    out_AA_best<<par_val[idx]<<"\t"<<par_err[idx]<<endl;
    out_AA_best.close();
  }
  
  /// https://www.fourmilab.ch/rpkp/experiments/analysis/chiCalc.html
  //min.SetErrorDef( 2.2977 );// 1 - 0.317 with ndf = 2
  //min.SetErrorDef( 6.2021 );// 1 - 0.045 with ndf = 2
  //min.SetErrorDef(11.6182 );// 1 - 0.003 with ndf = 2
  
#ifdef Countor
  const int Npoints = 100;
  unsigned int Npoints_usr = 100;

  double fit_AA_s239_1s[Npoints] = {0};
  double fit_AA_s235_1s[Npoints] = {0};

  double fit_AA_s239_2s[Npoints] = {0};
  double fit_AA_s235_2s[Npoints] = {0};

  double fit_AA_s239_3s[Npoints] = {0};
  double fit_AA_s235_3s[Npoints] = {0};

  min.SetErrorDef( 2.2977 ); min.Contour(0,1, Npoints_usr, fit_AA_s239_1s, fit_AA_s235_1s);
  min.SetErrorDef( 6.2021 ); min.Contour(0,1, Npoints_usr, fit_AA_s239_2s, fit_AA_s235_2s);
  min.SetErrorDef(11.6182 ); min.Contour(0,1, Npoints_usr, fit_AA_s239_3s, fit_AA_s235_3s);

  ofstream out_AA_1s("out_AA_1s.dat", ios::out | ios::app);
  ofstream out_AA_2s("out_AA_2s.dat", ios::out | ios::app);
  ofstream out_AA_3s("out_AA_3s.dat", ios::out | ios::app);
  for(int idx=0; idx<Npoints; idx++) {
    out_AA_1s<<fit_AA_s239_1s[idx]<<"\t"<<fit_AA_s235_1s[idx]<<endl;
    out_AA_2s<<fit_AA_s239_2s[idx]<<"\t"<<fit_AA_s235_2s[idx]<<endl;
    out_AA_3s<<fit_AA_s239_3s[idx]<<"\t"<<fit_AA_s235_3s[idx]<<endl;
  }

  out_AA_1s.close();
  out_AA_2s.close();
  out_AA_3s.close();
#endif  

  cout<<endl<<endl;
  
  for( int idx=0; idx<par_index; idx++ ) {
    cout<<TString::Format(" ---> idx %2d: %7.4f +/- %7.4f, ", 
			  idx+1, par_val[idx], par_err[idx])<<min.VariableName(idx)<<endl;
  }
  
  cout<<endl<<endl;

  return 0;
}
// ----------------------------------------------------- ----------------------------------------------------- 

double fcn(const double *par)
{
  double chi2_AA_tot = 0;
  double chi2_AA_main = 0;
  double chi2_AA_pull = 0;
  /*
  double ref_238 = 10.1;
  double ref_241 = 6.05;
  double err_238 = 0.1*ref_238;
  double err_241 = 0.1*ref_241;
  */
  
  double ref_238 = 10.1;
  double ref_241 = 6.05;
  double err_238 = 0.1*ref_238;
  double err_241 = 0.1*ref_241;
  

  int index_par = 0;
  sigmaf_AA_fit[index_par] = par[index_par]; index_par++;
  sigmaf_AA_fit[index_par] = par[index_par]; index_par++;
  sigmaf_AA_fit[index_par] = par[index_par]; index_par++;
  sigmaf_AA_fit[index_par] = par[index_par]; index_par++;

  //err_238 = 0.1 * sigmaf_AA_fit[2] ;
  //err_241 = 0.1* sigmaf_AA_fit[3] ;

  for( int idx=0; idx<iso_N; idx++ ) {
    sigmaf_AA[idx][0] = sigmaf_AA_fit[idx];
  }

  TMatrixD matrix_AA_diff = sigma_AA - f_AA * sigmaf_AA;
  TMatrixD matrix_AAT_diff(1, sigma_AA_N);
  matrix_AAT_diff.Transpose(matrix_AA_diff);
  TMatrixD cov_AA_invert = cov_AA;
  cov_AA_invert.Invert();

  TMatrixD matrin_AA_main = matrix_AAT_diff * cov_AA_invert * matrix_AA_diff;
  chi2_AA_main = matrin_AA_main[0][0];

  chi2_AA_pull += pow(sigmaf_AA_fit[2]-ref_238, 2)/err_238/err_238;
  chi2_AA_pull += pow(sigmaf_AA_fit[3]-ref_241, 2)/err_241/err_241;

  chi2_AA_tot = chi2_AA_main + chi2_AA_pull;

  return chi2_AA_tot;
}

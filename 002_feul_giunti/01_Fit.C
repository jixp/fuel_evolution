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

const int sigma_BB_N = 26;
const int iso_N = 4;

double sigmaf_BB_fit[iso_N] = {0};
TMatrixD f_BB(sigma_BB_N, iso_N);
TMatrixD sigmaf_BB(iso_N, 1); // sigma_BB - f_BB * sigmaf_BB;
TMatrixD cov_BB(sigma_BB_N, sigma_BB_N);
TMatrixD Ratio_exp_BB(sigma_BB_N, 1);

double a_sigmaf_BB_SH[iso_N] = {4.40, 6.69, 10.10, 6.03};// 239, 235, 238, 241: https://arxiv.org/pdf/1608.04096.pdf
TMatrixD sigmaf_BB_SH(iso_N, 1);

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
  
  sigmaf_BB_SH[0][0] = a_sigmaf_BB_SH[0];
  sigmaf_BB_SH[1][0] = a_sigmaf_BB_SH[1];
  sigmaf_BB_SH[2][0] = a_sigmaf_BB_SH[2];
  sigmaf_BB_SH[3][0] = a_sigmaf_BB_SH[3];

  ifstream file_BB_flux("./inputs/edit_table.dat");
  for( int i=1; i<=sigma_BB_N; i++ ) {
    int line;
    double f239, f235, f238, f241, Ratio, exp, corrA, corrB;    
    file_BB_flux >> line >> f235 >> f238>> f239>> f241>> Ratio>> exp>> corrA>> corrB;
    
    Ratio_exp_BB[i-1][0] = Ratio;

    f_BB[i-1][0] = f239;
    f_BB[i-1][1] = f235;
    f_BB[i-1][2] = f238;
    f_BB[i-1][3] = f241;

    for(int j=1; j<=sigma_BB_N; j++) {
      if( i == j ) {
	cov_BB[i-1][j-1] = exp;
      }
      else {
	if( (i==1||i==2) && (j==1||j==2) ) cov_BB[i-1][j-1] = corrA;

	if( (i==3||i==4) && (j==3||j==4) ) cov_BB[i-1][j-1] = corrA;
	if( (i==5||i==6||i==7) && (j==5||j==6||j==7) ) cov_BB[i-1][j-1] = corrA;
	if( (i==3||i==4) && (j==5||j==6||j==7) ) cov_BB[i-1][j-1] = corrB;
	if( (i==5||i==6||i==7) && (j==3||j==4) ) cov_BB[i-1][j-1] = corrB;

	if( (i==8||i==9||i==10) && (j==8||j==9||j==10) ) cov_BB[i-1][j-1] = corrA;

	if( (i==11||i==12||i==13) && (j==11||j==12||j==13) ) cov_BB[i-1][j-1] = corrA;
	if( (i==11||i==12||i==13) && (j==14) ) cov_BB[i-1][j-1] = corrB;
	if( (j==11||j==12||j==13) && (i==14) ) cov_BB[i-1][j-1] = corrB;

	if( (i==15||i==16) && (j==15||j==16) ) cov_BB[i-1][j-1] = corrA;
      }

    }

  }
  //cov_BB.Draw("colz text");
  
  // -----------------------------------------------------

  cout<<"Fit initialization"<<endl;
  
  ROOT::Minuit2::Minuit2Minimizer min( ROOT::Minuit2::kMigrad );
  min.SetPrintLevel(2);
  min.SetStrategy(1); //0- cursory, 1- default, 2- thorough yet no more successful
  min.SetMaxFunctionCalls(500000);
  min.SetMaxIterations(500000);
  min.SetTolerance(2e-8); // tolerance*2e-3 = edm precision
  //min.SetTolerance(2e-8); // tolerance*2e-3 = edm precision
  min.SetPrecision(1e-18); //precision in the target function

  ROOT::Math::Functor Chi2Functor( &fcn, 4 );
  min.SetFunction(Chi2Functor);

  int par_index = 0;
  min.SetVariable( par_index, "sigmaRCT_239", 4, 1e-4 ); par_index++;
  min.SetVariable( par_index, "sigmaRCT_235", 6, 1e-4 ); par_index++;
  min.SetVariable( par_index, "sigmaRCT_238", 11, 1e-4 ); par_index++;
  min.SetVariable( par_index, "sigmaRCT_241", 8, 1e-4 ); par_index++;
 
  min.Minimize();

  const double *par_val = min.X();
  const double *par_err = min.Errors();
  
  for( int idx=0; idx<par_index; idx++ ) {
    ofstream out_BB_best("out_BB_best.dat", ios::out | ios::app);
    out_BB_best<<par_val[idx]<<"\t"<<par_err[idx]<<endl;
    out_BB_best.close();
  }
  
  // -----------------------------------------------------

#ifdef Countor
  const int Npoints = 100;
  unsigned int Npoints_usr = 100;

  double fit_BB_s239_1s[Npoints] = {0};
  double fit_BB_s235_1s[Npoints] = {0};

  double fit_BB_s239_2s[Npoints] = {0};
  double fit_BB_s235_2s[Npoints] = {0};

  double fit_BB_s239_3s[Npoints] = {0};
  double fit_BB_s235_3s[Npoints] = {0};

  min.SetErrorDef( 2.2977 ); min.Contour(0,1, Npoints_usr, fit_BB_s239_1s, fit_BB_s235_1s);
  min.SetErrorDef( 6.2021 ); min.Contour(0,1, Npoints_usr, fit_BB_s239_2s, fit_BB_s235_2s);
  min.SetErrorDef(11.6182 ); min.Contour(0,1, Npoints_usr, fit_BB_s239_3s, fit_BB_s235_3s);

  ofstream out_BB_1s("out_BB_1s.dat", ios::out | ios::app);
  ofstream out_BB_2s("out_BB_2s.dat", ios::out | ios::app);
  ofstream out_BB_3s("out_BB_3s.dat", ios::out | ios::app);
  for(int idx=0; idx<Npoints; idx++) {
    out_BB_1s<<fit_BB_s239_1s[idx]<<"\t"<<fit_BB_s235_1s[idx]<<endl;
    out_BB_2s<<fit_BB_s239_2s[idx]<<"\t"<<fit_BB_s235_2s[idx]<<endl;
    out_BB_3s<<fit_BB_s239_3s[idx]<<"\t"<<fit_BB_s235_3s[idx]<<endl;
  }

  out_BB_1s.close();
  out_BB_2s.close();
  out_BB_3s.close();
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
  double chi2_BB_tot = 0;
  double chi2_BB_main = 0;
  double chi2_BB_pull = 0;
  
  double ref_238 = 10.10;
  double ref_241 = 6.03;
  double err_238 = 0.15*ref_238;
  double err_241 = 0.10*ref_241;

  int index_par = 0;
  sigmaf_BB_fit[index_par] = par[index_par]; index_par++;
  sigmaf_BB_fit[index_par] = par[index_par]; index_par++;
  sigmaf_BB_fit[index_par] = par[index_par]; index_par++;
  sigmaf_BB_fit[index_par] = par[index_par]; index_par++;

  for( int idx=0; idx<iso_N; idx++ ) {
    sigmaf_BB[idx][0] = sigmaf_BB_fit[idx];
  }
 
  TMatrixD matrix_sigma_BB = f_BB * sigmaf_BB;
  TMatrixD matrix_sigma_BB_SH =  f_BB * sigmaf_BB_SH;
  TMatrixD matrix_sigma_BB_exp(sigma_BB_N, 1);

  for( int idx=0; idx<sigma_BB_N; idx++ ) {
    matrix_sigma_BB_exp[idx][0] = Ratio_exp_BB[idx][0] * matrix_sigma_BB_SH[idx][0];
  }
  
  TMatrixD matrix_BB_diff = matrix_sigma_BB - matrix_sigma_BB_exp;
  TMatrixD matrix_BBT_diff(1, sigma_BB_N);
  matrix_BBT_diff.Transpose(matrix_BB_diff);

  TMatrixD edit_cov_BB = cov_BB;
  for( int i=0; i<sigma_BB_N; i++ ) {
    for( int j=0; j<sigma_BB_N; j++ ) {
      edit_cov_BB[i][j] = edit_cov_BB[i][j] * edit_cov_BB[i][j] * 0.01 * 0.01;
      edit_cov_BB[i][j] = edit_cov_BB[i][j] * matrix_sigma_BB_exp[i][0] * matrix_sigma_BB_exp[j][0];
    }
  }

  TMatrixD edit_cov_BB_invert = edit_cov_BB;
  edit_cov_BB_invert.Invert();

  TMatrixD matrin_BB_main = matrix_BBT_diff * edit_cov_BB_invert * matrix_BB_diff;
  chi2_BB_main = matrin_BB_main[0][0];

  chi2_BB_pull += pow(sigmaf_BB_fit[2]-ref_238, 2)/err_238/err_238;
  chi2_BB_pull += pow(sigmaf_BB_fit[3]-ref_241, 2)/err_241/err_241;

  chi2_BB_tot = chi2_BB_main + chi2_BB_pull; 


  return chi2_BB_tot;
}

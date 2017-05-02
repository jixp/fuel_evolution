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

#define Fit_DYB
#define Fit_RCT
//#define Countor

const int iso_N = 4;
int cycle = 0;

/*
// RCT only
double ref_238 = 10.10;
double ref_241 = 6.03;
double err_238 = 0.15*ref_238;
double err_241 = 0.10*ref_241;
  
// DYB only
double ref_238 = 10.1;
double ref_241 = 6.05;
double err_238 = 0.1*ref_238;
double err_241 = 0.1*ref_241;
*/

double ref_238 = 10.1;
double ref_241 = 6.05;
double err_238 = 0.10*ref_238;
double err_241 = 0.10*ref_241;
 
//////////////////////////////////////////////////////////
///////////////////////// RCT /////////////////////////
//////////////////////////////////////////////////////////

const int Ndata_RCT = 25;
double sigmaf_fit_RCT[iso_N] = {0};
TMatrixD f_RCT(Ndata_RCT, iso_N);
TMatrixD sigmaf_RCT(iso_N, 1); // sigma_RCT - f_RCT * sigmaf_RCT;
TMatrixD cov_RCT(Ndata_RCT, Ndata_RCT);
TMatrixD Ratio_exp_RCT(Ndata_RCT, 1);
double a_sigmaf_SH_RCT[iso_N] = {4.40, 6.69, 10.10, 6.03};// 239, 235, 238, 241: https://arxiv.org/pdf/1608.04096.pdf
TMatrixD sigmaf_SH_RCT(iso_N, 1);

//////////////////////////////////////////////////////////
///////////////////////// DYB /////////////////////////
//////////////////////////////////////////////////////////

const int Ndata_DYB = 8;
double sigmaf_fit_DYB[iso_N] = {0};
TMatrixD sigma_DYB(Ndata_DYB, 1);
TMatrixD f_DYB(Ndata_DYB, iso_N);
TMatrixD sigmaf_DYB(iso_N, 1); // sigma_DYB - f_DYB * sigmaf_DYB;
TMatrixD cov_DYB(Ndata_DYB, Ndata_DYB);

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

  //////////////////////////////////////////////////////////
  ///////////////////////// RCT /////////////////////////
  //////////////////////////////////////////////////////////

  sigmaf_SH_RCT[0][0] = a_sigmaf_SH_RCT[0];
  sigmaf_SH_RCT[1][0] = a_sigmaf_SH_RCT[1];
  sigmaf_SH_RCT[2][0] = a_sigmaf_SH_RCT[2];
  sigmaf_SH_RCT[3][0] = a_sigmaf_SH_RCT[3];

  ifstream file_flux_RCT("./inputs/edit_table.dat");// https://arxiv.org/pdf/1703.00860.pdf
  for( int i=1; i<=Ndata_RCT; i++ ) {
    int line;
    double f239, f235, f238, f241, Ratio, self, corrA, corrB;    
    file_flux_RCT >> line >> f235 >> f238>> f239>> f241>> Ratio>> self>> corrA>> corrB;
    
    Ratio_exp_RCT[i-1][0] = Ratio;

    f_RCT[i-1][0] = f239;
    f_RCT[i-1][1] = f235;
    f_RCT[i-1][2] = f238;
    f_RCT[i-1][3] = f241;

    for(int j=1; j<=Ndata_RCT; j++) {
      if( i == j ) {
	cov_RCT[i-1][j-1] = self;
      }
      else {
	if( (i==1||i==2) && (j==1||j==2) ) cov_RCT[i-1][j-1] = corrA;

	if( (i==3||i==4) && (j==3||j==4) ) cov_RCT[i-1][j-1] = corrA;
	if( (i==5||i==6||i==7) && (j==5||j==6||j==7) ) cov_RCT[i-1][j-1] = corrA;
	if( (i==3||i==4) && (j==5||j==6||j==7) ) cov_RCT[i-1][j-1] = corrB;
	if( (i==5||i==6||i==7) && (j==3||j==4) ) cov_RCT[i-1][j-1] = corrB;

	if( (i==8||i==9||i==10) && (j==8||j==9||j==10) ) cov_RCT[i-1][j-1] = corrA;

	if( (i==11||i==12||i==13) && (j==11||j==12||j==13) ) cov_RCT[i-1][j-1] = corrA;
	if( (i==11||i==12||i==13) && (j==14) ) cov_RCT[i-1][j-1] = corrB;
	if( (j==11||j==12||j==13) && (i==14) ) cov_RCT[i-1][j-1] = corrB;

	if( (i==15||i==16) && (j==15||j==16) ) cov_RCT[i-1][j-1] = corrA;
      }
    }
  }

  //////////////////////////////////////////////////////////
  ///////////////////////// DYB /////////////////////////
  //////////////////////////////////////////////////////////

  ifstream file_flux_DYB("./inputs/flux_data.dat");
  for(int idx=0; idx<8; idx++) {
    double lowF, hghF, f239, f235, f238, f241, sigma;    
    file_flux_DYB>>lowF>>hghF>> f239>> f235>> f238>> f241>> sigma;

    sigma_DYB[idx][0] = sigma;
    f_DYB[idx][0] = f239;
    f_DYB[idx][1] = f235;
    f_DYB[idx][2] = f238;
    f_DYB[idx][3] = f241;
  }

  ifstream file_covstat_DYB("./inputs/cov_stat.dat");
  ifstream file_covsyst_DYB("./inputs/cov_syst.dat");
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

      file_covstat_DYB >> temp1;
      file_covsyst_DYB >> temp2;
      
      cov_DYB[i][j] += temp1;
      cov_DYB[i][j] += temp2 * factor_a * factor_b;

      //cout<<i<<"\t"<<j<<"\t"<<temp1<<"\t"<<temp2<<endl;
    }
  // -----------------------------------------------------

  cout<<"Fit initialization"<<endl;
  
  ROOT::Minuit2::Minuit2Minimizer min( ROOT::Minuit2::kMigrad );
  min.SetPrintLevel(2);
  min.SetStrategy(1); //0- cursory, 1- default, 2- thorough yet no more successful
  min.SetMaxFunctionCalls(500000);
  min.SetMaxIterations(500000);
  min.SetTolerance(2e-8); // tolerance*2e-3 = edm precision
  min.SetPrecision(1e-18); //precision in the target function

  ROOT::Math::Functor Chi2Functor( &fcn, 8 );
  min.SetFunction(Chi2Functor);

  int par_index = 0;
  min.SetVariable( par_index, "sigmaRCT_239", 4, 1e-4 ); par_index++;
  min.SetVariable( par_index, "sigmaRCT_235", 6, 1e-4 ); par_index++;
  min.SetVariable( par_index, "sigmaRCT_238", 11, 1e-4 ); par_index++;
  min.SetVariable( par_index, "sigmaRCT_241", 8, 1e-4 ); par_index++;

  min.SetVariable( par_index, "sigmaDYB_239", 4, 1e-4 ); par_index++;
  min.SetVariable( par_index, "sigmaDYB_235", 6, 1e-4 ); par_index++;
  min.SetVariable( par_index, "sigmaDYB_238", 11, 1e-4 ); par_index++;
  min.SetVariable( par_index, "sigmaDYB_241", 8, 1e-4 ); par_index++;
 
  min.Minimize();

  const double *par_val = min.X();
  const double *par_err = min.Errors();
  
  for( int idx=0; idx<4; idx++ ) {
    ofstream out_RCT_best("out_RCT_best.dat", ios::out | ios::app);
    out_RCT_best<<par_val[idx]<<"\t"<<par_err[idx]<<endl;
    out_RCT_best.close();
  }
  
  for( int idx=4; idx<7; idx++ ) {
    ofstream out_DYB_best("out_DYB_best.dat", ios::out | ios::app);
    out_DYB_best<<par_val[idx]<<"\t"<<par_err[idx]<<endl;
    out_DYB_best.close();
  }
  
  // -----------------------------------------------------

#ifdef Countor
  const int Npoints = 100;
  unsigned int Npoints_usr = 100;

  //////////////////////////////////////////////////////////
  ///////////////////////// RCT /////////////////////////
  //////////////////////////////////////////////////////////

  double fit_RCT_s239_1s[Npoints] = {0};
  double fit_RCT_s235_1s[Npoints] = {0};

  double fit_RCT_s239_2s[Npoints] = {0};
  double fit_RCT_s235_2s[Npoints] = {0};

  double fit_RCT_s239_3s[Npoints] = {0};
  double fit_RCT_s235_3s[Npoints] = {0};

  min.SetErrorDef( 2.2977 ); min.Contour(0,1, Npoints_usr, fit_RCT_s239_1s, fit_RCT_s235_1s);
  min.SetErrorDef( 6.2021 ); min.Contour(0,1, Npoints_usr, fit_RCT_s239_2s, fit_RCT_s235_2s);
  min.SetErrorDef(11.6182 ); min.Contour(0,1, Npoints_usr, fit_RCT_s239_3s, fit_RCT_s235_3s);

  ofstream out_RCT_1s("out_RCT_1s.dat", ios::out | ios::app);
  ofstream out_RCT_2s("out_RCT_2s.dat", ios::out | ios::app);
  ofstream out_RCT_3s("out_RCT_3s.dat", ios::out | ios::app);
  for(int idx=0; idx<Npoints; idx++) {
    out_RCT_1s<<fit_RCT_s239_1s[idx]<<"\t"<<fit_RCT_s235_1s[idx]<<endl;
    out_RCT_2s<<fit_RCT_s239_2s[idx]<<"\t"<<fit_RCT_s235_2s[idx]<<endl;
    out_RCT_3s<<fit_RCT_s239_3s[idx]<<"\t"<<fit_RCT_s235_3s[idx]<<endl;
  }

  out_RCT_1s.close();
  out_RCT_2s.close();
  out_RCT_3s.close();

  //////////////////////////////////////////////////////////
  ///////////////////////// DYB /////////////////////////
  //////////////////////////////////////////////////////////

  double fit_DYB_s239_1s[Npoints] = {0};
  double fit_DYB_s235_1s[Npoints] = {0};

  double fit_DYB_s239_2s[Npoints] = {0};
  double fit_DYB_s235_2s[Npoints] = {0};

  double fit_DYB_s239_3s[Npoints] = {0};
  double fit_DYB_s235_3s[Npoints] = {0};

  min.SetErrorDef( 2.2977 ); min.Contour(4,5, Npoints_usr, fit_DYB_s239_1s, fit_DYB_s235_1s);
  min.SetErrorDef( 6.2021 ); min.Contour(4,5, Npoints_usr, fit_DYB_s239_2s, fit_DYB_s235_2s);
  min.SetErrorDef(11.6182 ); min.Contour(4,5, Npoints_usr, fit_DYB_s239_3s, fit_DYB_s235_3s);

  ofstream out_DYB_1s("out_DYB_1s.dat", ios::out | ios::app);
  ofstream out_DYB_2s("out_DYB_2s.dat", ios::out | ios::app);
  ofstream out_DYB_3s("out_DYB_3s.dat", ios::out | ios::app);
  for(int idx=0; idx<Npoints; idx++) {
    out_DYB_1s<<fit_DYB_s239_1s[idx]<<"\t"<<fit_DYB_s235_1s[idx]<<endl;
    out_DYB_2s<<fit_DYB_s239_2s[idx]<<"\t"<<fit_DYB_s235_2s[idx]<<endl;
    out_DYB_3s<<fit_DYB_s239_3s[idx]<<"\t"<<fit_DYB_s235_3s[idx]<<endl;
  }

  out_DYB_1s.close();
  out_DYB_2s.close();
  out_DYB_3s.close();
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
  cycle++;
  if( cycle%100==0 ) cout<<TString::Format(" ---> cycle %6d", cycle)<<endl;

  //////////////////////////////////////////////////////////
  ///////////////////////// RCT /////////////////////////
  //////////////////////////////////////////////////////////

  double chi2_tot_RCT = 0;
  double chi2_main_RCT = 0;
  double chi2_pull_RCT = 0;
  
  int index_par = 0;
  sigmaf_fit_RCT[index_par] = par[index_par]; index_par++;
  sigmaf_fit_RCT[index_par] = par[index_par]; index_par++;
  sigmaf_fit_RCT[index_par] = par[index_par]; index_par++;
  sigmaf_fit_RCT[index_par] = par[index_par]; index_par++;

  for( int idx=0; idx<iso_N; idx++ ) {
    sigmaf_RCT[idx][0] = sigmaf_fit_RCT[idx];
  }
 
  TMatrixD matrix_sigma_RCT = f_RCT * sigmaf_RCT;
  TMatrixD matrix_sigma_SH_RCT =  f_RCT * sigmaf_SH_RCT;
  TMatrixD matrix_sigma_exp_RCT(Ndata_RCT, 1);

  for( int idx=0; idx<Ndata_RCT; idx++ ) {
    matrix_sigma_exp_RCT[idx][0] = Ratio_exp_RCT[idx][0] * matrix_sigma_SH_RCT[idx][0];
  }
  
  TMatrixD matrix_diff_RCT = matrix_sigma_RCT - matrix_sigma_exp_RCT;
  TMatrixD matrixT_diff_RCT(1, Ndata_RCT);
  matrixT_diff_RCT.Transpose(matrix_diff_RCT);

  TMatrixD edit_cov_RCT = cov_RCT;
  for( int i=0; i<Ndata_RCT; i++ ) {
    for( int j=0; j<Ndata_RCT; j++ ) {
      edit_cov_RCT[i][j] = edit_cov_RCT[i][j] * edit_cov_RCT[i][j] * 0.01 * 0.01;
      edit_cov_RCT[i][j] = edit_cov_RCT[i][j] * matrix_sigma_exp_RCT[i][0] * matrix_sigma_exp_RCT[j][0];
    }
  }

  TMatrixD edit_cov_invert_RCT = edit_cov_RCT;
  edit_cov_invert_RCT.Invert();

  TMatrixD matrix_main_RCT = matrixT_diff_RCT * edit_cov_invert_RCT * matrix_diff_RCT;
  chi2_main_RCT = matrix_main_RCT[0][0];

  chi2_pull_RCT += pow(sigmaf_fit_RCT[2]-ref_238, 2)/err_238/err_238;
  chi2_pull_RCT += pow(sigmaf_fit_RCT[3]-ref_241, 2)/err_241/err_241;

  chi2_tot_RCT = chi2_main_RCT + chi2_pull_RCT; 

  //////////////////////////////////////////////////////////
  ///////////////////////// DYB /////////////////////////
  //////////////////////////////////////////////////////////

  double chi2_tot_DYB = 0;
  double chi2_main_DYB = 0;
  double chi2_pull_DYB = 0;
  
  //int index_par = 0;
  sigmaf_fit_DYB[index_par-4] = par[index_par]; index_par++;
  sigmaf_fit_DYB[index_par-4] = par[index_par]; index_par++;
  sigmaf_fit_DYB[index_par-4] = par[index_par]; index_par++;
  sigmaf_fit_DYB[index_par-4] = par[index_par]; index_par++;

  for( int idx=0; idx<iso_N; idx++ ) {
    sigmaf_DYB[idx][0] = sigmaf_fit_DYB[idx];
  }
  TMatrixD matrix_diff_DYB = sigma_DYB - f_DYB * sigmaf_DYB;
  TMatrixD matrixT_diff_DYB(1, Ndata_DYB);
  matrixT_diff_DYB.Transpose(matrix_diff_DYB);
  TMatrixD cov_invert_DYB = cov_DYB;
  cov_invert_DYB.Invert();
  TMatrixD matrix_main_DYB = matrixT_diff_DYB * cov_invert_DYB * matrix_diff_DYB;
  chi2_main_DYB = matrix_main_DYB[0][0];

  chi2_pull_DYB += pow(sigmaf_fit_DYB[2]-ref_238, 2)/err_238/err_238;
  chi2_pull_DYB += pow(sigmaf_fit_DYB[3]-ref_241, 2)/err_241/err_241;

  chi2_tot_DYB = chi2_main_DYB + chi2_pull_DYB;
  
  // -----------------------------------------------------

  return chi2_tot_RCT + chi2_tot_DYB;
}

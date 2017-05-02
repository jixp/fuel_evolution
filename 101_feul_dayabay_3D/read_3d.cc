#include<iostream>
#include<fstream>
using namespace std;

#include "TString.h"
#include "TMatrixD.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraph2D.h"
#include "TGraph.h"
// ----------------------------------------------------- Global variables

const int Ndata_DYB = 8;
const int iso_N = 4;

double sigmaf_fit_DYB[iso_N] = {0};
TMatrixD sigma_DYB(Ndata_DYB, 1);
TMatrixD f_DYB(Ndata_DYB, iso_N);
TMatrixD sigmaf_DYB(iso_N, 1); // sigma_DYB - f_DYB * sigmaf_DYB;
TMatrixD cov_DYB(Ndata_DYB, Ndata_DYB);

// 
double fcn(const double *par);

void read_3d()
{
  TString roostr = "";

  // -----------------------------------------------------
  
  ifstream file_flux_DYB("./inputs/flux_data.dat");
  for(int idx=0; idx<8; idx++) {
    double lowF, hghF, f239, f235, f238, f241, sigma;    
    file_flux_DYB >> lowF >> hghF >>  f239 >>  f235 >>  f238 >>  f241 >>  sigma;

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

      file_covstat_DYB  >>  temp1;
      file_covsyst_DYB  >>  temp2;
      
      cov_DYB[i][j] += temp1;
      cov_DYB[i][j] += temp2 * factor_a * factor_b;

      //cout<<i<<"\t"<<j<<"\t"<<temp1<<"\t"<<temp2<<endl;
    }

  //double a0_s239 = 4.2541 - 3.0*0.2561;
  //double an_s239 = 4.2541 + 3.0*0.2561;
  double a0_s239 = 3.5;
  double an_s239 = 5.0;
  int Nstep_s239 = 80;
  double step_s239 = (an_s239-a0_s239)/Nstep_s239;

  //double a0_s235 = 6.1758 - 3.0*0.1739;
  //double an_s235 = 6.1758 + 3.0*0.1739;;
  double a0_s235 = 5.2;
  double an_s235 = 7.2;
  int Nstep_s235 = 100;
  double step_s235 = (an_s235-a0_s235)/Nstep_s235;

  //double a0_s238 =  10.1168 - 3;
  //double an_s238 =  10.1168 + 3;
  double a0_s238 =  5;
  double an_s238 =  15;
  int Nstep_s238 = 100;
  double step_s238 = (an_s238-a0_s238)/Nstep_s238;

  //double a0_s241 = 6.1195 - 3;
  //double an_s241 = 6.1195 + 3;;
  double a0_s241 = 2;
  double an_s241 = 10;
  int Nstep_s241 = 100;
  double step_s241 = (an_s241-a0_s241)/Nstep_s241;

  double sigma_par[4] = {0};

  TGraph2D *gh2 = new TGraph2D();
  int gh2_line = 0;

  for( int i=1; i<=Nstep_s239; i++ ) {
    cout<<TString::Format(" ---> procesing %5d, total %5d", i, Nstep_s239)<<endl;
    for( int j=1; j<=Nstep_s235; j++ ) {
      for( int m=1; m<=Nstep_s238; m++ ) {
	for( int n=1; n<=Nstep_s241; n++ ) {
	  double val_s239 = a0_s239 + (i-1) * step_s239;
	  double val_s235 = a0_s235 + (j-1) * step_s235;
	  double val_s238 = a0_s238 + (m-1) * step_s238;
	  double val_s241 = a0_s241 + (n-1) * step_s241;
	  
	  double fbar_238 = 0.076;
	  double fbar_241 = 0.054;
	  double eff_factor = ( val_s238*fbar_238 + val_s241*fbar_241 )/( fbar_238 + fbar_241 );

	  sigma_par[0] = val_s239;
	  sigma_par[1] = val_s235;
	  sigma_par[2] = val_s238;
	  sigma_par[3] = val_s241;

	  double chi2_val = fcn( sigma_par );
	  double chi2_min = 3.7092;
	  if( chi2_val - chi2_min < 2.2977 ) {
	    gh2_line++;
	    gh2->SetPoint( gh2_line-1, val_s239, val_s235, eff_factor );
	  }

	}
      }
    }
  }
  
  TGraph2D *ghR = new TGraph2D();
  int ghR_line = 0;

  TGraph *ghC = new TGraph();
  int ghC_line = 0;

  for( int i=1; i<=1000; i++ )
    {
      double sbar_239 = 4.36;// Daya Bay paper
      double sbar_235 = 6.69;// Daya Bay paper
      double sbar_238 = 10.10;
      double sbar_241 = 6.03;

      double val_s239 =  sbar_239* ( 1 - i*1./1000 );
      double val_s235 =  sbar_235* ( 1 - i*1./1000 );
      double val_s238 =  sbar_238* ( 1 - i*1./1000 );
      double val_s241 =  sbar_241* ( 1 - i*1./1000 );
      double xxx = 1 - i*1./1000;

      double fbar_238 = 0.076;
      double fbar_241 = 0.054;
      double eff_factor = ( val_s238*fbar_238 + val_s241*fbar_241 )/( fbar_238 + fbar_241 );

      sigma_par[0] = val_s239;
      sigma_par[1] = val_s235;
      sigma_par[2] = val_s238;
      sigma_par[3] = val_s241;
      
      ghR_line++;
      ghR->SetPoint( ghR_line-1, val_s239, val_s235, eff_factor );

      double chi2_val = fcn( sigma_par );
      double chi2_min = 3.7092;
      ghC_line++;
      ghC->SetPoint( ghC_line-1, xxx,  chi2_val-chi2_min );
    }




  TCanvas *canv_gh2 = new TCanvas("canv_gh2", "canv_gh2", 800, 600);
  TH3D *h3_gh2 = new TH3D("h3_gh2", "h3_gh2", 100, 3.0, 5.5, 100, 5.2, 7.2, 100, 6, 12);
  //TH3D *h3_gh2 = new TH3D("h3_gh2", "h3_gh2", 100, a0_s239, an_s239, 100, a0_s235, an_s235, 100, 6, 12);
  h3_gh2->Draw();
  h3_gh2->SetStats(0);
  h3_gh2->SetXTitle(" #sigma_{239}");
  h3_gh2->SetYTitle(" #sigma_{235}");
  h3_gh2->SetZTitle(" #sigma_{eff}");

  gh2->Draw("p same");
  gh2->SetMarkerColor(kBlue);

  ghR->Draw("p same");
  ghR->SetMarkerStyle(7);

  canv_gh2->SaveAs("canv_gh2.root");


  TCanvas *canv_ghC = new TCanvas("canv_ghC", "canv_ghC", 800, 600);
  TH2D *h2_ghC = new TH2D("", "", 100, 0, 1, 100, 0, 20);
  h2_ghC->Draw();
  h2_ghC->SetStats(0);
  ghC->Draw("pl same");
  canv_ghC->SaveAs("canv_ghC.root");

}

// ----------------------------------------------------- ----------------------------------------------------- 

double fcn(const double *par)
{
  double chi2_tot_DYB = 0;
  double chi2_main_DYB = 0;
  double chi2_pull_DYB = 0;
 
  double ref_238 = 10.1;
  double ref_241 = 6.05;
  double err_238 = 0.1*ref_238;
  double err_241 = 0.1*ref_241;
  

  int index_par = 0;
  sigmaf_fit_DYB[index_par] = par[index_par]; index_par++;
  sigmaf_fit_DYB[index_par] = par[index_par]; index_par++;
  sigmaf_fit_DYB[index_par] = par[index_par]; index_par++;
  sigmaf_fit_DYB[index_par] = par[index_par]; index_par++;

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

  return chi2_tot_DYB;
}

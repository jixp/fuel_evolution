void read_contour()
{
  gROOT->ProcessLine(".x ./lhcbStyle.C");

  int color_BB_1sigma = kRed;
  int color_BB_2sigma = kRed;
  int color_BB_3sigma = kRed;

  int color_AA_1sigma = kBlue;
  int color_AA_2sigma = kBlue;
  int color_AA_3sigma = kBlue;

  const int Npoints = 101;





  double s239_BB_best    =  0;
  double s235_BB_best   = 0;
  double s239_BB_err  = 0;
  double s235_BB_err = 0;

  ifstream ReadFile_BB_best("out_BB_best.dat");
  for(int idx=1; idx<=1; idx++) {
    ReadFile_BB_best>>s239_BB_best>>s239_BB_err>>s235_BB_best>>s235_BB_err;
    cout<<s239_BB_best<<"\t"<<s239_BB_err<<"\t"<<s235_BB_best<<"\t"<<s235_BB_err<<endl;
  }
  
  /// 1sigma
  double s239_BB_1sigma[Npoints]   = {0};
  double s235_BB_1sigma[Npoints] = {0};
  ifstream ReadFile_BB_1sigma("out_BB_1s.dat");
  for(int idx=1; idx<=Npoints-1; idx++)
    {
      double aa, bb;
      ReadFile_BB_1sigma>>aa>>bb;

      s239_BB_1sigma[idx-1] = aa;
      s235_BB_1sigma[idx-1] = bb;
    }
  s239_BB_1sigma[Npoints-1]  = s239_BB_1sigma[0];    //cycle
  s235_BB_1sigma[Npoints-1] = s235_BB_1sigma[0];
  TGraph *gh_BB_1sigma = new TGraph(Npoints, s235_BB_1sigma, s239_BB_1sigma);
  
  /// 2sigma
  double s239_BB_2sigma[Npoints]   = {0};
  double s235_BB_2sigma[Npoints] = {0};
  ifstream ReadFile_BB_2sigma("out_BB_2s.dat");
  for(int idx=1; idx<=Npoints-1; idx++)
    {
      double aa, bb;
      ReadFile_BB_2sigma>>aa>>bb;

      s239_BB_2sigma[idx-1] = aa;
      s235_BB_2sigma[idx-1] = bb;
    }
  s239_BB_2sigma[Npoints-1]  = s239_BB_2sigma[0];    //cycle
  s235_BB_2sigma[Npoints-1] = s235_BB_2sigma[0];
  TGraph *gh_BB_2sigma = new TGraph(Npoints, s235_BB_2sigma, s239_BB_2sigma);
 
  /// 3sigma
  double s239_BB_3sigma[Npoints]   = {0};
  double s235_BB_3sigma[Npoints] = {0};
  ifstream ReadFile_BB_3sigma("out_BB_3s.dat");
  for(int idx=1; idx<=Npoints-1; idx++)
    {
      double aa, bb;
      ReadFile_BB_3sigma>>aa>>bb;

      s239_BB_3sigma[idx-1] = aa;
      s235_BB_3sigma[idx-1] = bb;
    }
  s239_BB_3sigma[Npoints-1]  = s239_BB_3sigma[0];    //cycle
  s235_BB_3sigma[Npoints-1] = s235_BB_3sigma[0];
  TGraph *gh_BB_3sigma = new TGraph(Npoints, s235_BB_3sigma, s239_BB_3sigma);
 

  /// best fit value
  TGraphErrors *gh_BB_best = new TGraphErrors(1, &s235_BB_best, &s239_BB_best, 0, 0);
 








  double s239_AA_best    =  0;
  double s235_AA_best   = 0;
  double s239_AA_err  = 0;
  double s235_AA_err = 0;

  ifstream ReadFile_AA_best("out_AA_best.dat");
  for(int idx=1; idx<=1; idx++) {
    ReadFile_AA_best>>s239_AA_best>>s239_AA_err>>s235_AA_best>>s235_AA_err;
    cout<<s239_AA_best<<"\t"<<s239_AA_err<<"\t"<<s235_AA_best<<"\t"<<s235_AA_err<<endl;
  }
  
  /// 1sigma
  double s239_AA_1sigma[Npoints]   = {0};
  double s235_AA_1sigma[Npoints] = {0};
  ifstream ReadFile_AA_1sigma("out_AA_1s.dat");
  for(int idx=1; idx<=Npoints-1; idx++)
    {
      double aa, bb;
      ReadFile_AA_1sigma>>aa>>bb;

      s239_AA_1sigma[idx-1] = aa;
      s235_AA_1sigma[idx-1] = bb;
    }
  s239_AA_1sigma[Npoints-1]  = s239_AA_1sigma[0];    //cycle
  s235_AA_1sigma[Npoints-1] = s235_AA_1sigma[0];
  TGraph *gh_AA_1sigma = new TGraph(Npoints, s235_AA_1sigma, s239_AA_1sigma);
  
  /// 2sigma
  double s239_AA_2sigma[Npoints]   = {0};
  double s235_AA_2sigma[Npoints] = {0};
  ifstream ReadFile_AA_2sigma("out_AA_2s.dat");
  for(int idx=1; idx<=Npoints-1; idx++)
    {
      double aa, bb;
      ReadFile_AA_2sigma>>aa>>bb;

      s239_AA_2sigma[idx-1] = aa;
      s235_AA_2sigma[idx-1] = bb;
    }
  s239_AA_2sigma[Npoints-1]  = s239_AA_2sigma[0];    //cycle
  s235_AA_2sigma[Npoints-1] = s235_AA_2sigma[0];
  TGraph *gh_AA_2sigma = new TGraph(Npoints, s235_AA_2sigma, s239_AA_2sigma);
 
  /// 3sigma
  double s239_AA_3sigma[Npoints]   = {0};
  double s235_AA_3sigma[Npoints] = {0};
  ifstream ReadFile_AA_3sigma("out_AA_3s.dat");
  for(int idx=1; idx<=Npoints-1; idx++)
    {
      double aa, bb;
      ReadFile_AA_3sigma>>aa>>bb;

      s239_AA_3sigma[idx-1] = aa;
      s235_AA_3sigma[idx-1] = bb;
    }
  s239_AA_3sigma[Npoints-1]  = s239_AA_3sigma[0];    //cycle
  s235_AA_3sigma[Npoints-1] = s235_AA_3sigma[0];
  TGraph *gh_AA_3sigma = new TGraph(Npoints, s235_AA_3sigma, s239_AA_3sigma);
 

  /// best fit value
  TGraphErrors *gh_AA_best = new TGraphErrors(1, &s235_AA_best, &s239_AA_best, 0, 0);









  //////
  TCanvas *canv_BB_contour = new TCanvas("canv_BB_contour", "canv_BB_contour", 800, 600);
  //canv_BB_contour->SetGridx();
  //canv_BB_contour->SetGridy();
  canv_BB_contour->SetLeftMargin(0.15);
  canv_BB_contour->SetBottomMargin(0.15);

  ///
  TH2D *h2_BB_basic = new TH2D("h1_BB_basic", "",  100, 5.4, 7.2, 100, 2.2, 5.5);
  h2_BB_basic->SetStats(0);
  h2_BB_basic->Draw();
  h2_BB_basic->SetXTitle("#sigma_{235}");
  h2_BB_basic->SetYTitle("#sigma_{239}");

  ///
  gh_BB_1sigma->Draw("l same");
  gh_BB_1sigma->SetLineColor(color_BB_1sigma);
  gh_BB_1sigma->SetLineWidth(2);
  gh_BB_1sigma->SetMarkerColor(color_BB_1sigma);

  ///
  gh_BB_2sigma->Draw("l same");
  gh_BB_2sigma->SetLineColor(color_BB_2sigma);
  gh_BB_2sigma->SetLineWidth(2);
  gh_BB_2sigma->SetLineStyle(7);
  gh_BB_2sigma->SetMarkerColor(color_BB_2sigma);

  ///
  gh_BB_3sigma->Draw("l same");
  gh_BB_3sigma->SetLineColor(color_BB_3sigma);
  gh_BB_3sigma->SetLineWidth(2);
  gh_BB_3sigma->SetLineStyle(3);
  gh_BB_3sigma->SetMarkerColor(color_BB_3sigma);

  ///
  gh_BB_best->Draw("p same");
  gh_BB_best->SetMarkerStyle(34);
  gh_BB_best->SetMarkerSize(1.5);
  gh_BB_best->SetMarkerColor(color_BB_3sigma);










  ///
  gh_AA_1sigma->Draw("l same");
  gh_AA_1sigma->SetLineColor(color_AA_1sigma);
  gh_AA_1sigma->SetLineWidth(2);
  gh_AA_1sigma->SetMarkerColor(color_AA_1sigma);

  ///
  gh_AA_2sigma->Draw("l same");
  gh_AA_2sigma->SetLineColor(color_AA_2sigma);
  gh_AA_2sigma->SetLineWidth(2);
  gh_AA_2sigma->SetLineStyle(7);
  gh_AA_2sigma->SetMarkerColor(color_AA_2sigma);

  ///
  gh_AA_3sigma->Draw("l same");
  gh_AA_3sigma->SetLineColor(color_AA_3sigma);
  gh_AA_3sigma->SetLineWidth(2);
  gh_AA_3sigma->SetLineStyle(3);
  gh_AA_3sigma->SetMarkerColor(color_AA_3sigma);

  ///
  gh_AA_best->Draw("p same");
  gh_AA_best->SetMarkerStyle(34);
  gh_AA_best->SetMarkerSize(1.5);
  gh_AA_best->SetMarkerColor(color_AA_3sigma);







  ///
  canv_BB_contour->SaveAs("canv_BB_contour.png");
}

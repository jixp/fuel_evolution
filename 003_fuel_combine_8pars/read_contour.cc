void read_contour()
{
  gROOT->ProcessLine(".x ./lhcbStyle.C");

  int color_RCT_1sigma = kRed;
  int color_RCT_2sigma = kRed;
  int color_RCT_3sigma = kRed;

  int color_DYB_1sigma = kBlue;
  int color_DYB_2sigma = kBlue;
  int color_DYB_3sigma = kBlue;

  const int Npoints = 101;





  double s239_RCT_best    =  0;
  double s235_RCT_best   = 0;
  double s239_RCT_err  = 0;
  double s235_RCT_err = 0;

  ifstream ReadFile_RCT_best("out_RCT_best.dat");
  for(int idx=1; idx<=1; idx++) {
    ReadFile_RCT_best>>s239_RCT_best>>s239_RCT_err>>s235_RCT_best>>s235_RCT_err;
    cout<<s239_RCT_best<<"\t"<<s239_RCT_err<<"\t"<<s235_RCT_best<<"\t"<<s235_RCT_err<<endl;
  }
  
  /// 1sigma
  double s239_RCT_1sigma[Npoints]   = {0};
  double s235_RCT_1sigma[Npoints] = {0};
  ifstream ReadFile_RCT_1sigma("out_RCT_1s.dat");
  for(int idx=1; idx<=Npoints-1; idx++)
    {
      double aa, bb;
      ReadFile_RCT_1sigma>>aa>>bb;

      s239_RCT_1sigma[idx-1] = aa;
      s235_RCT_1sigma[idx-1] = bb;
    }
  s239_RCT_1sigma[Npoints-1]  = s239_RCT_1sigma[0];    //cycle
  s235_RCT_1sigma[Npoints-1] = s235_RCT_1sigma[0];
  TGraph *gh_RCT_1sigma = new TGraph(Npoints, s235_RCT_1sigma, s239_RCT_1sigma);
  
  /// 2sigma
  double s239_RCT_2sigma[Npoints]   = {0};
  double s235_RCT_2sigma[Npoints] = {0};
  ifstream ReadFile_RCT_2sigma("out_RCT_2s.dat");
  for(int idx=1; idx<=Npoints-1; idx++)
    {
      double aa, bb;
      ReadFile_RCT_2sigma>>aa>>bb;

      s239_RCT_2sigma[idx-1] = aa;
      s235_RCT_2sigma[idx-1] = bb;
    }
  s239_RCT_2sigma[Npoints-1]  = s239_RCT_2sigma[0];    //cycle
  s235_RCT_2sigma[Npoints-1] = s235_RCT_2sigma[0];
  TGraph *gh_RCT_2sigma = new TGraph(Npoints, s235_RCT_2sigma, s239_RCT_2sigma);
 
  /// 3sigma
  double s239_RCT_3sigma[Npoints]   = {0};
  double s235_RCT_3sigma[Npoints] = {0};
  ifstream ReadFile_RCT_3sigma("out_RCT_3s.dat");
  for(int idx=1; idx<=Npoints-1; idx++)
    {
      double aa, bb;
      ReadFile_RCT_3sigma>>aa>>bb;

      s239_RCT_3sigma[idx-1] = aa;
      s235_RCT_3sigma[idx-1] = bb;
    }
  s239_RCT_3sigma[Npoints-1]  = s239_RCT_3sigma[0];    //cycle
  s235_RCT_3sigma[Npoints-1] = s235_RCT_3sigma[0];
  TGraph *gh_RCT_3sigma = new TGraph(Npoints, s235_RCT_3sigma, s239_RCT_3sigma);
 

  /// best fit value
  TGraphErrors *gh_RCT_best = new TGraphErrors(1, &s235_RCT_best, &s239_RCT_best, 0, 0);
 








  double s239_DYB_best    =  0;
  double s235_DYB_best   = 0;
  double s239_DYB_err  = 0;
  double s235_DYB_err = 0;

  ifstream ReadFile_DYB_best("out_DYB_best.dat");
  for(int idx=1; idx<=1; idx++) {
    ReadFile_DYB_best>>s239_DYB_best>>s239_DYB_err>>s235_DYB_best>>s235_DYB_err;
    cout<<s239_DYB_best<<"\t"<<s239_DYB_err<<"\t"<<s235_DYB_best<<"\t"<<s235_DYB_err<<endl;
  }
  
  /// 1sigma
  double s239_DYB_1sigma[Npoints]   = {0};
  double s235_DYB_1sigma[Npoints] = {0};
  ifstream ReadFile_DYB_1sigma("out_DYB_1s.dat");
  for(int idx=1; idx<=Npoints-1; idx++)
    {
      double aa, bb;
      ReadFile_DYB_1sigma>>aa>>bb;

      s239_DYB_1sigma[idx-1] = aa;
      s235_DYB_1sigma[idx-1] = bb;
    }
  s239_DYB_1sigma[Npoints-1]  = s239_DYB_1sigma[0];    //cycle
  s235_DYB_1sigma[Npoints-1] = s235_DYB_1sigma[0];
  TGraph *gh_DYB_1sigma = new TGraph(Npoints, s235_DYB_1sigma, s239_DYB_1sigma);
  
  /// 2sigma
  double s239_DYB_2sigma[Npoints]   = {0};
  double s235_DYB_2sigma[Npoints] = {0};
  ifstream ReadFile_DYB_2sigma("out_DYB_2s.dat");
  for(int idx=1; idx<=Npoints-1; idx++)
    {
      double aa, bb;
      ReadFile_DYB_2sigma>>aa>>bb;

      s239_DYB_2sigma[idx-1] = aa;
      s235_DYB_2sigma[idx-1] = bb;
    }
  s239_DYB_2sigma[Npoints-1]  = s239_DYB_2sigma[0];    //cycle
  s235_DYB_2sigma[Npoints-1] = s235_DYB_2sigma[0];
  TGraph *gh_DYB_2sigma = new TGraph(Npoints, s235_DYB_2sigma, s239_DYB_2sigma);
 
  /// 3sigma
  double s239_DYB_3sigma[Npoints]   = {0};
  double s235_DYB_3sigma[Npoints] = {0};
  ifstream ReadFile_DYB_3sigma("out_DYB_3s.dat");
  for(int idx=1; idx<=Npoints-1; idx++)
    {
      double aa, bb;
      ReadFile_DYB_3sigma>>aa>>bb;

      s239_DYB_3sigma[idx-1] = aa;
      s235_DYB_3sigma[idx-1] = bb;
    }
  s239_DYB_3sigma[Npoints-1]  = s239_DYB_3sigma[0];    //cycle
  s235_DYB_3sigma[Npoints-1] = s235_DYB_3sigma[0];
  TGraph *gh_DYB_3sigma = new TGraph(Npoints, s235_DYB_3sigma, s239_DYB_3sigma);
 

  /// best fit value
  TGraphErrors *gh_DYB_best = new TGraphErrors(1, &s235_DYB_best, &s239_DYB_best, 0, 0);









  //////
  TCanvas *canv_RCT_contour = new TCanvas("canv_RCT_contour", "canv_RCT_contour", 800, 600);
  //canv_RCT_contour->SetGridx();
  //canv_RCT_contour->SetGridy();
  canv_RCT_contour->SetLeftMargin(0.15);
  canv_RCT_contour->SetBottomMargin(0.15);

  ///
  TH2D *h2_RCT_basic = new TH2D("h1_RCT_basic", "",  100, 5.4, 7.2, 100, 2.2, 5.5);
  h2_RCT_basic->SetStats(0);
  h2_RCT_basic->Draw();
  h2_RCT_basic->SetXTitle("#sigma_{235}");
  h2_RCT_basic->SetYTitle("#sigma_{239}");

  ///
  gh_RCT_1sigma->Draw("l same");
  gh_RCT_1sigma->SetLineColor(color_RCT_1sigma);
  gh_RCT_1sigma->SetLineWidth(2);
  gh_RCT_1sigma->SetMarkerColor(color_RCT_1sigma);

  ///
  gh_RCT_2sigma->Draw("l same");
  gh_RCT_2sigma->SetLineColor(color_RCT_2sigma);
  gh_RCT_2sigma->SetLineWidth(2);
  gh_RCT_2sigma->SetLineStyle(7);
  gh_RCT_2sigma->SetMarkerColor(color_RCT_2sigma);

  ///
  gh_RCT_3sigma->Draw("l same");
  gh_RCT_3sigma->SetLineColor(color_RCT_3sigma);
  gh_RCT_3sigma->SetLineWidth(2);
  gh_RCT_3sigma->SetLineStyle(3);
  gh_RCT_3sigma->SetMarkerColor(color_RCT_3sigma);

  ///
  gh_RCT_best->Draw("p same");
  gh_RCT_best->SetMarkerStyle(34);
  gh_RCT_best->SetMarkerSize(1.5);
  gh_RCT_best->SetMarkerColor(color_RCT_3sigma);










  ///
  gh_DYB_1sigma->Draw("l same");
  gh_DYB_1sigma->SetLineColor(color_DYB_1sigma);
  gh_DYB_1sigma->SetLineWidth(2);
  gh_DYB_1sigma->SetMarkerColor(color_DYB_1sigma);

  ///
  gh_DYB_2sigma->Draw("l same");
  gh_DYB_2sigma->SetLineColor(color_DYB_2sigma);
  gh_DYB_2sigma->SetLineWidth(2);
  gh_DYB_2sigma->SetLineStyle(7);
  gh_DYB_2sigma->SetMarkerColor(color_DYB_2sigma);

  ///
  gh_DYB_3sigma->Draw("l same");
  gh_DYB_3sigma->SetLineColor(color_DYB_3sigma);
  gh_DYB_3sigma->SetLineWidth(2);
  gh_DYB_3sigma->SetLineStyle(3);
  gh_DYB_3sigma->SetMarkerColor(color_DYB_3sigma);

  ///
  gh_DYB_best->Draw("p same");
  gh_DYB_best->SetMarkerStyle(34);
  gh_DYB_best->SetMarkerSize(1.5);
  gh_DYB_best->SetMarkerColor(color_DYB_3sigma);







  ///
  canv_RCT_contour->SaveAs("canv_RCT_contour.png");
}

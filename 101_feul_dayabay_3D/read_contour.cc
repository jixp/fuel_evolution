void read_contour()
{
  gROOT->ProcessLine(".x ./lhcbStyle.C");

  int color_1sigma = kRed;
  int color_2sigma = kBlue;
  int color_3sigma = kGreen+1;

  const int Npoints = 101;

  double s239_best    =  0;
  double s235_best   = 0;
  double s239_err  = 0;
  double s235_err = 0;

  ifstream ReadFile_best("out_DYB_best.dat");
  for(int idx=1; idx<=1; idx++) {
    ReadFile_best>>s239_best>>s239_err>>s235_best>>s235_err;
    cout<<s239_best<<"\t"<<s239_err<<"\t"<<s235_best<<"\t"<<s235_err<<endl;
  }
  
  /// 1sigma
  double s239_1sigma[Npoints]   = {0};
  double s235_1sigma[Npoints] = {0};
  ifstream ReadFile_1sigma("out_DYB_1s.dat");
  for(int idx=1; idx<=Npoints-1; idx++)
    {
      double aa, bb;
      ReadFile_1sigma>>aa>>bb;

      s239_1sigma[idx-1] = aa;
      s235_1sigma[idx-1] = bb;
    }
  s239_1sigma[Npoints-1]  = s239_1sigma[0];    //cycle
  s235_1sigma[Npoints-1] = s235_1sigma[0];
  TGraph *gh_1sigma = new TGraph(Npoints, s235_1sigma, s239_1sigma);
  
  /// 2sigma
  double s239_2sigma[Npoints]   = {0};
  double s235_2sigma[Npoints] = {0};
  ifstream ReadFile_2sigma("out_DYB_2s.dat");
  for(int idx=1; idx<=Npoints-1; idx++)
    {
      double aa, bb;
      ReadFile_2sigma>>aa>>bb;

      s239_2sigma[idx-1] = aa;
      s235_2sigma[idx-1] = bb;
    }
  s239_2sigma[Npoints-1]  = s239_2sigma[0];    //cycle
  s235_2sigma[Npoints-1] = s235_2sigma[0];
  TGraph *gh_2sigma = new TGraph(Npoints, s235_2sigma, s239_2sigma);
 
  /// 3sigma
  double s239_3sigma[Npoints]   = {0};
  double s235_3sigma[Npoints] = {0};
  ifstream ReadFile_3sigma("out_DYB_3s.dat");
  for(int idx=1; idx<=Npoints-1; idx++)
    {
      double aa, bb;
      ReadFile_3sigma>>aa>>bb;

      s239_3sigma[idx-1] = aa;
      s235_3sigma[idx-1] = bb;
    }
  s239_3sigma[Npoints-1]  = s239_3sigma[0];    //cycle
  s235_3sigma[Npoints-1] = s235_3sigma[0];
  TGraph *gh_3sigma = new TGraph(Npoints, s235_3sigma, s239_3sigma);
 

  /// best fit value
  TGraphErrors *gh_best = new TGraphErrors(1, &s235_best, &s239_best, &s235_err, &s239_err);
 
  //////
  TCanvas *canv_contour = new TCanvas("canv_contour", "canv_contour", 800, 600);
  canv_contour->SetGridx();
  canv_contour->SetGridy();
  canv_contour->SetLeftMargin(0.15);
  canv_contour->SetBottomMargin(0.15);

  ///
  TH2D *h2_basic = new TH2D("h1_basic", "",  100, 5.2, 7.2, 100, 3.0, 5.5);
  h2_basic->SetStats(0);
  h2_basic->Draw();
  h2_basic->SetXTitle("#sigma_{235}");
  h2_basic->SetYTitle("#sigma_{239}");

  ///
  gh_1sigma->Draw("l same");
  gh_1sigma->SetLineColor(color_1sigma);
  gh_1sigma->SetLineWidth(2);
  gh_1sigma->SetMarkerColor(color_1sigma);

  ///
  gh_2sigma->Draw("l same");
  gh_2sigma->SetLineColor(color_2sigma);
  gh_2sigma->SetLineWidth(2);
  gh_2sigma->SetMarkerColor(color_2sigma);

  ///
  gh_3sigma->Draw("l same");
  gh_3sigma->SetLineColor(color_3sigma);
  gh_3sigma->SetLineWidth(2);
  gh_3sigma->SetMarkerColor(color_3sigma);

  ///
  gh_best->Draw("p same");
  gh_best->SetMarkerStyle(1);
  gh_best->SetMarkerSize(1.3);

  ///
  canv_contour->SaveAs("canv_contour.png");
}

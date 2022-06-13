#include "draw.C"

Double_t Resolution_Full(Double_t *x_val, Double_t *par)
{
    Double_t y;
    Double_t chi = x_val[0];
    Double_t arg = chi*chi/4.0;
    Double_t norm = TMath::Sqrt(TMath::Pi()/2.0)/2.0;

    y = norm*chi*TMath::Exp(-1.0*arg)*(TMath::BesselI0(arg)+TMath::BesselI1(arg));

    return y;
}
Double_t Resolution_12(Double_t *x_val, Double_t *par)
{
    Double_t y;
    Double_t chi = x_val[0];
    Double_t arg = chi*chi/4.0;
    Double_t norm = TMath::Sqrt(TMath::Pi())/(2.0*TMath::Sqrt(2.0));
    Double_t besselOneHalf = TMath::Sqrt(2.0*arg/TMath::Pi()) * TMath::SinH(arg)/arg;
    Double_t besselThreeHalf = TMath::Sqrt(2.0*arg/TMath::Pi()) * (TMath::CosH(arg)/arg - TMath::SinH(arg)/(arg*arg) );

    y = norm * chi * TMath::Exp(-1.0*arg) * (besselOneHalf + besselThreeHalf);

    return y;
}

void resSubEP()
{
    gROOT->Reset();
    gStyle->SetOptFit(11101);
    //setstyle();
    //TFile *infile=TFile::Open("/star/u/dchen/ana/3gev_2018/shift/test_shift.root");
    TFile *infile=TFile::Open("/star/data01/pwg/dchen/Ana/19p6GeV/shift/19gev_shifted.root");

    TProfile *p_r1_epd_ABCD_sub0_1_input = (TProfile*)infile->Get("p_r1_epd_ABCD_sub0_1");
    TProfile *p_r1_epd_ABCD_wt_sub0_1_input = (TProfile*)infile->Get("p_r1_epd_ABCD_wt_sub0_1");
    TProfile *p_r2_tpc_AB_sub0_1_input = (TProfile*)infile->Get("p_r2_tpc_AB_sub0_1");
    
    const int n= 9;
    double cent[n]={2.5,7.5,15,25,35,45,55,65,75};
    //double cent[n]={2.5,7.5,12.5,17.5,22.5,27.5,32.5,37.5,42.5,47.5,52.5,57.5,62.5,67.5,72.5,77.5};
    //double zero[n]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double zero[n]={0,0,0,0,0,0,0,0,0};

    /*double val_tAtB[n], val_tAeA[n], val_tAeB[n], val_tAeC[n], val_tAeD[n];
    double              val_tBeA[n], val_tBeB[n], val_tBeC[n], val_tBeD[n];
    double val_eAeB[n], val_eAeC[n], val_eAeD[n], val_eBeC[n], val_eBeD[n], val_eCeD[n];
    double val_tAeAB[n], val_tAeCD[n], val_tAeABCD[n];
    double val_tBeAB[n], val_tBeCD[n], val_tBeABCD[n];
    double val_eABeCD[n];
    double val_eAeCD[n], val_eBeCD[n];
    double val_eCeAB[n], val_eDeAB[n];*/
    
    double val_epA_epB[n], val_epA_epC[n], val_epB_epC[n];

    for(int i=0; i<n; i++)
    {
        /*val_tAeA[i] = p_tAeA->GetBinContent(n-i);
        val_tAeB[i] = p_tAeB->GetBinContent(n-i);
        val_tAeC[i] = p_tAeC->GetBinContent(n-i);
        val_tAeD[i] = p_tAeD->GetBinContent(n-i);
        val_tAeAB[i] = p_tAeAB->GetBinContent(n-i);
        val_tAeCD[i] = p_tAeCD->GetBinContent(n-i);
        val_tAeABCD[i] = p_tAeABCD->GetBinContent(n-i);

        val_tBeA[i] = p_tBeA->GetBinContent(n-i);
        val_tBeB[i] = p_tBeB->GetBinContent(n-i);
        val_tBeC[i] = p_tBeC->GetBinContent(n-i);
        val_tBeD[i] = p_tBeD->GetBinContent(n-i);
        val_tBeAB[i] = p_tBeAB->GetBinContent(n-i);
        val_tBeCD[i] = p_tBeCD->GetBinContent(n-i);
        val_tBeABCD[i] = p_tBeABCD->GetBinContent(n-i);

        val_eAeB[i] = p_eAeB->GetBinContent(n-i);
        val_eAeC[i] = p_eAeB->GetBinContent(n-i);
        val_eAeD[i] = p_eAeB->GetBinContent(n-i);
        val_eBeC[i] = p_eBeC->GetBinContent(n-i);
        val_eBeD[i] = p_eBeD->GetBinContent(n-i);
        val_eCeD[i] = p_eCeD->GetBinContent(n-i);
        val_eABeCD[i] = p_eABeCD->GetBinContent(n-i);
        val_eAeCD[i] = p_eAeCD->GetBinContent(n-i);
        val_eBeCD[i] = p_eBeCD->GetBinContent(n-i);
        val_eCeAB[i] = p_eCeAB->GetBinContent(n-i);
        val_eDeAB[i] = p_eDeAB->GetBinContent(n-i);

        val_tAtB[i] = p_tAtB->GetBinContent(n-i);*/

        val_epA_epB[i] = p_r1_epd_ABCD_sub0_1_input->GetBinContent(n-i);
        val_epA_epC[i] = p_r1_epd_ABCD_wt_sub0_1_input->GetBinContent(n-i);
        val_epB_epC[i] = p_r2_tpc_AB_sub0_1_input->GetBinContent(n-i);

        //cout << "res = " << sqrt(val_tAtB[i]) << endl;
    }

    double res1_epd_ABCD_sub0_1[n];
    double res1_epd_ABCD_wt_sub0_1[n];
    double res2_tpc_AB_sub0_1[n];
    
    /*double res1_eABeCDtA[n];
    double res1_eABeCDtB[n];
    double res1_eABtAtB[n];
    double res1_eCDeABtA[n];
    double res1_eCDeABtB[n];
    double res1_eCDtAtB[n];
    double res1_eABCDtAtB[n];
    double res1_eABeCtB[n];
    double res1_eABeDtB[n];
    double res1_eCDeAtB[n];
    double res1_eCDeBtB[n];(*/

    for(int i=0; i<n; i++)
    {
        res1_epd_ABCD_sub0_1[i] = sqrt( val_epA_epB[i]);// * val_epA_epC[i] / val_epB_epC[i] );
        res1_epd_ABCD_wt_sub0_1[i] = sqrt( val_epA_epC[i]);
        res2_tpc_AB_sub0_1[i] = sqrt( val_epB_epC[i]);

        /*res1_eABeCDtA[i] = sqrt( val_eABeCD[i] * val_tAeAB[i] / val_tAeCD[i] );
        res1_eABeCDtB[i] = sqrt( val_eABeCD[i] * val_tBeAB[i] / val_tBeCD[i] );
        res1_eABtAtB[i] = sqrt( val_tAeAB[i] * val_tBeAB[i] / val_tAtB[i] );
        res1_eABeCtB[i] = sqrt( val_eCeAB[i] * val_tBeAB[i] / val_tBeC[i] );
        res1_eABeDtB[i] = sqrt( val_eDeAB[i] * val_tBeAB[i] / val_tBeD[i] );

        //cout << "res = " << res1_eABeCDtB[i] << endl;

        res1_eCDeABtA[i] = sqrt( val_eABeCD[i] * val_tAeCD[i] / val_tAeAB[i] );
        res1_eCDeABtB[i] = sqrt( val_eABeCD[i] * val_tBeCD[i] / val_tBeAB[i] );
        res1_eCDtAtB[i] = sqrt( val_tAeCD[i] * val_tBeCD[i] / val_tAtB[i] );
        res1_eCDeAtB[i] = sqrt( val_eAeCD[i] * val_tBeCD[i] / val_tBeA[i] );
        res1_eCDeBtB[i] = sqrt( val_eBeCD[i] * val_tBeCD[i] / val_tBeB[i] );

        //cout << "res1 = " << res1_eCDeABtB[i] << endl;
        //cout << "res1 = " << res1_eABeCtB[i] << endl;
        //cout << "res1 = " << res1_eABeDtB[i] << endl;
        //cout << "res1 = " << res1_eABeCtB[i] << endl;*/
        cout << "res1 = " << res1_epd_ABCD_sub0_1[i] << endl;

        //res1_eABCDtAtB[i] = sqrt( val_tAeABCD[i] * val_tBeABCD[i] / val_tAtB[i] );
    }

    cout << "from 70-80% -...- 0-5% " << endl;
    cout << "res1 eC = ";
    //cout << "res1 eC = ";
    //for(int i=8; i>=0; i--){cout << res1_eABeCtB[i] << ", ";}cout << endl;
    for(int i=8; i>=0; i--){cout << res1_epd_ABCD_sub0_1[i] << ", ";}cout << endl;
    //cout << "res1 eD = ";
    //for(int i=8; i>=0; i--){cout << res1_eABeDtB[i] << ", ";}cout << endl;
    //cout << "res1 eCD = ";
    //for(int i=8; i>=0; i--){cout << res1_eABeCDtB[i] << ", ";}cout << endl;

    /*Double_t res12_eCDeABtB[n];
    Double_t res12_eABeCtB[n];
    Double_t res12_eABeDtB[n];
    Double_t res12_eABeCDtB[n];

     for(int i=0; i<n; i++)
     {
         TF1 *f_res1 = new TF1("f_res1",Resolution_Full,0,10,0);
         Float_t chi1_eCDeABtB = f_res1->GetX(res1_eCDeABtB[i]);
         Float_t chi1_eABeCtB = f_res1->GetX(res1_eABeCtB[i]);
         Float_t chi1_eABeDtB = f_res1->GetX(res1_eABeDtB[i]);
         Float_t chi1_eABeCDtB = f_res1->GetX(res1_eABeCDtB[i]);
         TF1 *f_res12 = new TF1("f_res12",Resolution_12,0,10,0);
         res12_eCDeABtB[i] = f_res12->Eval(chi1_eCDeABtB);
         res12_eABeCtB[i] = f_res12->Eval(chi1_eABeCtB);
         res12_eABeDtB[i] = f_res12->Eval(chi1_eABeDtB);
         res12_eABeCDtB[i] = f_res12->Eval(chi1_eABeCDtB);

         //cout << "res12 = " << res12_eCDeABtB[i] << endl;
         //cout << "res12 = " << res12_eABeCtB[i] << endl;
         //cout << "res12 = " << res12_eABeDtB[i] << endl;
         //cout << "res12 = " << res12_eABeCtB[i] << endl;
     }
    cout << "res12 eC = ";
    for(int i=8; i>=0; i--){cout << res12_eABeCtB[i] << ", ";}cout << endl;*/
    //cout << "res1 eD = ";
    //for(int i=8; i>=0; i--){cout << res12_eABeDtB[i] << ", ";}cout << endl;
    //cout << "res1 eCD = ";
    //for(int i=8; i>=0; i--){cout << res12_eABeCDtB[i] << ", ";}cout << endl;





    TGraphErrors *gr_epAepBepC =new TGraphErrors(n, cent, res1_epd_ABCD_sub0_1, zero, zero);
    TGraphErrors *gr_epBepCepA =new TGraphErrors(n, cent, res1_epd_ABCD_wt_sub0_1, zero, zero);
    TGraphErrors *gr_epCepAepB =new TGraphErrors(n, cent, res2_tpc_AB_sub0_1, zero, zero);
    /*TGraphErrors *gr_eABeCDtA =new TGraphErrors(n, cent, res1_eABeCDtA, zero, zero);
    TGraphErrors *gr_eABeCDtB =new TGraphErrors(n, cent, res1_eABeCDtB, zero, zero);
    TGraphErrors *gr_eABtAtB =new TGraphErrors(n, cent, res1_eABtAtB, zero, zero);

    TGraphErrors *gr_eABeCtB = new TGraphErrors(n, cent, res1_eABeCtB, zero, zero);
    TGraphErrors *gr_eABeDtB = new TGraphErrors(n, cent, res1_eABeDtB, zero, zero);
    TGraphErrors *gr_eCDeAtB = new TGraphErrors(n, cent, res1_eCDeAtB, zero, zero);
    TGraphErrors *gr_eCDeBtB = new TGraphErrors(n, cent, res1_eCDeBtB, zero, zero);

    TGraphErrors *gr_eCDeABtA =new TGraphErrors(n, cent, res1_eCDeABtA, zero, zero);
    TGraphErrors *gr_eCDeABtB =new TGraphErrors(n, cent, res1_eCDeABtB, zero, zero);
    TGraphErrors *gr_eCDtAtB =new TGraphErrors(n, cent, res1_eCDtAtB, zero, zero);
    TGraphErrors *gr_eABCDtAtB =new TGraphErrors(n, cent, res1_eABCDtAtB, zero, zero);*/

    TCanvas *c1=new TCanvas("c1", "", 600, 600);
    c1->cd()->Draw();
    c1->SetGridx(1);
    c1->SetGridy(2);
    TH2D *h_1=new TH2D("h_1","", 16, -5, 85, 16, -0.05, 0.90);
    h_1->SetStats(0);
    h_1->GetXaxis()->SetTitle("% Centrality");
    h_1->GetYaxis()->SetTitle("Resolution");
    h_1->Draw();

    //EPD-AB
    gr_epBepCepA->SetMarkerStyle(24);
    gr_epBepCepA->SetMarkerColor(2);
    gr_epBepCepA->SetLineColor(2);
    gr_epBepCepA->SetMarkerSize(1.8);
    gr_epBepCepA->Draw("PE");

    gr_epCepAepB->SetMarkerStyle(25);
    gr_epCepAepB->SetMarkerColor(4);
    gr_epCepAepB->SetLineColor(4);
    gr_epCepAepB->SetMarkerSize(1.8);
    gr_epCepAepB->Draw("PE");

    /*gr_eABtAtB->SetMarkerStyle(29);
    gr_eABtAtB->SetMarkerColor(4);
    gr_eABtAtB->SetLineColor(4);
    gr_eABtAtB->SetMarkerSize(1.8);
    //gr_eABtAtB->Draw("PE");

    gr_eCDtAtB->SetMarkerStyle(30);
    gr_eCDtAtB->SetMarkerColor(2);
    gr_eCDtAtB->SetLineColor(2);
    gr_eCDtAtB->SetMarkerSize(1.8);
    //gr_eCDtAtB->Draw("PE");

    gr_eABCDtAtB->SetMarkerStyle(26);
    gr_eABCDtAtB->SetMarkerColor(9);
    gr_eABCDtAtB->SetLineColor(9);
    gr_eABCDtAtB->SetMarkerSize(1.8);
    //gr_eABCDtAtB->Draw("PE");

    gr_eCDeABtB->SetMarkerStyle(20);
    gr_eCDeABtB->SetMarkerColor(4);
    gr_eCDeABtB->SetLineColor(4);
    gr_eCDeABtB->SetMarkerSize(1.8);
    //gr_eCDeABtB->Draw("PE");

    gr_eCDeAtB->SetMarkerStyle(24);
    gr_eCDeAtB->SetMarkerColor(4);
    gr_eCDeAtB->SetLineColor(4);
    gr_eCDeAtB->SetMarkerSize(1.8);
    //gr_eCDeAtB->Draw("PE");

    gr_eCDeBtB->SetMarkerStyle(25);
    gr_eCDeBtB->SetMarkerColor(4);
    gr_eCDeBtB->SetLineColor(4);
    gr_eCDeBtB->SetMarkerSize(1.8);
    //gr_eCDeBtB->Draw("PE");*/

    gr_epAepBepC->SetMarkerStyle(22);
    gr_epAepBepC->SetMarkerColor(2);
    gr_epAepBepC->SetLineColor(2);
    gr_epAepBepC->SetMarkerSize(1.8);
    gr_epAepBepC->Draw("PE");

    /*gr_eABeCDtB->SetMarkerStyle(22);
    gr_eABeCDtB->SetMarkerColor(2);
    gr_eABeCDtB->SetLineColor(2);
    gr_eABeCDtB->SetMarkerSize(1.8);
    gr_eABeCDtB->Draw("PE");

    gr_eABeDtB->SetMarkerStyle(26);
    gr_eABeDtB->SetMarkerColor(2);
    gr_eABeDtB->SetLineColor(2);
    gr_eABeDtB->SetMarkerSize(1.8);
    gr_eABeDtB->Draw("PE");

    gr_eABeCtB->SetMarkerStyle(32);
    gr_eABeCtB->SetMarkerColor(2);
    gr_eABeCtB->SetLineColor(2);
    gr_eABeCtB->SetMarkerSize(1.8);
    gr_eABeCtB->Draw("PE");*/

    TLine *l_zero=new TLine(-5,0.0,85,0.0);
    l_zero->SetLineStyle(2);
    l_zero->SetLineColor(1);
    l_zero->Draw();

    TLegend *leg=new TLegend(0.5,0.6,0.899,0.899);
    leg->SetFillColor(0);
    leg->AddEntry(gr_epAepBepC, "EPD sub r1","p");
    leg->AddEntry(gr_epBepCepA, "EPD sub r1 w/ eta weighting","p");
    leg->AddEntry(gr_epCepAepB, "TPC sub r2","p");
    //leg->AddEntry(gr_eABeCDtA, "EPD-AB vs. EPD-CD and TPC-A","p");
    //leg->AddEntry(gr_eCDeABtA, "EPD-CD vs. EPD-AB and TPC-A","p");
    //leg->AddEntry(gr_eABeCDtB, "EPD-AB vs. EPD-CD and TPC-B","p");
    //leg->AddEntry(gr_eABtAtB, "EPD-AB vs. TPC-A and TPC-B","p");
    //leg->AddEntry(gr_eCDtAtB, "EPD-CD vs. TPC-A and TPC-B","p");
    //leg->AddEntry(gr_eABCDtAtB, "EPD-ABCD vs. TPC-A and TPC-B","p");
    //leg->AddEntry(gr_eABeCtB, "EPD-AB vs. EPD-C and TPC-B", "p");
    //leg->AddEntry(gr_eABeDtB, "EPD-AB vs. EPD-D and TPC-B", "p");
    //leg->AddEntry(gr_eCDeABtB, "EPD-CD vs. EPD-AB and TPC-B","p");
    //leg->AddEntry(gr_eCDeAtB, "EPD-CD vs. EPD-A and TPC-B","p");
    //leg->AddEntry(gr_eCDeBtB, "EPD-CD vs. EPD-B and TPC-B","p");
    leg->SetLineColor(0);
    leg->SetFillColor(0);
    leg->Draw();

    c1->SaveAs("resSubEP.pdf");


}

#include "TH4D.h"

ClassImp(TH4D);

TH4D::TH4D(){
  fXaxis.SetName("xaxis");
}
TH4D::~TH4D(){
  for(auto hist:hists)
    if(hist) delete hist;
  hists.clear();
}
TH4D::TH4D(const char *name,const char *title,Int_t nbinsx,Double_t xlow,Double_t xup
	   ,Int_t nbinsy,Double_t ylow,Double_t yup
	   ,Int_t nbinsz,Double_t zlow,Double_t zup
	   ,Int_t nbinsu,Double_t ulow,Double_t uup){
  SetNameTitle(name,title);
  fXaxis.SetName("xaxis");
  if (nbinsx <= 0) {
    Warning("TH4D","nbinsx is <=0 - set to nbinsx = 1");
    nbinsx = 1;
  }
  fXaxis.Set(nbinsx,xlow,xup);
  for(int i=0;i<nbinsx+2;i++){
    hists.push_back(new TH3D(Form("%d",i),Form("%d",i),nbinsy,ylow,yup,nbinsz,zlow,zup,nbinsu,ulow,uup));
    hists.back()->Sumw2();
    hists.back()->SetDirectory(0);
  }
}
TH4D::TH4D(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins
	   ,Int_t nbinsy,const Double_t *ybins
	   ,Int_t nbinsz,const Double_t *zbins
	   ,Int_t nbinsu,const Double_t *ubins){
  SetNameTitle(name,title);
  fXaxis.SetName("xaxis");
  if (nbinsx <= 0) {
    Warning("TH4D","nbinsx is <=0 - set to nbinsx = 1");
    nbinsx = 1;
  }
  fXaxis.Set(nbinsx,xbins);
  for(int i=0;i<nbinsx+2;i++){
    hists.push_back(new TH3D(Form("%d",i),Form("%d",i),nbinsy,ybins,nbinsz,zbins,nbinsu,ubins));
    hists.back()->Sumw2();
    hists.back()->SetDirectory(0);
  }
}
TH4D::TH4D(TString path){
  if(path.EndsWith(".root")){
    ReadFile(path);
    return;
  }
  Error("TH4D","path should end with .root");
}

Int_t TH4D::Fill(Double_t x, Double_t y, Double_t z, Double_t u, Double_t w){
  return hists.at(fXaxis.FindBin(x))->Fill(y,z,u,w);
}
Int_t TH4D::Fill(Double_t x, Double_t y, Double_t z, Double_t u){
  return hists.at(fXaxis.FindBin(x))->Fill(y,z,u);
}
    
Int_t TH4D::WriteFile(TString path){
  TFile* f=NULL;
  if(path.EndsWith(".root")){
    f=new TFile(path,"recreate");
  }else{
    if(system("mkdir -p "+path)==0){
      f=new TFile(path+"/"+GetName()+".root","recreate");
    }else{
      Error("TH4D","cannot make directory");
      return false;
    }
  }
  fXaxis.Write();
  for(auto& hist:hists) hist->Write();
  delete f;
  return true;
}

Int_t TH4D::ReadFile(TString path){
  TFile f(path);
  TAxis* axis=(TAxis*)f.Get("xaxis");
  if(axis){
    axis->Copy(fXaxis);
    int nbinsx=fXaxis.GetNbins();
    for(int i=0;i<nbinsx+2;i++){
      TH3D* hist=(TH3D*)f.Get(Form("%d",i));
      hist->SetDirectory(0);
      hists.push_back(hist);
    }
  }else{
    Error("TH4D","no xaxis in file");
    return false;
  }
  return true;
}

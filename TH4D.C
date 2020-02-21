#include "TH4D.h"

ClassImp(TH4D);

TH4D::TH4D(){
  fXaxis.SetName("uaxis");
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
  fXaxis.SetName("uaxis");
  if (nbinsu <= 0) {
    Warning("TH4D","nbinsu is <=0 - set to nbinsu = 1");
    nbinsu = 1;
  }
  fXaxis.Set(nbinsu,ulow,uup);
  for(int i=0;i<nbinsu+2;i++){
    TH3D* hist=new TH3D(TString::Format("%s_ubin%d",GetName(),i),Form("%d",i),nbinsx,xlow,xup,nbinsy,ylow,yup,nbinsz,zlow,zup);
    hist->SetDirectory(GetDirectory());
    hists.push_back(hist);
    //hists.back()->Sumw2();
  }
}
TH4D::TH4D(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins
	   ,Int_t nbinsy,const Double_t *ybins
	   ,Int_t nbinsz,const Double_t *zbins
	   ,Int_t nbinsu,const Double_t *ubins){
  SetNameTitle(name,title);
  fXaxis.SetName("uaxis");
  if (nbinsu <= 0) {
    Warning("TH4D","nbinsu is <=0 - set to nbinsu = 1");
    nbinsu = 1;
  }
  fXaxis.Set(nbinsu,ubins);
  for(int i=0;i<nbinsu+2;i++){
    TH3D* hist=new TH3D(TString::Format("%s_ubin%d",GetName(),i),Form("%d",i),nbinsx,xbins,nbinsy,ybins,nbinsz,zbins);
    hist->SetDirectory(GetDirectory());
    hists.push_back(hist);
    //hists.back()->Sumw2();
  }
}

void TH4D::SetDirectory(TDirectory *dir){
  TH1::SetDirectory(dir);
  for(auto hist:hists) hist->SetDirectory(dir);
}

Int_t TH4D::Fill(Double_t x, Double_t y, Double_t z, Double_t u, Double_t w){
  return hists.at(GetUaxis()->FindBin(u))->Fill(x,y,z,w);
}
Int_t TH4D::Fill(Double_t x, Double_t y, Double_t z, Double_t u){
  return hists.at(GetUaxis()->FindBin(u))->Fill(x,y,z);
}

Double_t TH4D::IntegralAndError(Int_t binx1, Int_t binx2, Int_t biny1, Int_t biny2, Int_t binz1, Int_t binz2, Int_t binu1, Int_t binu2, Double_t & error, Option_t *option) const {
  Double_t integral=0;
  Double_t error2sum=0;
  for(int i=binu1;i<=binu2;i++){
    Double_t this_error=0;
    integral+=hists[i]->IntegralAndError(binx1,binx2,biny1,biny2,binz1,binz2,this_error,option);
    error2sum=this_error*this_error;
  }
  error=sqrt(error2sum);
  return integral;
}
Double_t TH4D::Integral(Int_t binx1, Int_t binx2, Int_t biny1, Int_t biny2, Int_t binz1, Int_t binz2, Int_t binu1, Int_t binu2, Option_t *option) const {
  Double_t dummy=0;
  return IntegralAndError(binx1,binx2,biny1,biny2,binz1,binz2,binu1,binu2,dummy,option);
}
Double_t TH4D::Integral(Option_t *option) const {
  return Integral(GetXaxis()->GetFirst(),GetXaxis()->GetLast(),
		  GetYaxis()->GetFirst(),GetYaxis()->GetLast(),
		  GetZaxis()->GetFirst(),GetZaxis()->GetLast(),
		  GetUaxis()->GetFirst(),GetUaxis()->GetLast(),option);
}

TH1D* TH4D::ProjectionX(const char *name, Int_t iymin, Int_t iymax, Int_t izmin, Int_t izmax, Int_t iumin, Int_t iumax, Option_t *option) const {
  TH1D* hist=NULL;
  if(iumax==-1) iumax=GetUaxis()->GetNbins()+1;
  for(int i=iumin;i<=iumax;i++){
    TH1D* this_hist=hists.at(i)->ProjectionX(name,iymin,iymax,izmin,izmax,option);
    if(hist){
      hist->Add(this_hist);
      delete this_hist;
    }
    else hist=this_hist;
  }
  return hist;
}
TH1D* TH4D::ProjectionY(const char *name, Int_t ixmin, Int_t ixmax, Int_t izmin, Int_t izmax, Int_t iumin, Int_t iumax, Option_t *option) const {
  TH1D* hist=NULL;
  if(iumax==-1) iumax=GetUaxis()->GetNbins()+1;
  for(int i=iumin;i<=iumax;i++){
    TH1D* this_hist=hists.at(i)->ProjectionY(name,ixmin,ixmax,izmin,izmax,option);
    if(hist){
      hist->Add(this_hist);
      delete this_hist;
    }
    else hist=this_hist;
  }
  return hist;
}
TH1D* TH4D::ProjectionZ(const char *name, Int_t ixmin, Int_t ixmax, Int_t iymin, Int_t iymax, Int_t iumin, Int_t iumax, Option_t *option) const {
  TH1D* hist=NULL;
  if(iumax==-1) iumax=GetUaxis()->GetNbins()+1;
  for(int i=iumin;i<=iumax;i++){
    TH1D* this_hist=hists.at(i)->ProjectionZ(name,ixmin,ixmax,iymin,iymax,option);
    if(hist){
      hist->Add(this_hist);
      delete this_hist;
    }
    else hist=this_hist;
  }
  return hist;
}
TH1D* TH4D::ProjectionU(const char *name, Int_t ixmin, Int_t ixmax, Int_t iymin, Int_t iymax, Int_t izmin, Int_t izmax, Option_t *option) const {
  TString hname = name;
  if (hname == "_px") hname = TString::Format("%s%s", GetName(), name);
  TString title =  TString::Format("%s ( Projection X )",GetTitle());
  TH1D* hist=new TH1D(hname,title,GetUaxis()->GetNbins(),GetUaxis()->GetXbins()->GetArray());
  for(int i=0;i<GetUaxis()->GetNbins()+2;i++){
    double content=0,error=0;
    content=hists.at(i)->IntegralAndError(ixmin,ixmax,iymin,iymax,izmin,izmax,error,option);
    hist->SetBinContent(i,content);
    hist->SetBinError(i,error);
  }
  return hist;
}
bool TH4D::CheckConsistency() const {
  int nbinsu=GetUaxis()->GetNbins();
  bool consistent=true;
  for(int i=1;i<nbinsu+2;i++){
    consistent&=TH1::CheckConsistency((TH1*)hists.at(0),(TH1*)hists.at(i));
  }
  if(!consistent){
    Error("TH4D","not consistent");
    return false;
  }
  return true;
}
bool TH4D::CheckConsistency(const TH4D* h1,const TH4D* h2){
  if(!h1->CheckConsistency()) return false;
  if(!h2->CheckConsistency()) return false;
  if(!CheckEqualAxes(h1->GetUaxis(),h2->GetUaxis())) return false; 
  return true;
}

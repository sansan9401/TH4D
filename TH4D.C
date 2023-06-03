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
    hist->SetDirectory(NULL);
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
    hist->SetDirectory(NULL);
    hists.push_back(hist);
    //hists.back()->Sumw2();
  }
}

Bool_t TH4D::Add(const TH1 *h1, Double_t c1){
  if (!h1) {
    Error("Add","Attempt to add a non-existing histogram");
    return kFALSE;
  }
  if(!CheckConsistency((TH4D*)h1)){
    Error("Add","Can't add. Not consistent");
    return kFALSE;
  }
  int nubins=GetUaxis()->GetNbins();
  for(int i=0;i<nubins+2;i++){
    hists.at(i)->Add((TH1*)((TH4D*)h1)->hists.at(i),c1);
  }
  return true;
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
bool TH4D::CheckConsistency(const TH4D* h1,const TH4D* h2) {
  if(!h1->CheckConsistency()) return false;
  if(!h2->CheckConsistency()) return false;
  if(!CheckEqualAxes(h1->GetUaxis(),h2->GetUaxis())) return false; 
  return true;
}

TObject* TH4D::Clone(const char* newname) const {
  TH4D* h=(TH4D*)TH1::Clone(newname);
  int nubins=GetUaxis()->GetNbins();
  for(int i=0;i<nubins+2;i++){
    h->hists.push_back((TH3D*)hists.at(i)->Clone(TString::Format("%s_ubin%d",h->GetName(),i)));
    h->hists.back()->SetDirectory(NULL);
  }
  return (TObject*)h;
}

void TH4D::Draw(Option_t *option){
  TString opt=option; opt.ToLower();
  if(opt.Contains("px")) ProjectionX()->Draw(option);
  else if(opt.Contains("py")) ProjectionY()->Draw(option);
  else if(opt.Contains("pz")) ProjectionZ()->Draw(option);
  else ProjectionU()->Draw(option);
}

Bool_t TH4D::Divide(const TH1 *h1){
  if (!h1) {
    Error("Divide","Attempt to divide a non-existing histogram");
    return kFALSE;
  }
  if(!CheckConsistency((TH4D*)h1)){
    Error("Divide","Can't divide. Not consistent");
    return kFALSE;
  }
  int nubins=GetUaxis()->GetNbins();
  for(int i=0;i<nubins+2;i++){
    hists.at(i)->Divide((TH1*)((TH4D*)h1)->hists.at(i));
  }
  return true;  
}

Bool_t TH4D::Divide(const TH1 *h1, const TH1 *h2, Double_t c1, Double_t c2, Option_t *option){
  if (!h1 || !h2) {
    Error("Divide","Attempt to divide a non-existing histogram");
    return kFALSE;
  }
  if(!CheckConsistency((TH4D*)h1)||!CheckConsistency((TH4D*)h2)){
    Error("Divide","Can't divide. Not consistent");
    return kFALSE;
  }
  int nubins=GetUaxis()->GetNbins();
  for(int i=0;i<nubins+2;i++){
    hists.at(i)->Divide((TH1*)((TH4D*)h1)->hists.at(i),(TH1*)((TH4D*)h2)->hists.at(i),c1,c2,option);
  }
  return true;  
}

Int_t TH4D::GetBin(Int_t binx, Int_t biny, Int_t binz, Int_t binu) const {
  Int_t ofx = GetXaxis()->GetNbins() + 1;
  if (binx < 0) binx = 0;
  if (binx > ofx) binx = ofx;
 
  Int_t ofy = GetYaxis()->GetNbins() + 1;
  if (biny < 0) biny = 0;
  if (biny > ofy) biny = ofy;

  Int_t ofz = GetZaxis()->GetNbins() + 1;
  if (binz < 0) binz = 0;
  if (binz > ofz) binz = ofz;

  Int_t ofu = GetUaxis()->GetNbins() + 1;
  if (binu < 0) binu = 0;
  if (binu > ofu) binu = ofu;

  return binx + (GetXaxis()->GetNbins() + 2) * (biny + (GetYaxis()->GetNbins() + 2) * (binz + (GetZaxis()->GetNbins() + 2) * binu ));
}
void TH4D::GetBinXYZU(Int_t binglobal, Int_t &binx, Int_t &biny, Int_t &binz, Int_t &binu) const
{
  Int_t nx  = GetXaxis()->GetNbins()+2;
  Int_t ny  = GetYaxis()->GetNbins()+2;
  Int_t nz  = GetZaxis()->GetNbins()+2;

  binx = binglobal%nx;
  biny = ((binglobal-binx)/nx)%ny;
  binz = (((binglobal-binx)/nx -biny)/ny)%nz;
  binu = (((binglobal-binx)/nx -biny)/ny-binz)/nz;
}

Double_t TH4D::GetBinContent(Int_t bin) const {
  Int_t ix,iy,iz,iu;
  GetBinXYZU(bin,ix,iy,iz,iu);
  return hists.at(iu)->GetBinContent(ix,iy,iz);
}
Double_t TH4D::GetBinError(Int_t bin) const {
  Int_t ix,iy,iz,iu;
  GetBinXYZU(bin,ix,iy,iz,iu);
  return hists.at(iu)->GetBinError(ix,iy,iz);
}  
Int_t TH4D::GetNcells() const {
  if(hists.size()&&hists.at(0)) return hists.size()*hists.at(0)->GetNcells();
  else return 0;
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
Double_t TH4D::IntegralAndError(Int_t binx1, Int_t binx2, Int_t biny1, Int_t biny2, Int_t binz1, Int_t binz2, Int_t binu1, Int_t binu2, Double_t & error, Option_t *option) const {
  Double_t integral=0;
  Double_t error2sum=0;

  Int_t nu = GetNbinsU() + 2;
  if (binu1 < 0) binu1 = 0;
  if (binu2 >= nu || binu2 < binu1) binu2 = nu - 1;

  for(int i=binu1;i<=binu2;i++){
    Double_t this_error=0;
    integral+=hists[i]->IntegralAndError(binx1,binx2,biny1,biny2,binz1,binz2,this_error,option);
    error2sum=this_error*this_error;
  }
  error=sqrt(error2sum);
  return integral;
}

Long64_t TH4D::Merge(TCollection *list){
  TIter iter(list);
  while(TObject *obj=iter()){
    if(obj) Add((TH1*)obj);
  }    
  return true;
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
    else{
      hist=this_hist;
      TString hname = name;
      if (hname == "_px") hname = TString::Format("%s%s", GetName(), name);
      hist->SetName(hname);
      hist->SetTitle(TString::Format("%s ( Projection X )",GetTitle()));
    }
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
    else{
      hist=this_hist;
      TString hname = name;
      if (hname == "_py") hname = TString::Format("%s%s", GetName(), name);
      hist->SetName(hname);
      hist->SetTitle(TString::Format("%s ( Projection Y )",GetTitle()));
    }
  }
  return hist;
}
TH1D* TH4D::ProjectionZ(const char *name, Int_t ixmin, Int_t ixmax, Int_t iymin, Int_t iymax, Int_t iumin, Int_t iumax, Option_t *option) const {
  TH1D* hist=NULL;
  if(iumax==-1) iumax=GetUaxis()->GetNbins()+1;
  for(int i=iumin;i<=iumax;i++){
    TH1D* this_hist=hists.at(i)->ProjectionZ("_pz",ixmin,ixmax,iymin,iymax,option);
    if(hist){
      hist->Add(this_hist);
      delete this_hist;
    }
    else{
      hist=this_hist;
      TString hname = name;
      if (hname == "_pz") hname = TString::Format("%s%s", GetName(), name);
      hist->SetName(hname);
      hist->SetTitle(TString::Format("%s ( Projection Z )",GetTitle()));
    }
  }
  return hist;
}
TH1D* TH4D::ProjectionU(const char *name, Int_t ixmin, Int_t ixmax, Int_t iymin, Int_t iymax, Int_t izmin, Int_t izmax, Option_t *option) const {
  TString hname = name;
  if (hname == "_pu") hname = TString::Format("%s%s", GetName(), name);
  TString title =  TString::Format("%s ( Projection U )",GetTitle());
  TH1D* hist=NULL;
  if(GetUaxis()->GetXbins()->GetArray()) hist=new TH1D(hname,title,GetUaxis()->GetNbins(),GetUaxis()->GetXbins()->GetArray());
  else hist=new TH1D(hname,title,GetUaxis()->GetNbins(),GetUaxis()->GetXmin(),GetUaxis()->GetXmax());
  for(int i=0;i<GetUaxis()->GetNbins()+2;i++){
    double content=0,error=0;
    content=hists.at(i)->IntegralAndError(ixmin,ixmax,iymin,iymax,izmin,izmax,error,option);
    hist->SetBinContent(i,content);
    hist->SetBinError(i,error);
  }
  return hist;
}

void TH4D::Scale(Double_t c1, Option_t *option){
  int nubins=GetUaxis()->GetNbins();
  for(int i=0;i<nubins+2;i++){
    hists.at(i)->Scale(c1,option);
  }
}

void TH4D::SetBinContent(Int_t bin, Double_t content){
  Int_t ix,iy,iz,iu;
  GetBinXYZU(bin,ix,iy,iz,iu);
  return hists.at(iu)->SetBinContent(ix,iy,iz,content);
}  

void TH4D::SetBinError(Int_t bin, Double_t error){
  Int_t ix,iy,iz,iu;
  GetBinXYZU(bin,ix,iy,iz,iu);
  return hists.at(iu)->SetBinError(ix,iy,iz,error);
}  

void TH4D::SetDirectory(TDirectory *dir){
  TH1::SetDirectory(dir);
  for(auto hist:hists) hist->SetDirectory(NULL);
}

void TH4D::Sumw2(Bool_t flag){
  TH1::Sumw2(flag);
  for(auto hist:hists) hist->TH1::Sumw2(flag);
}

void TH4D::Reset(Option_t *option){
  TH1::Reset(option);
  int nubins=GetUaxis()->GetNbins();
  for(int i=0;i<nubins+2;i++){
    hists.at(i)->Reset(option);
  } 
}

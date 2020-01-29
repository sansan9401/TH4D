#ifndef HS_TH4D
#define HS_TH4D
#include "TH3.h"
#include "TFile.h"

class TH4D : public TNamed{
 public:
  TAxis fXaxis;
  vector<TH3D*> hists;
  
  TH4D();
  ~TH4D();
  TH4D(const char *name,const char *title,Int_t nbinsx,Double_t xlow,Double_t xup
      ,Int_t nbinsy,Double_t ylow,Double_t yup
      ,Int_t nbinsz,Double_t zlow,Double_t zup
      ,Int_t nbinsu,Double_t ulow,Double_t uup);
  TH4D(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins
      ,Int_t nbinsy,const Double_t *ybins
      ,Int_t nbinsz,const Double_t *zbins
      ,Int_t nbinsu,const Double_t *ubins);
  TH4D(TString path);
  
  virtual Int_t Fill(Double_t x, Double_t y, Double_t z, Double_t u);
  virtual Int_t Fill(Double_t x, Double_t y, Double_t z, Double_t u, Double_t w);
  
  virtual Int_t WriteFile(TString path=".");
  virtual Int_t ReadFile(TString path);

  ClassDef(TH4D,1);
};
  
#endif

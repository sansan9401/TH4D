#ifndef __HS_TH4D__
#define __HS_TH4D__
#include <vector>
#include "TH3.h"
#include "TH2D.h"
#include "TFile.h"

class TH4D : public TH1{
 protected:
  virtual Bool_t Add(const TH1 *h, const TH1 *h2, Double_t c1=1, Double_t c2=1){return false;} //NotImplemented
  virtual Bool_t Add(TF1 *h1, Double_t c1=1, Option_t *option=""){return false;} //NotImplemented
  virtual Int_t Fill(Double_t){return -1;} // MayNotUse
  virtual Int_t Fill(Double_t, Double_t){return -1;} // MayNotUse
  virtual Int_t Fill(const char*, Double_t){return -1;} // MayNotUse
  virtual Int_t GetBin(Int_t binx, Int_t biny=0, Int_t binz=0) const {return -1;} //MayNotUse
  virtual Double_t GetBinContent(Int_t, Int_t) const { return -1; } //MayNotUse
  virtual Double_t GetBinContent(Int_t, Int_t, Int_t) const { return -1; } //MayNotUse
  virtual Double_t GetBinError(Int_t, Int_t) const { return -1; } //MayNotUse
  virtual Double_t GetBinError(Int_t, Int_t, Int_t) const { return -1; } //MayNotUse
  virtual Double_t IntegralAndError(Int_t, Int_t, Double_t & err, Option_t *option="") const{return -1;} // MayNotUse
  virtual Double_t Integral(Int_t, Int_t, Option_t*) const{return -1;} // MayNotUse
  virtual void SetBinContent(Int_t, Int_t, Double_t) { return; } //MayNotUse
  virtual void SetBinContent(Int_t, Int_t, Int_t, Double_t) { return; } //MayNotUse
  virtual void SetBinError(Int_t, Int_t, Double_t) { return; } //MayNotUse
  virtual void SetBinError(Int_t, Int_t, Int_t, Double_t) { return; } //MayNotUse

 public:
  std::vector<TH3D*> hists;
  
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


  virtual Bool_t Add(const TH1 *h1, Double_t c1=1);

  virtual bool CheckConsistency() const;
  virtual bool CheckConsistency(const TH4D* h1) const { return CheckConsistency(h1,this); }
  static bool CheckConsistency(const TH4D* h1,const TH4D* h2);

  virtual Int_t Fill(Double_t x, Double_t y, Double_t z, Double_t u, Double_t w) { return hists.at(GetUaxis()->FindBin(u))->Fill(x,y,z,w); }
  virtual Int_t Fill(Double_t x, Double_t y, Double_t z, Double_t u) { return hists.at(GetUaxis()->FindBin(u))->Fill(x,y,z); };

  virtual Int_t GetBin(Int_t binx, Int_t biny, Int_t binz, Int_t binu) const;
  virtual void GetBinXYZU(Int_t binglobal, Int_t &binx, Int_t &biny, Int_t &binz, Int_t &binu) const;
  virtual Double_t GetBinContent(Int_t bin) const;
  virtual Double_t GetBinContent(Int_t binx, Int_t biny, Int_t binz, Int_t binu) const { return GetBinContent( GetBin(binx, biny, binz, binu) ); }
  virtual Double_t GetBinError(Int_t bin) const;
  virtual Double_t GetBinError(Int_t binx, Int_t biny, Int_t binz, Int_t binu) const { return GetBinError(GetBin(binx, biny, binz, binu)); }

  virtual const TAxis* GetXaxis() const {if(hists.size()) return hists[0]->GetXaxis();else return NULL;}
  virtual const TAxis* GetYaxis() const {if(hists.size()) return hists[0]->GetYaxis();else return NULL;}
  virtual const TAxis* GetZaxis() const {if(hists.size()) return hists[0]->GetZaxis();else return NULL;}
  virtual const TAxis* GetUaxis() const {return &fXaxis;}
  
  virtual Double_t Integral(Int_t binx1, Int_t binx2, Int_t biny1, Int_t biny2, Int_t binz1, Int_t binz2, Int_t binu1, Int_t binu2, Option_t *option="") const;
  virtual Double_t Integral(Option_t *option="") const;
  virtual Double_t IntegralAndError(Int_t binx1, Int_t binx2, Int_t biny1, Int_t biny2, Int_t binz1, Int_t binz2, Int_t binu1, Int_t binu2, Double_t & error, Option_t *option="") const;

  virtual Long64_t Merge(TCollection *list);
  virtual TH1D* ProjectionX(const char *name="_px", Int_t iymin=0, Int_t iymax=-1, Int_t izmin=0, Int_t izmax=-1, Int_t iumin=0, Int_t iumax=-1, Option_t *option="") const;
  virtual TH1D* ProjectionY(const char *name="_py", Int_t ixmin=0, Int_t ixmax=-1, Int_t izmin=0, Int_t izmax=-1, Int_t iumin=0, Int_t iumax=-1, Option_t *option="") const;
  virtual TH1D* ProjectionZ(const char *name="_pz", Int_t ixmin=0, Int_t ixmax=-1, Int_t iymin=0, Int_t iymax=-1, Int_t iumin=0, Int_t iumax=-1, Option_t *option="") const;
  virtual TH1D* ProjectionU(const char *name="_pu", Int_t ixmin=0, Int_t ixmax=-1, Int_t iymin=0, Int_t iymax=-1, Int_t izmin=0, Int_t izmax=-1, Option_t *option="") const;
  virtual TH2D* ProjectionXZ(const char *name="_pxz", Int_t iymin=0, Int_t iymax=-1, Int_t iumin=0, Int_t iumax=-1, Option_t *option="") const;

  virtual void SetBinContent(Int_t bin, Double_t content);
  virtual void SetBinContent(Int_t binx, Int_t biny, Int_t binz, Int_t binu, Double_t content) { SetBinContent(GetBin(binx, biny, binz, binu), content); }
  virtual void SetBinError(Int_t bin, Double_t error);
  virtual void SetBinError(Int_t binx, Int_t biny, Int_t binz, Int_t binu, Double_t error) { SetBinError(GetBin(binx, biny, binz, binu), error); }
  

  virtual void SetDirectory(TDirectory *dir);

  virtual void Sumw2(Bool_t flag=kTRUE);

  //virtual Int_t WriteFile(TString path=".");
  //virtual Int_t ReadFile(TString path);

  ClassDef(TH4D,3);
};
  
#endif

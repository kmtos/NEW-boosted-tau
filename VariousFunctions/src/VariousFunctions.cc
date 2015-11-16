#include "BoostedTauAnalysis/VariousFunctions/interface/VariousFunctions.h"
#include <math.h>
#include "TH1F.h"
#include "THistPainter.h"
#include "DataFormats/Math/interface/deltaPhi.h"

//
//Format and Draw 1D Histograms
//
void VariousFunctions::formatAndDrawCanvasAndHist1D(TCanvas& canvas, TH1F* hist, const Int_t grid, const Int_t logY, const Int_t logZ, const Color_t color, const Size_t size, const Style_t style,  const char* xAxisTitle,
						  const Float_t xTitleSize, const Float_t xLabelSize, const Float_t xTitleOffset, const char* yAxisTitle, const Float_t yTitleSize, const Float_t yLabelSize,
						  const Float_t yTitleOffset, const bool errorBars)
{
  canvas.SetFillStyle(0);
  canvas.SetFillColor(0);
  canvas.SetGrid(grid, grid);
  canvas.SetLogy(logY);
  canvas.SetLogz(logZ);
  hist->SetMarkerColor(color);
  hist->SetMarkerSize(size);
  hist->SetMarkerStyle(style);
  hist->SetLineColor(color);
  hist->SetLineWidth(2);
  
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetXaxis()->SetLabelSize(xLabelSize);
  hist->GetXaxis()->SetTitle(xAxisTitle);
  hist->GetXaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetTitleSize(xTitleSize);
  hist->GetXaxis()->SetTitleOffset(xTitleOffset);
  hist->GetXaxis()->SetTitle(xAxisTitle);
  
  hist->GetYaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelSize(yLabelSize);
  hist->GetYaxis()->SetTitle(yAxisTitle);
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetYaxis()->SetTitleSize(yTitleSize);
  hist->GetYaxis()->SetTitleOffset(yTitleOffset);
  hist->GetYaxis()->SetTitle(yAxisTitle);

  canvas.cd();
  if(errorBars)
    hist->Draw("E1");
  else
    hist->Draw();
}//VariousFunctions::formatAndDrawCanvasAndHist1D

//
//Format and Draw 2D Histograms
//
void VariousFunctions::formatAndDrawCanvasAndHist2D(TCanvas& canvas, TH2F* hist, const Int_t grid, const Int_t logY, const Int_t logZ, const Color_t color, const Size_t size, const Style_t style,
 						   const char* xAxisTitle, const Float_t xTitleSize, const Float_t xLabelSize, const Float_t xTitleOffset,
						   const char* yAxisTitle, const Float_t yTitleSize, const Float_t yLabelSize, const Float_t yTitleOffset,
                                                   const char* zAxisTitle, const Float_t zTitleSize, const Float_t zLabelSize, const Float_t zTitleOffset)
{
  canvas.SetFillStyle(0);
  canvas.SetFillColor(0);
  canvas.SetGrid(grid, grid);
  canvas.SetLogy(logY);
  canvas.SetLogz(logZ);
  canvas.cd()->SetLeftMargin(.2);
  canvas.cd()->SetTopMargin(.2);
  canvas.cd()->SetRightMargin(.2);
  canvas.cd()->SetBottomMargin(.2);

  hist->SetMarkerColor(color);
  hist->SetMarkerSize(size);
  hist->SetMarkerStyle(style);
  hist->SetLineColor(color);
  hist->SetLineWidth(1);
  hist->SetFillStyle(0);

  hist->GetXaxis()->SetLabelFont(42);
  hist->GetXaxis()->SetLabelOffset(0.007);
  hist->GetXaxis()->SetLabelSize(xLabelSize);
  hist->GetXaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetTitleSize(xTitleSize);
  hist->GetXaxis()->SetTitleOffset(xTitleOffset);
  hist->GetXaxis()->SetTitle(xAxisTitle);

  hist->GetYaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelOffset(0.007);
  hist->GetYaxis()->SetLabelSize(yLabelSize);
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetYaxis()->SetTitleSize(yTitleSize);
  hist->GetYaxis()->SetTitleOffset(yTitleOffset);
  hist->GetYaxis()->SetTitle(yAxisTitle);

  hist->GetZaxis()->SetLabelFont(42);
  hist->GetZaxis()->SetLabelOffset(0.007);
  hist->GetZaxis()->SetLabelSize(zLabelSize);
  hist->GetZaxis()->SetTitleFont(42);
  hist->GetZaxis()->SetTitleSize(zTitleSize);
  hist->GetZaxis()->SetTitleOffset(zTitleOffset);
  hist->GetZaxis()->SetTitle(zAxisTitle);

  canvas.cd();
  hist->Draw("COLZ");

}//VariousFunctions::formatAndDrawCanvasAndHist1D


//
//This finds and returns the GenParticleRef a daughter with the givne pdgId of a provided mother
//
reco::GenParticleRef  VariousFunctions::findDaughterInDaughters(const reco::GenParticleRef& mother, const double pdgId, const bool abs)
{
  unsigned int iDaughters = 0;
  reco::GenParticleRef  found;
  reco::GenParticleRef childRef;
  if(abs)
  {
    while (iDaughters < mother->numberOfDaughters())
    {
      childRef = mother->daughterRef(iDaughters);
      if (fabs(childRef->pdgId()) == pdgId)
      {
	found = childRef;
        break;
      }//if fabs childRef
      iDaughters++;
    }//for iGenDaughters
  }//ifabs
 else 
  { 
    while (iDaughters < mother->numberOfDaughters())
    { 
      childRef = mother->daughterRef(iDaughters);
      if (childRef->pdgId() == pdgId)
      { 
        found = childRef;
        break;
      }//if fabs childRef
      iDaughters++;
    }//for iGenDaughters
  }//else
  return found;
}//findDaughter



bool VariousFunctions::findIfInDaughters(const reco::GenParticleRef& mother, const double pdgId, const bool abs)
{
  unsigned int iDaughters = 0;
  bool found = false;
  reco::GenParticleRef childRef;
  if(abs)
  {
    while (iDaughters < mother->numberOfDaughters())
    {
      childRef = mother->daughterRef(iDaughters);
      if (fabs(childRef->pdgId()) == pdgId)
      {
        found = true;
        break;
      }//if fabs childRef
      iDaughters++;
    }//for iGenDaughters
  }//ifabs
 else
  {
    while (iDaughters < mother->numberOfDaughters())
    {
      childRef = mother->daughterRef(iDaughters);
      if (childRef->pdgId() == pdgId)
      {
        found = true;
        break;
      }//if fabs childRef
      iDaughters++;
    }//for iGenDaughters
  }//else
  return found;
}//findDaughter

//
//Takes a tau that has had or mu or ele decay products
//Tau Decay Mode: 1= 1prong , 2= 1pr + 1pi_0 , 3= 1pr + 2pi_0 , 4=3pr , 5=other , 6=electronic , 7=muonic
//
int VariousFunctions::tauDecayMode(const reco::GenParticleRef& tauRef)
{
  if(tauRef->numberOfDaughters() == 3)
  {
  reco::GenParticleRef daughter1 = tauRef->daughterRef(0);
  reco::GenParticleRef daughter2 = tauRef->daughterRef(1);
  reco::GenParticleRef daughter3 = tauRef->daughterRef(2);

  //Leptonic checks
  if(fabs(daughter1->pdgId()) == 13 || fabs(daughter2->pdgId()) == 13 || fabs(daughter3->pdgId()) == 13) 
    return 7;
  if(fabs(daughter1->pdgId()) == 11 || fabs(daughter2->pdgId()) == 11 || fabs(daughter3->pdgId()) == 11) 
    return 6;
  }//numof Daughters == 3

  unsigned int iDaughters = 0, nProngs = 0, nPi0=0;
  while(iDaughters < tauRef->numberOfDaughters())
  {
    reco::GenParticleRef interDaughter = tauRef->daughterRef(iDaughters);
    if(interDaughter->pdgId() == 111 || interDaughter->pdgId() == 310 || interDaughter->pdgId() == 130) //checks for pi_0, KS_0, or KL_0
      nPi0++;
    if(fabs(interDaughter->pdgId()) == 211 || fabs(interDaughter->pdgId()) == 321)  //checks for +-pion or +- Kaon
      nProngs++;
    if(fabs(interDaughter->pdgId()) == 213 || fabs(interDaughter->pdgId()) == 223 || fabs(interDaughter->pdgId()) == 20213)  //checks for +-rho and omega_0
    {
      unsigned int iDaughters2 = 0;
      while(iDaughters2 < interDaughter->numberOfDaughters())
      {
	reco::GenParticleRef iGrandDaughter = interDaughter->daughterRef(iDaughters2);
        if(iGrandDaughter->pdgId() == 111 || iGrandDaughter->pdgId() == 310 || iGrandDaughter->pdgId() == 130) //checks for rho-> pi_0, KS_0, or KL_0
          nPi0++;
        if(fabs(iGrandDaughter->pdgId()) == 211 || fabs(iGrandDaughter->pdgId()) == 321)//checks for rho->+-pi or +- kaon
          nProngs++;
    	iDaughters2++;
      }//while
    }//iffabs
    iDaughters++;
  }//while

  if(nProngs == 1 && nPi0 == 0)
    return 1;
  if(nProngs == 1 && nPi0 == 1)
    return 2;
  if(nProngs == 1 && nPi0 == 2)
    return 3;
  if(nProngs == 3 && nPi0 == 0)
    return 4;
  else
    return 5;
}//tauDecayMode

//
//
//Takes a tau Ref and returns the Decay Products
//
//
reco::LeafCandidate::LorentzVector VariousFunctions::sumTauP4(const reco::GenParticleRef& particleRef, const int decayMode, const bool pi0_decay)
{
 reco::LeafCandidate::LorentzVector tauP4; 
  
  if(decayMode == 1)
  {
    reco::GenParticleRef prong1;
    if(fabs(particleRef->daughter(0)->pdgId()) == 211 || fabs(particleRef->daughter(0)->pdgId()) == 321)
      prong1 = particleRef->daughterRef(0);
    else
      prong1 = particleRef->daughterRef(1);
    tauP4 = prong1->p4(); 
  }//if decayMode ==1 

  if(decayMode == 2)
  {
    reco::GenParticleRef neutPart, prong1;
    unsigned int iDaughters = 0;
    while(iDaughters < particleRef->numberOfDaughters())
    {
      if(fabs(particleRef->daughter(iDaughters)->pdgId()) == 213 || fabs(particleRef->daughter(iDaughters)->pdgId() ) == 223 || fabs(particleRef->daughter(iDaughters)->pdgId() ) == 20213)
      {
	unsigned int iDaughters2 = 0;
	reco::GenParticleRef newRef = particleRef->daughterRef(iDaughters);
	while(iDaughters2 < newRef->numberOfDaughters())
	{
	  if(newRef->daughter(iDaughters2)->pdgId() == 111 || newRef->daughter(iDaughters2)->pdgId() == 130 || newRef->daughter(iDaughters2)->pdgId() == 310)
	    neutPart = newRef->daughterRef(iDaughters2);
	  if(fabs(newRef->daughter(iDaughters2)->pdgId()) == 211 || fabs(newRef->daughter(iDaughters2)->pdgId()) == 321)
	    prong1 = newRef->daughterRef(iDaughters2);
	  iDaughters2++;
	}//while Daughters2
      }//if rho or omega
      if(particleRef->daughter(iDaughters)->pdgId() == 111 || particleRef->daughter(iDaughters)->pdgId() == 130 || particleRef->daughter(iDaughters)->pdgId() == 310)
        neutPart = particleRef->daughterRef(iDaughters);
      if(fabs(particleRef->daughter(iDaughters)->pdgId()) == 211 || fabs(particleRef->daughter(iDaughters)->pdgId()) == 321)
        prong1 = particleRef->daughterRef(iDaughters);
      iDaughters++;
    }//while iDaughters

    if(pi0_decay && neutPart->pdgId() == 111)
      tauP4 = neutPart->daughter(0)->p4() + neutPart->daughter(1)->p4() + prong1->p4();
    else
      tauP4 = neutPart->p4() + prong1->p4();
  }//if decayMode == 2

  if(decayMode == 3)
  {
    reco::GenParticleRef neutPart1, neutPart2, prong1;
    unsigned int iDaughters = 0, check = 0; 
    while(iDaughters < particleRef->numberOfDaughters())
    {
      if(fabs(particleRef->daughter(iDaughters)->pdgId()) == 213 || fabs(particleRef->daughter(iDaughters)->pdgId() ) == 223 || fabs(particleRef->daughter(iDaughters)->pdgId() ) == 20213)
      {
        unsigned int iDaughters2 = 0;
        reco::GenParticleRef newRef = particleRef->daughterRef(iDaughters);
        while(iDaughters2 < newRef->numberOfDaughters())
        {
          if((newRef->daughter(iDaughters2)->pdgId() == 111 || newRef->daughter(iDaughters2)->pdgId() == 130 || newRef->daughter(iDaughters2)->pdgId() == 310) && check == 1 )
 	  {
            neutPart2 = newRef->daughterRef(iDaughters2);
	    check++;
	  }//if check==0
          if((newRef->daughter(iDaughters2)->pdgId() == 111 || newRef->daughter(iDaughters2)->pdgId() == 130 || newRef->daughter(iDaughters2)->pdgId() == 310) && check == 0 )
	  {
            neutPart1 = newRef->daughterRef(iDaughters2);
	    check++;
	  }//check==0
          if(fabs(newRef->daughter(iDaughters2)->pdgId()) == 211 || fabs(newRef->daughter(iDaughters2)->pdgId()) == 321)
            prong1 = newRef->daughterRef(iDaughters2);
          iDaughters2++;
        }//while Daughters2
      }//if rho or omega
      if( (particleRef->daughter(iDaughters)->pdgId() == 111 || particleRef->daughter(iDaughters)->pdgId() == 130 || particleRef->daughter(iDaughters)->pdgId() == 310 )  && check == 1 )
      {
        neutPart2 = particleRef->daughterRef(iDaughters);
	check++;
      }//check==0
      if( (particleRef->daughter(iDaughters)->pdgId() == 111 || particleRef->daughter(iDaughters)->pdgId() == 130 || particleRef->daughter(iDaughters)->pdgId() == 310 )  && check == 0 )
      { 
        neutPart1 = particleRef->daughterRef(iDaughters);
        check++;
      }//check==0
      if(fabs(particleRef->daughter(iDaughters)->pdgId()) == 211 || fabs(particleRef->daughter(iDaughters)->pdgId()) == 321)
        prong1 = particleRef->daughterRef(iDaughters);
      iDaughters++;
    }//while iDaughters

    if(pi0_decay)
    {
      if(neutPart1->pdgId() == 111 && neutPart2->pdgId() != 111)
        tauP4 = prong1->p4() + neutPart1->daughter(0)->p4() + neutPart1->daughter(1)->p4() + neutPart2->p4();
      if(neutPart2->pdgId() == 111 && neutPart1->pdgId() != 111)
	tauP4 = prong1->p4() + neutPart1->p4() + neutPart2->daughter(0)->p4() + neutPart2->daughter(1)->p4();
      if(neutPart2->pdgId() == 111 && neutPart1->pdgId() == 111)
        tauP4 = prong1->p4() + neutPart1->daughter(0)->p4() + neutPart1->daughter(1)->p4() + neutPart2->daughter(0)->p4() + neutPart2->daughter(1)->p4(); 
      else
        tauP4 = neutPart1->p4() + neutPart2->p4() + prong1->p4();
    }//if
    else
      tauP4 = neutPart1->p4() + neutPart2->p4() + prong1->p4();
  }//if decayMode ==3 


  if(decayMode == 4)
  {
    reco::GenParticleRef prong1, prong2, prong3;
    unsigned int iDaughters = 0, check = 0;
    while(iDaughters < particleRef->numberOfDaughters())
    {
      if(fabs(particleRef->daughter(iDaughters)->pdgId()) == 213 || fabs(particleRef->daughter(iDaughters)->pdgId() ) == 223  || fabs(particleRef->daughter(iDaughters)->pdgId() ) == 20213)
      {
	unsigned int iDaughters2 = 0;
	reco::GenParticleRef newRef = particleRef->daughterRef(iDaughters);
	while(iDaughters2 < newRef->numberOfDaughters())
	{  
          if( (fabs(newRef->daughter(iDaughters2)->pdgId()) == 211 || fabs(newRef->daughter(iDaughters2)->pdgId()) == 321 ) && check == 0)
          {
            prong1 = newRef->daughterRef(iDaughters2);
	    check++;
          }//if check ==0
          if( (fabs(newRef->daughter(iDaughters2)->pdgId()) == 211 || fabs(newRef->daughter(iDaughters2)->pdgId()) == 321 ) && check == 1)
          {
            prong2 = newRef->daughterRef(iDaughters2);
	    check++;
          }//check == 1
          if( (fabs(newRef->daughter(iDaughters2)->pdgId()) == 211 || fabs(newRef->daughter(iDaughters2)->pdgId()) == 321 ) && check == 2)
   	    prong3 = newRef->daughterRef(iDaughters2);
          iDaughters2++;
	}//while iDaughters2
      }//if not final state particle
      if( (fabs(particleRef->daughter(iDaughters)->pdgId()) == 211 || fabs(particleRef->daughter(iDaughters)->pdgId()) == 321 ) && check == 0)
      {
	prong1 = particleRef->daughterRef(iDaughters);
	check++;
      }// if check == 0
      if( (fabs(particleRef->daughter(iDaughters)->pdgId()) == 211 || fabs(particleRef->daughter(iDaughters)->pdgId()) == 321 ) && check == 1)
      {
	prong2 = particleRef->daughterRef(iDaughters);
	check++;
      }//check == 1
      if( (fabs(particleRef->daughter(iDaughters)->pdgId()) == 211 || fabs(particleRef->daughter(iDaughters)->pdgId()) == 321 ) && check == 2)
      {
	prong3 = particleRef->daughterRef(iDaughters);
	check++;
      }//if 2
      iDaughters++;
    }//while iDaughters

    tauP4 = prong1->p4() + prong2->p4() + prong3->p4();
  }//if decayMode == 4


  if(decayMode == 6)
  {
    reco::GenParticleRef electron;
    unsigned int nDaughters = particleRef->numberOfDaughters(), iDaughters=0;
    while (iDaughters < nDaughters)
    {
      reco::GenParticleRef iChild = particleRef->daughterRef(iDaughters);
      if(fabs(iChild->pdgId()) == 11)
	electron = iChild;
      iDaughters++;
    }//while  

    tauP4 = electron->p4();
  }//of decayMode == 6


  if(decayMode == 7)
  {
    reco::GenParticleRef muon;
    unsigned int nDaughters = particleRef->numberOfDaughters(), iDaughters=0;
    while (iDaughters < nDaughters) 
    {
      reco::GenParticleRef iChild = particleRef->daughterRef(iDaughters);
      if(fabs(iChild->pdgId()) == 13)
        muon = iChild;
      iDaughters++;
    }//while  
 
    tauP4 = muon->p4();
  }//if decayMode == 7

  return tauP4; 
}//sumDaughtersPt

//
//Calculates the dR of the diTau object
//
double VariousFunctions::getDiTauDR(const reco::GenParticleRef& tau1Ref, const reco::GenParticleRef& tau2Ref, const bool piDecay)
{
  int tau1DecayMode = VariousFunctions::tauDecayMode(tau1Ref), tau2DecayMode = VariousFunctions::tauDecayMode(tau2Ref);
  reco::LeafCandidate::LorentzVector tau1P4 = sumTauP4(tau1Ref, tau1DecayMode, piDecay), tau2P4 = VariousFunctions::sumTauP4(tau2Ref, tau2DecayMode, piDecay);
  double dPhi = reco::deltaPhi(tau1P4.Phi(), tau2P4.Phi() );
  double aDR = sqrt( (tau1P4.Eta() - tau2P4.Eta() )*(tau1P4.Eta() - tau2P4.Eta() )  +  (dPhi )*(dPhi ) );
  std::cout << "eta1= " << tau1P4.Eta() << " eta2= " << tau2P4.Eta() << " phi1= " << tau1P4.Phi() << " phi2= " << tau2P4.Phi() << " dr= " << aDR << std::endl;
  return aDR;
}//VariousFunctions::getDiTauDR

//
//Calculates the dR between a and b
//
double VariousFunctions::getABDR(const double bEta, const double bPhi, const reco::GenParticleRef& tau1Ref, const reco::GenParticleRef& tau2Ref, const bool piDecay)
{
  int tau1DecayMode = VariousFunctions::tauDecayMode(tau1Ref), tau2DecayMode = VariousFunctions::tauDecayMode(tau2Ref);
  reco::LeafCandidate::LorentzVector tau1P4 = VariousFunctions::sumTauP4(tau1Ref, tau1DecayMode, piDecay),  tau2P4 = VariousFunctions::sumTauP4(tau2Ref, tau2DecayMode, piDecay);
  reco::LeafCandidate::LorentzVector aP4 = tau1P4 + tau2P4;
  double dPhi = reco::deltaPhi(bPhi, aP4.Phi() );
  double aBDR = sqrt( (bEta - aP4.Eta() )*(bEta - aP4.Eta() )  +  (dPhi )*(dPhi ) );
  return aBDR;
}//VariousFunctions::getABDR

double* VariousFunctions::orderFour(const double o1, double const o2, const double o3, const double o4)
{
  static double a[4];
  double  o1o1 = 0, o1o2 = 0, o2o1 = 0, o2o2 = 0;
  if (o1 > o2 )
  {
    o1o1 = o1;
    o1o2 = o2;
  }//if 11 > 12
  else
  {
    o1o1 = o2;
    o1o2 = o1;
  }//else 12 > 11
  if (o3 > o4 )
  {
    o2o1 = o3;
    o2o2 = o4;
  }//if 21 > 22
  else
  {
    o2o1 = o4;
    o2o2 = o3;
  }//else 22 > 21
  if (o1o1 > o2o1)
  {
    a[0] = o1o1;
    if (o1o2 > o2o2)
    {
      a[3] = o2o2;
      if( o1o2 > o2o1 )
      {
        a[1] = o1o2;
        a[2] = o2o1;
      }//if o1o2 > o2o1
      else
      {
        a[1] = o2o1;
        a[2] = o1o2;
      }//else o2o1 > o1o2
    }//if o1o2 > o2o2
    else
    {
      a[3] = o1o2;
      a[1] = o2o1;
      a[2] = o2o2;
    }//else o2o2 > o1o2
  }// if o1o1 > o2o1
  else
  {
    a[0] = o2o1;
    if (o1o2 > o2o2)
    {
      a[3] = o2o2;
      a[1] = o1o1;
      a[2] = o1o2;
    }// o1o2 > o2o2
    else
    {
      a[3] = o1o2;
      if (o1o1 > o2o2)
      {
        a[1] = o1o1;
        a[2] = o2o2;
      }//if o1o1 >  o2o2
      else
      {
        a[1] = o2o2;
        a[2] = o1o1;
      }//else o2o2 > o1o1
    }//else o2o1 > o1o2
  }//else o1o1 < o2o1

 return a;
}//VariousFunctions::orderFour


void VariousFunctions::findAndPlotBMuons(const reco::GenParticleRef& bRef, const int treeDepth, TH1F* hist, const bool displayTree)
{
  for(unsigned int i = 0; i < bRef->numberOfDaughters(); i++)
  {
    reco::GenParticleRef childRef = bRef->daughterRef(i);
    for(int j = 0; j < treeDepth; j++)
      std::cout << "\t";
    if(displayTree)
      std::cout << "child #" << i << " has pdgId= " << childRef->pdgId() << " and has #" << childRef->numberOfDaughters() << " daughters." << std::endl;
    if(childRef->numberOfDaughters() > 1 || VariousFunctions::findIfInDaughters(bRef, 5, true))
      VariousFunctions::findAndPlotBMuons(childRef, treeDepth + 1, hist, displayTree);
    if(fabs(childRef->pdgId() ) == 15 )
    {
      std::cout << "\t\t<----BMOUN----->" << std::endl;
      hist->Fill(childRef->pt() );
    }//if pdgid == 15
  }//for i


}//VariousFunctions::FindAndPlotBMuons

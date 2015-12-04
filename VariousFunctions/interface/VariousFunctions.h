#ifndef BoostedTauAnalysis_VariousFunctions_interface_VariousFunctions_h
#define BoostedTauAnalysis_VariousFunctions_interface_VariousFunctions_h

#include <vector>
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TAxis.h"
#include <string>
#include "TH2F.h"


class VariousFunctions { 

  public:

    //Format and Draw Histograms
    static void formatAndDrawCanvasAndHist1D(TCanvas&, TH1F*, const Int_t, const Int_t, const Int_t, const Color_t, const Size_t, const Style_t,
					   const char*, const Float_t, const Float_t, const Float_t, const char*, const Float_t, const Float_t, const Float_t, const bool); 

    //Format and Draw Histograms
    static void formatAndDrawCanvasAndHist2D(TCanvas&, TH2F*, const Int_t, const Int_t, const Int_t, const Color_t, const Size_t, const Style_t,
                                           const char*, const Float_t, const Float_t, const Float_t, 
				           const char*, const Float_t, const Float_t, const Float_t,
  					   const char*, const Float_t, const Float_t, const Float_t);

    //Search for Specific Daughter with PDG ID
    static reco::GenParticleRef findDaughterInDaughters(const reco::GenParticleRef& , const double, const bool);
  
    //Search for Specific Daughter with PDG ID
    static bool findIfInDaughters(const reco::GenParticleRef& , const double, const bool);

    //Tau Decay Mode: 1= 1prong , 2= 1pr + 1pi_0 , 3= 1pr + 2pi_0 , 4=3pr , 5=other , 6=electronic , 7=muonic
    static int tauDecayMode(const reco::GenParticleRef&);

    //Sums the tau daughter's pt
    static reco::LeafCandidate::LorentzVector sumTauP4(const reco::GenParticleRef&, const int, const bool);

    //Calculates and returns the DR between the two taus
    static double getDiTauDR(const reco::GenParticleRef&, const reco::GenParticleRef&, const bool);

    //Calculates the dR between the A and the B
    static double getABDR(const double, const double, const reco::GenParticleRef&, const reco::GenParticleRef&, const bool);

    //Orders 4 given values
    static double* orderFour(const double, double const, const double, const double);

    //Finds all muons in a Ref and plots their pt in the given histogram
    static void findAndPlotBMuons(const reco::GenParticleRef&, const int, TH1F*, const bool);

};



#endif

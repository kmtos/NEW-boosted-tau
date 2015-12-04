// -*- C++ -*-
//
// Package:    TriggerObjectFilter_BBA
// Class:      TriggerObjectFilter_BBA
// 
/**\class TriggerObjectFilter_BBA TriggerObjectFilter_BBA.cc BoostedTauAnalysis/TriggerObjectFilter_BBA/src/TriggerObjectFilter_BBA.cc

 Description: matching the reco object to the trigger object

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Francesca Ricci-Tam,6 R-025,+41227672274,
//         Created:  Mon Aug 19 11:19:05 CEST 2013
// $Id$
//
//


// system include files
#include <memory>
#include <cmath>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "BoostedTauAnalysis/Common/interface/Common.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "HLTrigger/HLTanalyzers/interface/HLTInfo.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTConfigData.h"
#include "FWCore/Common/interface/TriggerNames.h"

using namespace std;
using namespace edm;
using namespace reco;
using namespace trigger;

//
// class declaration
//

template<class T>
class TriggerObjectFilter_BBA : public edm::EDFilter {
   public:
      explicit TriggerObjectFilter_BBA(const edm::ParameterSet&);
      ~TriggerObjectFilter_BBA();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginRun(const edm::Run& iRun, edm::EventSetup const& iSetup);
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
//      virtual bool endRun(edm::Run&, edm::EventSetup const&);
//      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
//      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------

  edm::InputTag mRecoObjTag_;
  edm::InputTag bRecoObjTag_;
  std::string outFileName_;
  edm::InputTag genParticleTag_;
  edm::InputTag triggerEventTag_;
  edm::InputTag triggerResultsTag_;
  double delRMatchingCut_;
  std::vector<edm::InputTag> hltTags_;
  HLTConfigProvider hltConfig_;
  edm::InputTag theRightHLTTag_;
  edm::InputTag mSubFilter_;
  edm::InputTag bSubFilter_;
  bool highestPTOnly_;
  unsigned int minNumObjsToPassFilter_;

  //pointer to output file object
  TFile* out_;

  //Histograms
  TH1F* NEvents_;
  TH1F* DiTaudR_;


};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
template<class T>
TriggerObjectFilter_BBA<T>::TriggerObjectFilter_BBA(const edm::ParameterSet& iConfig):hltConfig_()
{
   //now do what ever initialization is needed
  mRecoObjTag_ = iConfig.getParameter<edm::InputTag>("mRecoObjTag");
  bRecoObjTag_ = iConfig.getParameter<edm::InputTag>("bRecoObjTag");
  outFileName_ = iConfig.getParameter<std::string>("outFileName");
  genParticleTag_ = iConfig.getParameter<edm::InputTag>("genParticleTag");
  const edm::InputTag dTriggerEventTag("hltTriggerSummaryAOD","","HLT");
  triggerEventTag_ = iConfig.getUntrackedParameter<edm::InputTag>("triggerEventTag",dTriggerEventTag);
  const edm::InputTag dTriggerResults("TriggerResults","","HLT"); // By default, trigger results are labeled "TriggerResults" with process name "HLT" in the event.
  triggerResultsTag_ = iConfig.getUntrackedParameter<edm::InputTag>("triggerResultsTag",dTriggerResults);
  delRMatchingCut_ = iConfig.getUntrackedParameter<double>("triggerDelRMatch", 0.30);
  hltTags_ = iConfig.getParameter<std::vector<edm::InputTag> >("hltTags");
  theRightHLTTag_ = iConfig.getParameter<edm::InputTag>("theRightHLTTag");
  mSubFilter_ = iConfig.getParameter<edm::InputTag>("mSubFilter");
  bSubFilter_ = iConfig.getParameter<edm::InputTag>("bSubFilter"); //Whether using HLT trigger path name or the actual trigger filter name. Trigger path is default.
  minNumObjsToPassFilter_ = iConfig.getParameter<unsigned int>("minNumObjsToPassFilter");

  produces<edm::RefVector<std::vector<T> > >();
}

template<class T>
TriggerObjectFilter_BBA<T>::~TriggerObjectFilter_BBA()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
template<class T>
bool TriggerObjectFilter_BBA<T>::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::cout << "<--------------NEW EVENT-------------->" << std::endl;
  NEvents_->Fill(0);
  //create pointer to output collection
  std::auto_ptr<edm::RefVector<std::vector<T> > > RecoObjColl(new edm::RefVector<std::vector<T> >);

  int mIndex = 0, bIndex = 0;

  //get muon reco objects
  edm::Handle<edm::RefVector<std::vector<T> > > mRecoObjs;
  iEvent.getByLabel(mRecoObjTag_, mRecoObjs);

  //get jet reco objects
  edm::Handle<edm::RefVector<std::vector<T> > > bRecoObjs;
  iEvent.getByLabel(bRecoObjTag_, bRecoObjs);

  //Get gen particle collection
  edm::Handle<reco::GenParticleCollection> pGenParticles;
  iEvent.getByLabel(genParticleTag_, pGenParticles);

   // Trigger Info
  edm::Handle<trigger::TriggerEvent> trgEvent;
  iEvent.getByLabel(triggerEventTag_,trgEvent);
  edm::Handle<edm::TriggerResults> pTrgResults;
  iEvent.getByLabel(triggerResultsTag_, pTrgResults);
  std::map<std::string, bool> triggerInMenu;
  std::string myHLTFilter = "";

  // get names of active HLT paths in this event
  std::vector<std::string> activeHLTPathsInThisEvent = hltConfig_.triggerNames();
  // loop over active HLT paths to search for desired path
  for (std::vector<std::string>::const_iterator iHLT = activeHLTPathsInThisEvent.begin();iHLT != activeHLTPathsInThisEvent.end(); ++iHLT) 
  { // active paths loop
    for (std::vector<edm::InputTag>::const_iterator iMyHLT = hltTags_.begin(); iMyHLT != hltTags_.end(); ++iMyHLT) 
    {
      if ((*iMyHLT).label() == *iHLT) 
      {
	std::cout << "\tin  if ((*iMyHLT).label() == *iHLT) " << std::endl;
	myHLTFilter = (*iMyHLT).label();
	triggerInMenu[(*iMyHLT).label()] = true;
      }//if
    }//for iMyHLT
  } // active paths loop
  
  edm::InputTag filterTag;
  // loop over these objects to see whether they match
  const trigger::TriggerObjectCollection& TOC( trgEvent->getObjects() );

  //choose the right sub-filter depending on the HLT path name for 
  std::vector<std::string> filters;
  try { filters = hltConfig_.moduleLabels( theRightHLTTag_.label() ); }
  catch (std::exception ex) { cout << "\tbad trigger\n"; }


  //loop over filterTags of the trgEvent for the muon subfilter
  //store the position of the one that matches the right sub-filter
  for(int i=0; i != trgEvent->sizeFilters(); ++i) 
  {
    std::string label(trgEvent->filterTag(i).label() );
    if (label.find(mSubFilter_.label() ) != std::string::npos )
      mIndex = i;
  }//for i

  //loop over filterTags of the trgEvent for the bjet subfilter
  //store the position of the one that matches the right sub-filter
  for(int i=0; i != trgEvent->sizeFilters(); ++i)
  {
    std::string label(trgEvent->filterTag(i).label() );
    if (label.find(bSubFilter_.label() ) != std::string::npos )
      bIndex = i;
  }//for i

  // find how many objects there are
  const trigger::Keys& mKEYS(trgEvent->filterKeys(mIndex)), bKEYS(trgEvent->filterKeys(bIndex));
  const size_type mnK(mKEYS.size() ), bnK(bKEYS.size() );
  std::cout << "\tmIndex= " << mIndex << "  bIndex= " << bIndex << "  mnK= " << mnK << "  bnK= " << bnK << std::endl;

  std::cout << "\ttrgEvent->filterTag(mIndex).label()= " << trgEvent->filterTag(mIndex).label() << " \ttrgEvent->filterTag(bIndex).label()= " << trgEvent->filterTag(bIndex).label() << std::endl;
  //did this event fire the HLT?
  const edm::TriggerNames &trgNames = iEvent.triggerNames(*pTrgResults);
//   for (unsigned int i = 0; i < trgNames.size(); i++)
//     cout << "trgName = " << trgNames.triggerName(i) << endl; 
  const unsigned int trgIndex = trgNames.triggerIndex(myHLTFilter);
  bool firedHLT = (trgIndex < trgNames.size()) && (pTrgResults->accept(trgIndex));
  cout << "\ttrgIndex= " << trgIndex << "  and trgNames size= " << trgNames.size() << "  firedHLT= " << firedHLT << endl;

  //store ref key of each passing reco object so we can check that no reco object is written into the produced collection more than once
  std::vector<unsigned int> mPassingRecoObjRefKeys, bPassingRecoObjRefKeys;

/////////////////////////
// global Tracker Muon
/////////////////////////
   // If the event fired the HLT, loop over the trigger object collection 
//   if (firedHLT)
//   { // firedHLT, Get cut decision for each candidate, Did this candidate cause an HLT trigger?
    for(int ipart = 0; ipart != mnK; ++ipart) 
    { 
      const trigger::TriggerObject& TO = TOC[mKEYS[ipart]];	
      //save RECO objects matched to trigger objects
      for (typename edm::RefVector<std::vector<T> >::const_iterator iRecoObj =	mRecoObjs->begin(); iRecoObj != mRecoObjs->end(); ++iRecoObj) 
      {
	if ((deltaR(**iRecoObj, TO) < delRMatchingCut_) && (std::find(mPassingRecoObjRefKeys.begin(), mPassingRecoObjRefKeys.end(), iRecoObj->key()) == mPassingRecoObjRefKeys.end())) 
	{
	  RecoObjColl->push_back(*iRecoObj);
	  mPassingRecoObjRefKeys.push_back(iRecoObj->key());
	}//if
      }//for typename
    }//for ipart        
//   } //firedHLT

///////////////////////
// for FullTrigger
///////////////////////
   // If the event fired the HLT, loop over the trigger object collection 
//   if (firedHLT)
//   { // firedHLT, Get cut decision for each candidate, Did this candidate cause an HLT trigger?
    for(int ipart = 0; ipart != bnK; ++ipart)
    { 
      const trigger::TriggerObject& TO = TOC[bKEYS[ipart]];
      //save RECO objects matched to trigger objects
      for (typename edm::RefVector<std::vector<T> >::const_iterator iRecoObj = mRecoObjs->begin(); iRecoObj != mRecoObjs->end(); ++iRecoObj)
      { 
        if ((deltaR(**iRecoObj, TO) < delRMatchingCut_) && (std::find(bPassingRecoObjRefKeys.begin(), bPassingRecoObjRefKeys.end(), iRecoObj->key()) == bPassingRecoObjRefKeys.end()))
        { 
          RecoObjColl->push_back(*iRecoObj);
          bPassingRecoObjRefKeys.push_back(iRecoObj->key());
        }//if
      }//for typename
    }//for ipart        
//   } //firedHLT

  if (firedHLT)
  {
    NEvents_->Fill(1);
    for (reco::GenParticleCollection::const_iterator iGenParticle = pGenParticles->begin(); iGenParticle != pGenParticles->end(); ++iGenParticle)
    { 
      if (iGenParticle->pdgId() == 36 && fabs(iGenParticle->daughter(0)->pdgId() ) == 13 )
      {
        double dPhi = reco::deltaPhi(iGenParticle->daughter(0)->phi(), iGenParticle->daughter(1)->phi() );
        double genDiTaudR = sqrt( (iGenParticle->daughter(0)->eta() * iGenParticle->daughter(1)->eta()) * (iGenParticle->daughter(0)->eta() * iGenParticle->daughter(1)->eta() )  + dPhi*dPhi);
        DiTaudR_->Fill(genDiTaudR );
      }//
    }//for
  }//if (firedHLT)

  if (mPassingRecoObjRefKeys.size() > 0)
    NEvents_->Fill(4);
  if (mPassingRecoObjRefKeys.size() > 1)
  {
    NEvents_->Fill(4);
    NEvents_->Fill(5);
  }//if > 0
  if (bPassingRecoObjRefKeys.size() > 0)
    NEvents_->Fill(2);
  if (bPassingRecoObjRefKeys.size() > 1)
  {
    NEvents_->Fill(2);
    NEvents_->Fill(3);
  }//if > 0


  //put collection of RECO objects matched to trigger objects into the event
  iEvent.put(RecoObjColl);
  std::cout << "\tbPassingRecoObjRefKeys.size()= " << bPassingRecoObjRefKeys.size() << " \tmPassingRecoObjRefKeys.size()= " << mPassingRecoObjRefKeys.size() << std::endl;

  bool diMatchedTrigger = (mPassingRecoObjRefKeys.size() + bPassingRecoObjRefKeys.size()  >= minNumObjsToPassFilter_ && bPassingRecoObjRefKeys.size() > 1 && mPassingRecoObjRefKeys.size() > 1) ? true : false; 
  return (diMatchedTrigger);
}

// ------------ method called once each job just before starting event loop  ------------
template<class T>
void TriggerObjectFilter_BBA<T>::beginJob()
{
  std::cout << "\nTriggerObjectFilter_BBA<T>::beginJob\n" << std::endl;
  //Open output file
  out_ = new TFile(outFileName_.c_str(), "RECREATE");

  //Book histograms
  NEvents_ = new TH1F("NEvents", "", 11, -0.5, 10.5);
      NEvents_->GetXaxis()->SetBinLabel(1, "TotalEvents"); 
      NEvents_->GetXaxis()->SetBinLabel(2, "Fired HLT");
      NEvents_->GetXaxis()->SetBinLabel(3, "> 0 bMatch");
      NEvents_->GetXaxis()->SetBinLabel(4, "> 1 bMatch");
      NEvents_->GetXaxis()->SetBinLabel(5, "> 0 #muMatch");
      NEvents_->GetXaxis()->SetBinLabel(6, "> 1 #muMatch");
  DiTaudR_ = new TH1F("DiTaudR", "", 50, -3.5, 3.5);

}

// ------------ method called once each job just after ending the event loop  ------------
template<class T>
void TriggerObjectFilter_BBA<T>::endJob() 
{
  std::cout << "\nTriggerObjectFilter_BBA<T>::endJob\n" << std::endl;
  TCanvas NEventsCanvas_("NEventsCanvas","",600,600);
  TCanvas DiTaudRCanvas_("DiTaudRCanvas","",600,600);

  NEventsCanvas_.cd();
  NEvents_->Draw();
  DiTaudRCanvas_.cd();
  DiTaudR_->Draw();
  out_->cd();
  NEventsCanvas_.Write();
  DiTaudRCanvas_.Write();
  out_->Write();
  out_->Close();
}

// ------------ method called when starting to processes a run  ------------
template<class T>
void TriggerObjectFilter_BBA<T>::beginRun(const edm::Run& iRun, edm::EventSetup const& iSetup)
{ 
  bool changed_ = true;
  std::cout << "\nTriggerObjectFilter_BBA<T>::beginRun\n" << std::endl;
  if ( !hltConfig_.init(iRun,iSetup,hltTags_[0].process(),changed_) )
  {
    edm::LogError("TriggerObjectFilter_BBA") << "Error! Can't initialize HLTConfigProvider";
    throw cms::Exception("HLTConfigProvider::init() returned non 0");
  }
}

//// ------------ method called when ending the processing of a run  ------------
//template<class T>
//bool TriggerObjectFilter_BBA<T>::endRun(edm::Run&, edm::EventSetup const&)
//{
//  std::cout << "\n\n\n\n TriggerObjectFilter_BBA<T>::endRun\n\n\n" << std::endl;
//  return true;
//}

// ------------ method called when starting to processes a luminosity block  ------------
//template<class T>
//bool TriggerObjectFilter_BBA<T>::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
//{
// return true;
//}

// ------------ method called when ending the processing of a luminosity block  ------------
//template<class T>
//bool TriggerObjectFilter_BBA<T>::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
//{
//  return true;
//}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
template<class T>
void TriggerObjectFilter_BBA<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
typedef TriggerObjectFilter_BBA<reco::Muon> BBATriggerObjectFilter;
DEFINE_FWK_MODULE(BBATriggerObjectFilter);

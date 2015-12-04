// -*- C++ -*-
//
// Package:    TriggerObjectFilter
// Class:      TriggerObjectFilter
// 
/**\class TriggerObjectFilter TriggerObjectFilter.cc BoostedTauAnalysis/TriggerObjectFilter/src/TriggerObjectFilter.cc

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
class TriggerObjectFilter : public edm::EDFilter {
   public:
      explicit TriggerObjectFilter(const edm::ParameterSet&);
      ~TriggerObjectFilter();

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
  edm::InputTag triggerEventTag_;
  edm::InputTag triggerResultsTag_;
  double delRMatchingCut_;
  std::vector<edm::InputTag> hltTags_;
  HLTConfigProvider hltConfig_;
  edm::InputTag theRightHLTTag_;
  edm::InputTag globalTrkMuonSubFilter_;
  edm::InputTag diMuonSubFilter_;
  edm::InputTag mu16SubFilter_;
  bool highestPTOnly_;
  unsigned int minNumObjsToPassFilter_;
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
TriggerObjectFilter<T>::TriggerObjectFilter(const edm::ParameterSet& iConfig):hltConfig_()
{
   //now do what ever initialization is needed
  mRecoObjTag_ = iConfig.getParameter<edm::InputTag>("mRecoObjTag");
  bRecoObjTag_ = iConfig.getParameter<edm::InputTag>("bRecoObjTag");
  const edm::InputTag dTriggerEventTag("hltTriggerSummaryAOD","","HLT");
  triggerEventTag_ = iConfig.getUntrackedParameter<edm::InputTag>("triggerEventTag",dTriggerEventTag);
  const edm::InputTag dTriggerResults("TriggerResults","","HLT"); // By default, trigger results are labeled "TriggerResults" with process name "HLT" in the event.
  triggerResultsTag_ = iConfig.getUntrackedParameter<edm::InputTag>("triggerResultsTag",dTriggerResults);
  delRMatchingCut_ = iConfig.getUntrackedParameter<double>("triggerDelRMatch", 0.30);
  hltTags_ = iConfig.getParameter<std::vector<edm::InputTag> >("hltTags");
  theRightHLTTag_ = iConfig.getParameter<edm::InputTag>("theRightHLTTag");
  globalTrkMuonSubFilter_ = iConfig.getParameter<edm::InputTag>("globalTrkMuonSubFilter");
  diMuonSubFilter_ = iConfig.getParameter<edm::InputTag>("diMuonSubFilter"); //Whether using HLT trigger path name or the actual trigger filter name. Trigger path is default.
  mu16SubFilter_ = iConfig.getParameter<edm::InputTag>("mu16SubFilter"); 
  minNumObjsToPassFilter_ = iConfig.getParameter<unsigned int>("minNumObjsToPassFilter");

  produces<edm::RefVector<std::vector<T> > >();
}

template<class T>
TriggerObjectFilter<T>::~TriggerObjectFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
template<class T>
bool TriggerObjectFilter<T>::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::cout << "\n\n\n\n TriggerObjectFilter<T>::filter\n\n\n" << std::endl;
  //create pointer to output collection
  std::auto_ptr<edm::RefVector<std::vector<T> > > glbTrkRecoObjColl(new edm::RefVector<std::vector<T> >);
  std::auto_ptr<edm::RefVector<std::vector<T> > > diMuonRecoObjColl(new edm::RefVector<std::vector<T> >);
  std::auto_ptr<edm::RefVector<std::vector<T> > > mu16RecoObjColl(new edm::RefVector<std::vector<T> >);

  int glbTrkIndex = 0, diMuonIndex = 0, mu16Index = 0;

  //get muon reco objects
  edm::Handle<edm::RefVector<std::vector<T> > > mRecoObjs;
  iEvent.getByLabel(mRecoObjTag_, mRecoObjs);

  //get jet reco objects
  edm::Handle<edm::RefVector<std::vector<T> > > bRecoObjs;
  iEvent.getByLabel(bRecoObjTag_, bRecoObjs);


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
   catch (std::exception ex) { cout << "bad trigger\n"; }


   //loop over filterTags of the trgEvent for the muon subfilter
   //store the position of the one that matches the right sub-filter
   for(int i=0; i != trgEvent->sizeFilters(); ++i) 
   {
     std::string label(trgEvent->filterTag(i).label() );
     if (label.find(globalTrkMuonSubFilter_.label() ) != std::string::npos )
       glbTrkIndex = i;
   }//for i

   //loop over filterTags of the trgEvent for the bjet subfilter
   //store the position of the one that matches the right sub-filter
   for(int i=0; i != trgEvent->sizeFilters(); ++i)
   {
     std::string label(trgEvent->filterTag(i).label() );
     if (label.find(diMuonSubFilter_.label() ) != std::string::npos )
       diMuonIndex = i;
   }//for i

   //l(op over filterTags of the trgEvent for the bjet subfilter
   //store the position of the one that matches the right sub-filter
   for(int i=0; i != trgEvent->sizeFilters(); ++i)
   {
     std::string label(trgEvent->filterTag(i).label() );
     if (label.find(mu16SubFilter_.label() ) != std::string::npos )
       mu16Index = i; 
   }//for i

   // find how many objects there are
   const trigger::Keys& glbTrkKEYS(trgEvent->filterKeys(glbTrkIndex)), diMuonKEYS(trgEvent->filterKeys(diMuonIndex)), mu16KEYS(trgEvent->filterKeys(mu16Index));
   const size_type glbTrknK(glbTrkKEYS.size() ), diMuonnK(diMuonKEYS.size() ), mu16nK(mu16KEYS.size() );

   //did this event fire the HLT?
   const edm::TriggerNames &trgNames = iEvent.triggerNames(*pTrgResults);
   for (unsigned int i = 0; i < trgNames.size(); i++)
     cout << "trgName = " << trgNames.triggerName(i) << endl; 
   const unsigned int trgIndex = trgNames.triggerIndex(myHLTFilter);
   bool firedHLT = (trgIndex < trgNames.size()) && (pTrgResults->accept(trgIndex));
   cout << "trgIndex = " << trgIndex << "  and trgNames size = " << trgNames.size() << "  firedHLT = " << firedHLT << endl;

   //store ref key of each passing reco object so we can check that no reco object is written into the produced collection more than once
   std::vector<unsigned int> glbTrkPassingRecoObjRefKeys, diMuonPassingRecoObjRefKeys, mu16PassingRecoObjRefKeys;

/////////////////////////
// global Tracker Muon
/////////////////////////
   // If the event fired the HLT, loop over the trigger object collection 
//   if (firedHLT)
//   { // firedHLT, Get cut decision for each candidate, Did this candidate cause an HLT trigger?
     for(int ipart = 0; ipart != glbTrknK; ++ipart) 
     { 
       const trigger::TriggerObject& TO = TOC[glbTrkKEYS[ipart]];	
       //save RECO objects matched to trigger objects
       for (typename edm::RefVector<std::vector<T> >::const_iterator iRecoObj =	mRecoObjs->begin(); iRecoObj != mRecoObjs->end(); ++iRecoObj) 
       {
	 if ((deltaR(**iRecoObj, TO) < delRMatchingCut_) && (std::find(glbTrkPassingRecoObjRefKeys.begin(), glbTrkPassingRecoObjRefKeys.end(), iRecoObj->key()) == glbTrkPassingRecoObjRefKeys.end())) 
	 {
	   glbTrkRecoObjColl->push_back(*iRecoObj);
	   glbTrkPassingRecoObjRefKeys.push_back(iRecoObj->key());
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
     for(int ipart = 0; ipart != diMuonnK; ++ipart)
     { 
       const trigger::TriggerObject& TO = TOC[diMuonKEYS[ipart]];
       //save RECO objects matched to trigger objects
       for (typename edm::RefVector<std::vector<T> >::const_iterator iRecoObj = mRecoObjs->begin(); iRecoObj != mRecoObjs->end(); ++iRecoObj)
       { 
         if ((deltaR(**iRecoObj, TO) < delRMatchingCut_) && (std::find(diMuonPassingRecoObjRefKeys.begin(), diMuonPassingRecoObjRefKeys.end(), iRecoObj->key()) == diMuonPassingRecoObjRefKeys.end()))
         { 
           diMuonRecoObjColl->push_back(*iRecoObj);
           diMuonPassingRecoObjRefKeys.push_back(iRecoObj->key());
         }//if
       }//for typename
     }//for ipart        
//   } //firedHLT

///////////////////////////////////////
// for Single Mu 16 or Single Mu 14er
///////////////////////////////////////
   // If the event fired the HLT, loop over the trigger object collection 
//   if (firedHLT)
//   { // firedHLT, Get cut decision for each candidate, Did this candidate cause an HLT trigger?
     for(int ipart = 0; ipart != mu16nK; ++ipart)
     {
       const trigger::TriggerObject& TO = TOC[mu16KEYS[ipart]];
       //save RECO objects matched to trigger objects
       for (typename edm::RefVector<std::vector<T> >::const_iterator iRecoObj = mRecoObjs->begin(); iRecoObj != mRecoObjs->end(); ++iRecoObj)
       {
         if ((deltaR(**iRecoObj, TO) < delRMatchingCut_) && (std::find(mu16PassingRecoObjRefKeys.begin(), mu16PassingRecoObjRefKeys.end(), iRecoObj->key()) == mu16PassingRecoObjRefKeys.end()))
         {
           mu16RecoObjColl->push_back(*iRecoObj);
           mu16PassingRecoObjRefKeys.push_back(iRecoObj->key());
         }//if
       }//for typename
     }//for ipart        
//   } //firedHLT


   //put collection of RECO objects matched to trigger objects into the event
   iEvent.put(diMuonRecoObjColl);


   bool diMatchedTrigger = false;
   return (diMatchedTrigger );
}

// ------------ method called once each job just before starting event loop  ------------
template<class T>
void TriggerObjectFilter<T>::beginJob()
{
  std::cout << "\n\n\n\n TriggerObjectFilter<T>::beginJob\n\n\n" << std::endl;
}

// ------------ method called once each job just after ending the event loop  ------------
template<class T>
void TriggerObjectFilter<T>::endJob() {
std::cout << "\n\n\n\n TriggerObjectFilter<T>::endJob\n\n\n" << std::endl;

}

// ------------ method called when starting to processes a run  ------------
template<class T>
void TriggerObjectFilter<T>::beginRun(const edm::Run& iRun, edm::EventSetup const& iSetup)
{ 
  bool changed_ = true;
  std::cout << "\n\n\n\n TriggerObjectFilter<T>::beginRun\n\n\n" << std::endl;
  if ( !hltConfig_.init(iRun,iSetup,hltTags_[0].process(),changed_) )
  {
    edm::LogError("TriggerObjectFilter") << "Error! Can't initialize HLTConfigProvider";
    throw cms::Exception("HLTConfigProvider::init() returned non 0");
  }
}

// ------------ method called when ending the processing of a run  ------------
//template<class T>
//bool TriggerObjectFilter<T>::endRun(edm::Run&, edm::EventSetup const&)
//{
//  std::cout << "\n\n\n\n TriggerObjectFilter<T>::endRun\n\n\n" << std::endl;
//  return true;
//}

// ------------ method called when starting to processes a luminosity block  ------------
//template<class T>
//bool TriggerObjectFilter<T>::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
//{
//  return true;
//}

// ------------ method called when ending the processing of a luminosity block  ------------
//template<class T>
//bool TriggerObjectFilter<T>::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
//{
//  return true;
//}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
template<class T>
void TriggerObjectFilter<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
typedef TriggerObjectFilter<reco::Muon> MuonTriggerObjectFilter;
typedef TriggerObjectFilter<reco::Photon> PhotonTriggerObjectFilter;
DEFINE_FWK_MODULE(MuonTriggerObjectFilter);
DEFINE_FWK_MODULE(PhotonTriggerObjectFilter);

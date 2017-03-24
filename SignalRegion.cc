// -*- C++ -*-
//
// Package:    test/SignalRegion
// Class:      SignalRegion
// 
/**\class SignalRegion SignalRegion.cc test/SignalRegion/plugins/SignalRegion.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mohit Gola
//         Created:  Mon, 20 Mar 2017 11:56:11 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "HLTrigger/JetMET/interface/HLTMinDPhiMETFilter.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "TTree.h"
#include "TFile.h"
#include "DataFormats/Provenance/interface/ProductID.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/Common/interface/RefHolder.h"
#include "DataFormats/Common/interface/RefVectorHolder.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <regex>

using namespace std;

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class SignalRegion : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit SignalRegion(const edm::ParameterSet&);
      ~SignalRegion();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------

  edm::EDGetTokenT<pat::JetCollection> fatjetToken_;
  edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  edm::EDGetTokenT<pat::METCollection> metToken_;
  edm::EDGetTokenT<pat::TauCollection> tauToken_;
  edm::EDGetTokenT<pat::JetCollection> jetToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleLooseIdMapToken_;
  edm::EDGetToken electronsMiniAODToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapToken_;

  edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > eleMediumIdFullInfoMapToken_;


  edm::EDGetTokenT <edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT <edm::TriggerResults> triggerObjects_;
  edm::EDGetTokenT <edm::TriggerResults> triggerPrescales_;



  TFile* file;
  TTree* tree;

  vector<double> AK8_pt;
  vector<double> AK8_eta;
  vector<double> AK8_mass;
  vector<bool> AK8_TightID;
  vector<bool> AK8_Subjet;

  vector<double> AK4_size;
  vector<double> AK4_bDisc;
  vector<bool> AK4_pileUpID_loose;
  vector <double> AK4_pt;
  vector<double> AK4_eta;
  vector<double> AK8_MET_dPhi;
  vector<bool> AK4_looseID;
  vector<bool> passEleId;
  vector<bool> Muon_Veto;
  vector<bool> Tau_Veto;
  vector<double> MET;

  vector<bool> pfMET170;
  vector<bool> pfMET90;
  vector<bool> singleEl;

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
SignalRegion::SignalRegion(const edm::ParameterSet& iConfig):
  fatjetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets"))),
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
  tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
  jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits")))

{
   //now do what ever initialization is needed
   usesResource("TFileService");

}


SignalRegion::~SignalRegion()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
SignalRegion::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;

  AK8_pt.clear();
  AK8_eta.clear();
  AK8_mass.clear();
  AK8_TightID.clear();
  AK8_Subjet.clear();
  AK4_bDisc.clear();
  AK4_pileUpID_loose.clear();
  AK4_size.clear();
  AK8_MET_dPhi.clear();
  AK4_looseID.clear();
  Muon_Veto.clear();
  Tau_Veto.clear();
  MET.clear();
  AK4_pt.clear();
  AK4_eta.clear();

  edm::Handle<pat::JetCollection> fatjets;
  iEvent.getByToken(fatjetToken_, fatjets);
  for (const pat::Jet &j : *fatjets)
    {
      bool jet_AK8_TightID = false;
      bool jet_AK8_Subjet = false;
      AK8_pt.push_back(j.pt());
      AK8_eta.push_back(j.eta());
      AK8_mass.push_back(j.mass());


      if ((j.neutralHadronEnergyFraction() < 0.90 && j.neutralEmEnergyFraction() < 0.90 && (j.chargedMultiplicity()+j.neutralMultiplicity()) > 1) &&(abs(j.eta()) <= 2.4 && j.chargedHadronEnergyFraction() > 0 && j.chargedMultiplicity() > 0 && j.chargedEmEnergyFraction() < 0.99))
	{
	  jet_AK8_TightID = true;
	}
      AK8_TightID.push_back(jet_AK8_TightID);

      
      auto wSubjets = j.subjets("SoftDrop");
      for ( auto const & iw : wSubjets )
	{
	  if(iw->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.46)
	    {

	      jet_AK8_Subjet = true;

	    }

	}
      AK8_Subjet.push_back(jet_AK8_Subjet);
    }

      edm::Handle<reco::VertexCollection> vertices;
      iEvent.getByToken(vtxToken_, vertices);
      if (vertices->empty()) return;
      const reco::Vertex &PV = vertices->front();




      edm::Handle<pat::MuonCollection> muons;
      iEvent.getByToken(muonToken_, muons);

      for (const pat::Muon &mu : *muons)
	{
	  bool Muon_Veto_Var = false;
	  if (!(mu.pt() > 10)) continue;
	  if(!(fabs(mu.eta()) < 2.4)) continue;
	  if (!((mu.pfIsolationR04().sumChargedHadronPt + max(0., mu.pfIsolationR04().sumNeutralHadronEt + mu.pfIsolationR04().sumPhotonEt - 0.5*mu.pfIsolationR04().sumPUPt))/mu.pt() < 0.4)) continue;


	  if(mu.isGlobalMuon() && mu.numberOfMatchedStations() > 1 && mu.innerTrack()->hitPattern().numberOfValidPixelHits() > 0 && mu.muonBestTrack()->dxy(PV.position()) < 0.2 && mu.muonBestTrack()->dz(PV.position())<0.5 && mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 && ((mu.innerTrack()->ptError())/(mu.innerTrack()->pt())) < 0.3 &&  mu.globalTrack()->hitPattern().numberOfValidMuonHits() > 0)
	    {

	      Muon_Veto_Var = true;

	    }

	  Muon_Veto.push_back(Muon_Veto_Var);
	}




      edm::Handle<pat::METCollection> mets;
      iEvent.getByToken(metToken_, mets);
      const pat::MET &met = mets->front();
      MET.push_back(met.pt());




      edm::Handle<pat::TauCollection> taus;
      iEvent.getByToken(tauToken_, taus);
      for (const pat::Tau &tau : *taus)
	{
	  bool Tau_Veto_Var = false;
	  bool decaymodefinding = tau.tauID("decayModeFinding") ;
	  bool byLooseCombinedIsolation =tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") ;

	  if (!(tau.pt() > 20)) continue;
	  if (!(fabs(tau.eta()) < 2.3)) continue;

	  if(decaymodefinding == true && byLooseCombinedIsolation == true)
	    {

	      Tau_Veto_Var = true;


	    }

	  Tau_Veto.push_back(Tau_Veto_Var);

	}


      edm::Handle<pat::JetCollection> jets;
      iEvent.getByToken(jetToken_, jets);
      for (const pat::Jet &j : *jets)
	{
	  AK4_pt.push_back(j.pt());
	  AK4_eta.push_back(j.eta());


	  AK4_bDisc.push_back(j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));



	  bool PileupID_loose = (j.userInt("pileupJetId:fullId")) & (1 << 2) ;
	  AK4_pileUpID_loose.push_back(PileupID_loose);


	  AK4_size.push_back(jets->size());


	  bool AK4_loose_var = false;
	  if(((j.neutralHadronEnergyFraction() < 0.99 && j.neutralEmEnergyFraction() < 0.99 && (j.chargedMultiplicity()+j.neutralMultiplicity()) > 1) && ((fabs(j.eta()) < 2.4 && j.chargedHadronEnergyFraction() > 0 && j.chargedMultiplicity() > 0 && j.chargedEmEnergyFraction()<0.99) || fabs(j.eta()) > 2.4) && fabs(j.eta()) <= 2.7) || (j.neutralHadronEnergyFraction() < 0.98 && j.neutralEmEnergyFraction() > 0.01 && j.neutralMultiplicity() > 2 && fabs(j.eta()) > 2.7 && fabs(j.eta()) <= 3.0) || (j.neutralEmEnergyFraction() < 0.9 && j.neutralMultiplicity() > 10 && fabs(j.eta()) > 3.0))
	    {
	      AK4_loose_var = true;
	    }
	  AK4_looseID.push_back(AK4_loose_var);
	  
	  
	}
  
	  if(jets->size() > 0 && mets->size() > 0)
	    {
	      double metPhi = met.phi();
	      for(size_t i = 0; i < jets->size(); ++i)
		{
		  const auto j = jets->at(i);
		  double jetPhi = j.phi();
		  double deltaPhi = std::abs(reco::deltaPhi(metPhi, jetPhi));
		  AK8_MET_dPhi.push_back(deltaPhi);
		}
	    }
	

      std::regex Pf90("(HLT_PFMET90_PFMHT90_IDTight_v)(.*)");
      std::regex Pf170("(HLT_PFMET170_NoiseCleaned_v)(.*)");

      edm::Handle <edm::TriggerResults> triggerBits;


      iEvent.getByToken(triggerBits_, triggerBits);


      const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

      pfMET90.clear();
      pfMET170.clear();

      for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i)
	{


	  bool passPFMET170 = false;
	  bool passPFMET90 = false;

	  if (std::regex_match(names.triggerName(i),Pf90))
	    {
	      if(triggerBits->accept(i))
		{

		  passPFMET90 = true;
		}
	      else
		{

		  passPFMET90 = false;
		}

	      pfMET90.push_back(passPFMET90);


	    }

	  if (std::regex_match(names.triggerName(i),Pf170))
	    {
	      if(triggerBits->accept(i))
		{
		  passPFMET170 = true;
		}
	      else
		{

		  passPFMET170 = false;
		}
	      
	      pfMET170.push_back(passPFMET170);

	    }

	}

      tree->Fill();
      




#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
SignalRegion::beginJob()
{
  file = new TFile("M_zp_600_MA0_300.root","RECREATE");
  tree = new TTree("T","Signal Region");
  tree->Branch("AK8_TightID",&AK8_TightID);
  tree->Branch("AK8_Subjet",&AK8_Subjet);
  tree->Branch("AK8_pt",&AK8_pt);
  tree->Branch("AK8_eta",&AK8_eta);
  tree->Branch("AK8_mass",&AK8_mass);
  tree->Branch("AK8_MET_dPhi",&AK8_MET_dPhi);
  tree->Branch("AK4_looseID",&AK4_looseID);
  tree->Branch("AK4_size",&AK4_size);
  tree->Branch("AK4_pileUpID_loose",&AK4_pileUpID_loose);
  tree->Branch("AK4_bDisc",&AK4_bDisc);
  tree->Branch("Tau_Veto",&Tau_Veto);
  tree->Branch("Muon_Veto",&Muon_Veto);
  tree->Branch("MET",&MET);
  tree->Branch("AK4_pt",&AK4_pt);
  tree->Branch("AK4_eta",&AK4_eta);
  tree->Branch("passEleId",&passEleId);
  tree->Branch("pfMET90",&pfMET90);
  tree->Branch("pfMET170",&pfMET170);
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SignalRegion::endJob() 
{
  file->Write();
  file->Close();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SignalRegion::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SignalRegion);

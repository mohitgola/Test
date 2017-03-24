import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use                                                                                                                     
    fileNames = cms.untracked.vstring(
#        'file:myfile.root'                                                                                                                                                          

#        '/store/mc/RunIISpring16MiniAODv2/ZprimeToA0hToA0chichihbb_2HDM_MZp-1000_MA0-300_13TeV-madgraph/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/1C0089C1-163B-E611-9931-FA163EA5733B.root',                                                                                                                                            
#        '/store/mc/RunIISpring16MiniAODv2/ZprimeToA0hToA0chichihbb_2HDM_MZp-1000_MA0-300_13TeV-madgraph/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/DE57921B-173B-E611-BFF8-FA163EAE0FC1.root',                                                                                                                                            
#        '/store/mc/RunIISpring16MiniAODv2/ZprimeToA0hToA0chichihbb_2HDM_MZp-1000_MA0-300_13TeV-madgraph/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/F0428724-173B-E611-9DAF-FA163E49EEF9.root'                                                                                                                                             

        '/store/mc/RunIISpring16MiniAODv2/ZprimeToA0hToA0chichihbb_2HDM_MZp-600_MA0-300_13TeV-madgraph/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/10000/22700D6F-C43A-E611-A4FB-008CFA0516BC.root',
        '/store/mc/RunIISpring16MiniAODv2/ZprimeToA0hToA0chichihbb_2HDM_MZp-600_MA0-300_13TeV-madgraph/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/10000/28BC187B-C43A-E611-9901-549F35AF4517.root',
        '/store/mc/RunIISpring16MiniAODv2/ZprimeToA0hToA0chichihbb_2HDM_MZp-600_MA0-300_13TeV-madgraph/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/10000/3CE70E51-C83A-E611-ADCA-001E67D5D8EF.root',
        '/store/mc/RunIISpring16MiniAODv2/ZprimeToA0hToA0chichihbb_2HDM_MZp-600_MA0-300_13TeV-madgraph/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/10000/88E7D87A-C43A-E611-8091-008CFA582D78.root'





    )
)

process.load("RecoBTag.SecondaryVertex.pfBoostedDoubleSecondaryVertexAK8BJetTags_cfi")
process.load("RecoBTag.Configuration.RecoBTag_cff")


process.demo = cms.EDAnalyzer('SignalRegion',
                              fatjets = cms.InputTag("slimmedJetsAK8"),
                              vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                              mets = cms.InputTag("slimmedMETs"),
                              muons = cms.InputTag("slimmedMuons"),
                              taus = cms.InputTag("slimmedTaus"),
                              jets = cms.InputTag("slimmedJets"),
                              bits = cms.InputTag("TriggerResults","","HLT2")


)


process.p = cms.Path(process.demo)

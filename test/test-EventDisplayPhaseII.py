import FWCore.ParameterSet.Config as cms

process = cms.Process("test")

process.load('Configuration.StandardSequences.Services_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.EventContent.EventContent_cff')
process.MessageLogger.categories = cms.untracked.vstring('L1EGRateStudies', 'FwkReport')
process.MessageLogger.cerr.FwkReport = cms.untracked.PSet(
   reportEvery = cms.untracked.int32(1)
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring('file:/hdfs/store/mc/TTI2023Upg14D/SingleTauOneProngFlatPt10To100/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/06ABA3EB-FCE6-E311-85A2-0026189438BD.root')
)

# All this stuff just runs the various EG algorithms that we are studying
                         
# ---- Global Tag :
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'PH2_1K_FB_V3::All', '')

process.load('Configuration.Geometry.GeometryExtended2023TTIReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
process.load('Geometry.TrackerGeometryBuilder.StackedTrackerGeometry_cfi')
process.load('IOMC.EventVertexGenerators.VtxSmearedHLLHC_cfi')
process.load('IOMC.EventVertexGenerators.VtxSmearedHLLHC_cfi')

process.load('Configuration/StandardSequences/L1HwVal_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load("SLHCUpgradeSimulations.L1CaloTrigger.SLHCCaloTrigger_cff")

# bug fix for missing HCAL TPs in MC RAW
from SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff import HcalTPGCoderULUT
HcalTPGCoderULUT.LUTGenerationMode = cms.bool(True)
process.valRctDigis.hcalDigis             = cms.VInputTag(cms.InputTag('valHcalTriggerPrimitiveDigis'))
process.L1CaloTowerProducer.HCALDigis =  cms.InputTag("valHcalTriggerPrimitiveDigis")

process.slhccalo = cms.Path( process.RawToDigi + process.valHcalTriggerPrimitiveDigis+process.SLHCCaloTrigger)

# run L1Reco to produce the L1EG objects corresponding
# to the current trigger
process.load('Configuration.StandardSequences.L1Reco_cff')
process.L1Reco = cms.Path( process.l1extraParticles )

# producer for UCT2015 / Stage-1 trigger objects
process.load("L1Trigger.UCT2015.emulationMC_cfi")
process.load("L1Trigger.UCT2015.uctl1extraparticles_cfi")
process.pUCT = cms.Path(
    process.emulationSequence *
    process.uct2015L1Extra
)

# --------------------------------------------------------------------------------------------
#
# ----    Produce the L1EGCrystal clusters (code of Sasha Savin)

# first you need the ECAL RecHIts :
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.reconstruction_step = cms.Path( process.calolocalreco )

process.L1EGammaCrystalsProducer = cms.EDProducer("L1EGCrystalClusterProducer",
   EtminForStore = cms.double(0.),
   debug = cms.untracked.bool(False),
   useECalEndcap = cms.bool(True)
)
process.pSasha = cms.Path( process.L1EGammaCrystalsProducer )

# --------------------------------------------------------------------------------------------
#
# ----  Match the L1EG stage-2 objects created by the SLHCCaloTrigger sequence above
#	with the crystal-level clusters.
#	This produces a new collection of L1EG objects, starting from the original
#	L1EG collection. The eta and phi of the L1EG objects is corrected using the
#	information of the xtal level clusters.

process.l1ExtraCrystalProducer = cms.EDProducer("L1ExtraCrystalPosition",
   eGammaSrc = cms.InputTag("SLHCL1ExtraParticles","EGamma"),
   eClusterSrc = cms.InputTag("L1EGammaCrystalsProducer","EGCrystalCluster")
)
process.egcrystal_producer = cms.Path(process.l1ExtraCrystalProducer)


# ----------------------------------------------------------------------------------------------
# 
# Do offline reconstruction step to get cluster pt

process.load('RecoEcal.Configuration.RecoEcal_cff')
process.ecalClusters = cms.Path(process.ecalClustersNoPFBox)


# ---------------------------------------------------------------------------
#
# --- Create the collection of special tracks for electrons
#

process.load("SLHCUpgradeSimulations.L1TrackTrigger.L1TrackingSequence_cfi")
process.pTracking = cms.Path( process.ElectronTrackingSequence )

# --- Produce the L1TkPrimaryVertex

process.load("SLHCUpgradeSimulations.L1TrackTrigger.L1TkPrimaryVertexProducer_cfi")
process.pL1TkPrimaryVertex = cms.Path( process.L1TkPrimaryVertex )


############################################################
# pixel stuff, uncomment if needed
############################################################
# from ./SLHCUpgradeSimulations/L1TrackTrigger/test/L1TrackNtupleMaker_cfg.py

#BeamSpotFromSim =cms.EDProducer("VtxSmeared")
#process.pBeamSpot = cms.Path( process.BeamSpotFromSim )
#
## pixel additions
#process.load('Configuration.StandardSequences.RawToDigi_cff')
#process.load('Configuration.StandardSequences.Reconstruction_cff')
#process.load('SimGeneral.MixingModule.mixNoPU_cfi')
#
#
#from RecoLocalTracker.SiPixelRecHits.SiPixelRecHits_cfi import *
#process.siPixelRecHits = siPixelRecHits
#
#process.L1PixelTrackFit = cms.EDProducer("L1PixelTrackFit")
#process.pixTrk = cms.Path(process.siPixelRecHits +  process.L1PixelTrackFit)
#
#process.pixRec = cms.Path(
#    process.RawToDigi+
#    process.siPixelRecHits
#)
#
#process.raw2digi_step = cms.Path(process.RawToDigi)



# ----------------------------------------------------------------------------------------------
# 
# Analyzer starts here
process.tauAnalyzer = cms.EDAnalyzer('L1TEventDisplayGenerator',
                                     L1TrackInputTag = cms.InputTag("TTTracksFromPixelDigisLargerPhi","Level1TTTracks"),
                                     L1TrackPrimaryVertexTag = cms.InputTag("L1TkPrimaryVertex"),
                                     ecalDigis = cms.InputTag("ecalDigis:EcalTriggerPrimitives"),#simEcalTriggerPrimitiveDigis
                                     hcalDigis = cms.InputTag("valHcalTriggerPrimitiveDigis"),
                                     vertices = cms.InputTag("vertices"),
                                     genParticles = cms.InputTag("genParticles"),
                                     L1CrystalClustersInputTag = cms.InputTag("L1EGammaCrystalsProducer","EGCrystalCluster"),
                                     emCrystals = cms.InputTag("L1EGammaCrystalsProducer","EGammaCrystal"),
                                     genMatchDeltaRcut = cms.untracked.double(0.25),
                                     recoPtCut = cms.double(5),
                                     folderName = cms.untracked.string("dummy")
                                     )

process.analyzer = cms.EDAnalyzer('L1EGRateStudies',
   L1EGammaInputTags = cms.VInputTag(
      # Old stage-2 trigger
      cms.InputTag("SLHCL1ExtraParticles","EGamma"),
      # 'dynamic clustering'
      cms.InputTag("SLHCL1ExtraParticlesNewClustering","EGamma"),
      # Run 1 algo.
      cms.InputTag("l1extraParticles", "Isolated"),
      cms.InputTag("l1extraParticles", "NonIsolated"),
      # UCT alg.
      cms.InputTag("l1extraParticlesUCT", "Isolated"),
      cms.InputTag("l1extraParticlesUCT", "NonIsolated"),
      # Crystal-level algo.
      cms.InputTag("L1EGammaCrystalsProducer","EGammaCrystal")
   ),
   L1CrystalClustersInputTag = cms.InputTag("L1EGammaCrystalsProducer","EGCrystalCluster"),
   OfflineRecoClustersInputTag = cms.InputTag("correctedHybridSuperClusters"),
   L1TrackInputTag = cms.InputTag("TTTracksFromPixelDigisLargerPhi","Level1TTTracks"),
   #L1PixelTrackInputTag = cms.InputTag("L1PixelTrackFit","Level1PixelTracks"),
   L1TrackPrimaryVertexTag = cms.InputTag("L1TkPrimaryVertex"),
   doEfficiencyCalc = cms.untracked.bool(True),
   useOfflineClusters = cms.untracked.bool(False),
   useEndcap = cms.untracked.bool(False),
   turnOnThresholds = cms.untracked.vint32(20, 30, 16),
   histogramBinCount = cms.untracked.int32(60),
   histogramRangeLow = cms.untracked.double(0),
   histogramRangeHigh = cms.untracked.double(50),
   histogramEtaBinCount = cms.untracked.int32(20),
   genMatchDeltaRcut = cms.untracked.double(0.25),
   genMatchRelPtcut = cms.untracked.double(0.5)
)

process.load("SLHCUpgradeSimulations.L1EGRateStudies.L1EGCrystalsHeatMap_cff")
process.L1EGCrystalsHeatMap.saveAllClusters = cms.untracked.bool(True)
#process.panalyzer = cms.Path(process.analyzer+process.L1EGCrystalsHeatMap)

process.panalyzer = cms.Path(process.tauAnalyzer)

process.TFileService = cms.Service("TFileService", 
   fileName = cms.string("eventDisplay-PhaseII.root"), 
   closeFileFast = cms.untracked.bool(True)
)


process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("fullDump-eventDisplay.root"),
    outputCommands = cms.untracked.vstring('keep *') #'keep *_*_*_L1TCaloSummaryTest')
    #outputCommands = cms.untracked.vstring('drop *', 'keep *_l1tCaloLayer1Digis_*_*, keep *_*_*_L1TCaloSummaryTest' )
)

process.e = cms.EndPath(process.out)

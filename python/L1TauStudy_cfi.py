import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('L1TauStudy',
                      L1TrackInputTag = cms.InputTag("TTTracksFromPixelDigisLargerPhi","Level1TTTracks"),
                      L1TrackPrimaryVertexTag = cms.InputTag("L1TkPrimaryVertex"),
                      ecalDigis = cms.InputTag("hcalDigis"),
                      hcalDigis = cms.InputTag("ecalDigis"),
                      vertices = cms.InputTag("vertices"),
                      genParticles = cms.InputTag("genParticles"),
                      genMatchDeltaRcut = cms.untracked.double(0.25)
                      )

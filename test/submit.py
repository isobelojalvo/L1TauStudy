process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        $inputFileNames
                                      )
    )


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("$outputFileName")
)

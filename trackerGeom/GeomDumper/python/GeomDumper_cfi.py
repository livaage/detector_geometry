import FWCore.ParameterSet.Config as cms 

ntuplizer = cms.EDFilter(
    "Geom",
    verbose = cms.int32(0),
)

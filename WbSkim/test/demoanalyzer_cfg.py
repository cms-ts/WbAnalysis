import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.cerr.threshold = cms.untracked.string("WARNING")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
	# replace 'myfile.root' with the source file you want to use
        fileNames = cms.untracked.vstring(
		'file:patTuple.root'
        )
)

process.demoEle = cms.EDProducer('WbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT  = cms.untracked.string("ee"),
	lepton  = cms.untracked.string("electron")
)

process.demoEleQCD = cms.EDProducer('WbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT  = cms.untracked.string("ee"),
	lepton  = cms.untracked.string("electronQCD")
)

process.demoEleFWD = cms.EDProducer('WbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT  = cms.untracked.string("ee"),
	lepton  = cms.untracked.string("electronFWD")
)

process.demoElePum = cms.EDProducer('WbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT  = cms.untracked.string("ee_pum"),
	lepton  = cms.untracked.string("electron")
)

process.demoElePup = cms.EDProducer('WbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT  = cms.untracked.string("ee_pup"),
	lepton  = cms.untracked.string("electron")
)

process.demoEleUp = cms.EDProducer('WbAnalyzer',
	pileupMC = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("ee"),
	lepton  = cms.untracked.string("electron"),
	JEC     = cms.untracked.double(1)
)

process.demoEleDown = cms.EDProducer('WbAnalyzer',
	pileupMC = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("ee"),
	lepton 	= cms.untracked.string("electron"),
	JEC     = cms.untracked.double(-1)
)

process.demoMuo = cms.EDProducer('WbAnalyzer',
	pileupMC = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm"),
	lepton  = cms.untracked.string("muon")
)

process.demoMuoQCD = cms.EDProducer('WbAnalyzer',
	pileupMC = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm"),
	lepton  = cms.untracked.string("muonQCD")
)

process.demoMuoFWD = cms.EDProducer('WbAnalyzer',
	pileupMC = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm"),
	lepton  = cms.untracked.string("muonFWD")
)

process.demoMuoPum = cms.EDProducer('WbAnalyzer',
	pileupMC = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm_pum"),
	lepton  = cms.untracked.string("muon")
)

process.demoMuoPup = cms.EDProducer('WbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm_pup"),
	lepton  = cms.untracked.string("muon")
)

process.demoMuoUp = cms.EDProducer('WbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm"),
	lepton 	= cms.untracked.string("muon"),
	JEC    	= cms.untracked.double(1)
)

process.demoMuoDown = cms.EDProducer('WbAnalyzer',
	pileupMC    = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm"),
	lepton  = cms.untracked.string("muon"),
	JEC     = cms.untracked.double(-1)
)

process.demoEleBtag = cms.EDProducer('WbAnalyzer',
        pileupMC    = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
	usePartonFlavour = cms.untracked.bool(True)
)

process.demoMuoBtag = cms.EDProducer('WbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
	usePartonFlavour = cms.untracked.bool(True)
)

process.demoElePur = cms.EDProducer('WbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT  = cms.untracked.string("ee"),
	lepton  = cms.untracked.string("electron"),
	pcut = cms.untracked.bool(True)
)

process.demoMuoPur = cms.EDProducer('WbAnalyzer',
	pileupMC = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm"),
	lepton  = cms.untracked.string("muon"),
	pcut = cms.untracked.bool(True)
)

process.demoEleDR = cms.EDProducer('WbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT  = cms.untracked.string("ee"),
	lepton  = cms.untracked.string("electron"),
	useDeltaR = cms.untracked.bool(True)
)

process.demoMuoDR = cms.EDProducer('WbAnalyzer',
	pileupMC = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm"),
	lepton  = cms.untracked.string("muon"),
	useDeltaR = cms.untracked.bool(True)
)

process.demoEleJerUp = cms.EDProducer('WbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
        JER     = cms.untracked.double(1)
)

process.demoEleJerDown = cms.EDProducer('WbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
        JER     = cms.untracked.double(-1)
)

process.demoMuoJerUp = cms.EDProducer('WbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
        JER     = cms.untracked.double(1)
)

process.demoMuoJerDown = cms.EDProducer('WbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
        JER     = cms.untracked.double(-1)
)

process.TFileService = cms.Service("TFileService",
	fileName = cms.string('WbTree.root')
)
process.p = cms.Path(process.demoEle*process.demoEleQCD*process.demoEleFWD*process.demoMuo*process.demoMuoQCD*process.demoMuoFWD)
#process.p = cms.Path(process.demoEle*process.demoElePum*process.demoElePup*process.demoEleUp*process.demoEleDown*process.demoMuo*process.demoMuoPum*process.demoMuoPup*process.demoMuoUp*process.demoMuoDown*process.demoEleBtag*process.demoMuoBtag*process.demoElePur*process.demoMuoPur*process.demoEleDR*process.demoMuoDR*process.demoEleJerUp*process.demoEleJerDown*process.demoMuoJerUp*process.demoMuoJerDown)

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

process.demoEle2 = cms.EDAnalyzer('ZJetsAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("ee"),
	lepton  = cms.untracked.string("electron")
)

process.demoMuo2 = cms.EDAnalyzer('ZJetsAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm"),
	lepton  = cms.untracked.string("muon")
)

process.demoEleMuo = cms.EDProducer('WbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("em"),
	lepton  = cms.untracked.string("electron+muon")
)

process.demoEleMuoUp = cms.EDProducer('WbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("em"),
	lepton  = cms.untracked.string("electron+muon"),
	JEC     = cms.untracked.double(1)
)

process.demoEleMuoDown = cms.EDProducer('WbAnalyzer',
	pileupMC   = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("em"),
	lepton  = cms.untracked.string("electron+muon"),
	JEC     = cms.untracked.double(-1)
)

process.demoEleMuoPum = cms.EDProducer('WbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("em_pum"),
        lepton  = cms.untracked.string("electron+muon")
)

process.demoEleMuoPup = cms.EDProducer('WbAnalyzer',
        pileupMC   = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("em_pup"),
        lepton  = cms.untracked.string("electron+muon")
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

process.demoEleMuoPur = cms.EDProducer('WbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT  = cms.untracked.string("em"),
	lepton  = cms.untracked.string("electron+muon"),
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

process.demoEleMuoDR = cms.EDProducer('WbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT  = cms.untracked.string("em"),
	lepton  = cms.untracked.string("electron+muon"),
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

process.demoEleMuoJerUp = cms.EDProducer('WbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT  = cms.untracked.string("em"),
	lepton  = cms.untracked.string("electron+muon"),
	JER     = cms.untracked.double(1)
)

process.demoEleMuoJerDown = cms.EDProducer('WbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT  = cms.untracked.string("em"),
	lepton  = cms.untracked.string("electron+muon"),
	JER     = cms.untracked.double(-1)
)

process.TFileService = cms.Service("TFileService",
	fileName = cms.string('WbTree.root')
)
process.p = cms.Path(process.demoEle*process.demoMuo)
#process.p = cms.Path(process.demoEle*process.demoElePum*process.demoElePup*process.demoEleUp*process.demoEleDown*process.demoMuo*process.demoMuoPum*process.demoMuoPup*process.demoMuoUp*process.demoMuoDown*process.demoEleBtag*process.demoMuoBtag*process.demoEle2*process.demoMuo2*process.demoEleMuo*process.demoEleMuoUp*process.demoEleMuoDown*process.demoEleMuoPum*process.demoEleMuoPup*process.demoElePur*process.demoMuoPur*process.demoEleMuoPur*process.demoEleDR*process.demoMuoDR*process.demoEleMuoDR*process.demoEleJerUp*process.demoEleJerDown*process.demoMuoJerUp*process.demoMuoJerDown*process.demoEleMuoJerUp*process.demoEleMuoJerDown)
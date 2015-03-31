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

process.demoEleTOP = cms.EDProducer('WbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT  = cms.untracked.string("ee"),
	lepton  = cms.untracked.string("electronTOP")
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

process.demoEleQCDPum = cms.EDProducer('WbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("ee_pum"),
        lepton  = cms.untracked.string("electronQCD")
)

process.demoEleQCDPup = cms.EDProducer('WbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("ee_pup"),
        lepton  = cms.untracked.string("electronQCD")
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

process.demoEleQCDUp = cms.EDProducer('WbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electronQCD"),
        JEC     = cms.untracked.double(1)
)

process.demoEleQCDDown = cms.EDProducer('WbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electronQCD"),
        JEC     = cms.untracked.double(-1)
)

process.demoEleTOPUp = cms.EDProducer('WbAnalyzer',
	pileupMC = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("ee"),
	lepton  = cms.untracked.string("electronTOP"),
	JEC     = cms.untracked.double(1)
)

process.demoEleTOPDown = cms.EDProducer('WbAnalyzer',
	pileupMC = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("ee"),
	lepton 	= cms.untracked.string("electronTOP"),
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

process.demoMuoTOP = cms.EDProducer('WbAnalyzer',
	pileupMC = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm"),
	lepton  = cms.untracked.string("muonTOP")
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

process.demoMuoQCDPum = cms.EDProducer('WbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm_pum"),
        lepton  = cms.untracked.string("muonQCD")
)

process.demoMuoQCDPup = cms.EDProducer('WbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm_pup"),
        lepton  = cms.untracked.string("muonQCD")
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

process.demoMuoQCDUp = cms.EDProducer('WbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muonQCD"),
        JEC     = cms.untracked.double(1)
)

process.demoMuoQCDDown = cms.EDProducer('WbAnalyzer',
        pileupMC    = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muonQCD"),
        JEC     = cms.untracked.double(-1)
)

process.demoMuoTOPUp = cms.EDProducer('WbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm"),
	lepton 	= cms.untracked.string("muonTOP"),
	JEC    	= cms.untracked.double(1)
)

process.demoMuoTOPDown = cms.EDProducer('WbAnalyzer',
	pileupMC    = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm"),
	lepton  = cms.untracked.string("muonTOP"),
	JEC     = cms.untracked.double(-1)
)

process.demoEleBtag = cms.EDProducer('WbAnalyzer',
        pileupMC    = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
	usePartonFlavour = cms.untracked.bool(True)
)

process.demoEleQCDBtag = cms.EDProducer('WbAnalyzer',
        pileupMC    = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electronQCD"),
        usePartonFlavour = cms.untracked.bool(True)
)

process.demoMuoBtag = cms.EDProducer('WbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
	usePartonFlavour = cms.untracked.bool(True)
)

process.demoMuoQCDBtag = cms.EDProducer('WbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muonQCD"),
        usePartonFlavour = cms.untracked.bool(True)
)

process.demoElePur = cms.EDProducer('WbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT  = cms.untracked.string("ee"),
	lepton  = cms.untracked.string("electron"),
	pcut = cms.untracked.bool(True)
)

process.demoEleQCDPur = cms.EDProducer('WbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electronQCD"),
        pcut = cms.untracked.bool(True)
)

process.demoMuoPur = cms.EDProducer('WbAnalyzer',
	pileupMC = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm"),
	lepton  = cms.untracked.string("muon"),
	pcut = cms.untracked.bool(True)
)

process.demoMuoQCDPur = cms.EDProducer('WbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muonQCD"),
        pcut = cms.untracked.bool(True)
)

process.demoEleDR = cms.EDProducer('WbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT  = cms.untracked.string("ee"),
	lepton  = cms.untracked.string("electron"),
	useDeltaR = cms.untracked.bool(True)
)

process.demoEleQCDDR = cms.EDProducer('WbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electronQCD"),
        useDeltaR = cms.untracked.bool(True)
)

process.demoMuoDR = cms.EDProducer('WbAnalyzer',
	pileupMC = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm"),
	lepton  = cms.untracked.string("muon"),
	useDeltaR = cms.untracked.bool(True)
)

process.demoMuoQCDDR = cms.EDProducer('WbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muonQCD"),
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

process.demoEleQCDJerUp = cms.EDProducer('WbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electronQCD"),
        JER     = cms.untracked.double(1)
)

process.demoEleQCDJerDown = cms.EDProducer('WbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electronQCD"),
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

process.demoMuoQCDJerUp = cms.EDProducer('WbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muonQCD"),
        JER     = cms.untracked.double(1)
)

process.demoMuoQCDJerDown = cms.EDProducer('WbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muonQCD"),
        JER     = cms.untracked.double(-1)
)

process.TFileService = cms.Service("TFileService",
	fileName = cms.string('WbTree.root')
)
process.p = cms.Path(process.demoEle*process.demoEleQCD*process.demoEleFWD*process.demoEleTOP*process.demoElePum*process.demoElePup*process.demoEleQCDPum*process.demoEleQCDPup*process.demoEleUp*process.demoEleDown*process.demoEleQCDUp*process.demoEleQCDDown*process.demoEleTOPUp*process.demoEleTOPDown*process.demoMuo*process.demoMuoQCD*process.demoMuoFWD*process.demoMuoTOP*process.demoMuoPum*process.demoMuoPup*process.demoMuoQCDPum*process.demoMuoQCDPup*process.demoMuoUp*process.demoMuoDown*process.demoMuoQCDUp*process.demoMuoQCDDown*process.demoMuoTOPUp*process.demoMuoTOPDown*process.demoEleBtag*process.demoEleQCDBtag*process.demoMuoBtag*process.demoMuoQCDBtag*process.demoElePur*process.demoEleQCDPur*process.demoMuoPur*process.demoMuoQCDPur*process.demoEleDR*process.demoEleQCDDR*process.demoMuoDR*process.demoMuoQCDDR*process.demoEleJerUp*process.demoEleJerDown*process.demoEleQCDJerUp*process.demoEleQCDJerDown*process.demoMuoJerUp*process.demoMuoJerDown*process.demoMuoQCDJerUp*process.demoMuoQCDJerDown)

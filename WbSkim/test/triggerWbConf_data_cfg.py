import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import *
from PhysicsTools.PatAlgos.tools.trigTools import *
from PhysicsTools.PatAlgos.tools.pfTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *
from CommonTools.ParticleFlow.ParticleSelectors.pfSelectedMuons_cfi import pfSelectedMuons 
from PhysicsTools.PatAlgos.selectionLayer1.muonSelector_cfi import * 
from CommonTools.ParticleFlow.ParticleSelectors.pfSelectedElectrons_cfi import pfSelectedElectrons 
from PhysicsTools.PatAlgos.selectionLayer1.electronSelector_cfi import *
from PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi import *
from PhysicsTools.PatAlgos.triggerLayer1.triggerEventProducer_cfi import *
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector

switchOnTrigger(process,sequence='patDefaultSequence',hltProcess='*')

process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
process.load("JetMETCorrections.Type1MET.pfMETsysShiftCorrections_cfi")

process.goodOfflinePrimaryVertices = cms.EDFilter("PrimaryVertexObjectFilter",
		filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) ),
		src=cms.InputTag('offlinePrimaryVertices')
)

process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

process.load("PhysicsTools.PatUtils.patPFMETCorrections_cff")

process.producePatPFMETCorrections.replace(
    process.pfCandMETcorr,
    process.type0PFMEtCorrection *
    process.patPFMETtype0Corr *
    process.pfCandMETcorr
)

########### Run PF2PAT

postfix = "PFlow"
usePF2PAT(process,
	  runPF2PAT=True,
          jetAlgo='AK5', 
	  runOnMC=False, 
	  postfix=postfix,
          jetCorrections=('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']),
	  typeIMetCorrections=True,
	  pvCollection=cms.InputTag('goodOfflinePrimaryVertices')
)

getattr(process,"pfPileUp"+postfix).checkClosestZVertex = False

########### Initialize lepton and PU removal from jets

usePFnoPU     = True
useNoMuon     = True # before electron top projection
useNoElectron = True # before jet top projection
useNoJet      = True # before tau top projection
useNoTau      = True # before MET top projection

########## to turn on MET type-0+1 corrections

getattr(process,'patMETs'+postfix).metSource = cms.InputTag('patType1p2CorrectedPFMet'+postfix)

getattr(process,'patType1CorrectedPFMet'+postfix).srcType1Corrections = cms.VInputTag(
    cms.InputTag("patPFJetMETtype1p2Corr"+postfix,"type1"),
    cms.InputTag("patPFMETtype0Corr"+postfix)
)

########## to turn on MET x/y shift corrections

process.pfMEtSysShiftCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runAvsNvtx_data

process.pfType1CorrectedMet.srcType1Corrections = cms.VInputTag(
		    cms.InputTag('pfJetMETcorr', 'type1'),
		        cms.InputTag('pfMEtSysShiftCorr')  
)
process.pfType1p2CorrectedMet.srcType1Corrections = cms.VInputTag(
		    cms.InputTag('pfJetMETcorr', 'type1'),
		    cms.InputTag('pfMEtSysShiftCorr')       
)

########### mu Trigger Matching

pathTriggerMu = 'path("HLT_IsoMu24_eta2p1*")'

process.selectedTriggeredPatMuons = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                                   src     = cms.InputTag('selectedPatMuons'+postfix),
                                                   matched = cms.InputTag('patTrigger'),    # selections of trigger objects
                                                   matchedCuts = cms.string(pathTriggerMu),   # selection of matches
                                                   maxDPtRel   = cms.double(0.5), 
                                                   maxDeltaR   = cms.double(0.3),
                                                   resolveAmbiguities    = cms.bool(True),
                                                   resolveByMatchQuality = cms.bool(True)
)

process.selectedPatMuonsTriggerMatch = cms.EDProducer("PATTriggerMatchMuonEmbedder",
		     src     = cms.InputTag('selectedPatMuons'+postfix),
	             matches = cms.VInputTag('selectedTriggeredPatMuons')
)

########### mu Trigger Matching - em version

pathTriggerEleMu = 'path("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL*") || path("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*")'

process.selectedTriggeredPatMuonsEM = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                                   src     = cms.InputTag('selectedPatMuons'+postfix),
                                                   matched = cms.InputTag('patTrigger'),    # selections of trigger objects
                                                   matchedCuts = cms.string(pathTriggerEleMu),   # selection of matches
                                                   maxDPtRel   = cms.double(0.5),
                                                   maxDeltaR   = cms.double(0.3),
                                                   resolveAmbiguities    = cms.bool(True),
                                                   resolveByMatchQuality = cms.bool(True)
)

process.selectedPatMuonsTriggerMatchEM = cms.EDProducer("PATTriggerMatchMuonEmbedder",
                   src     = cms.InputTag('selectedPatMuons'+postfix),
                   matches = cms.VInputTag('selectedTriggeredPatMuonsEM')
)

############## Making Jets

process.goodJets = selectedPatJets.clone(
		src = cms.InputTag('selectedPatJets'+postfix),
                     cut = cms.string(
			'pt > 25. & abs(eta) < 2.5 &'
			'numberOfDaughters > 1 &'
			'neutralHadronEnergyFraction < 0.99 &'
			'neutralEmEnergyFraction < 0.99 &'
			'chargedEmEnergyFraction < 0.99 &'
			'chargedHadronEnergyFraction > 0 &'
			'chargedMultiplicity > 0'
		)
)

############## Making Z to mumu

process.matchedMuons0 = cms.EDProducer("MuScleFitPATMuonCorrector",
                         src = cms.InputTag("selectedPatMuonsTriggerMatch"),
                         debug = cms.bool(False),
                         identifier = cms.string("Data2012_53X_ReReco"),
                         applySmearing = cms.bool(False),
                         fakeSmearing = cms.bool(False)
)

process.matchedMuons = selectedPatMuons.clone(
		src = cms.InputTag('matchedMuons0'),
		cut = cms.string(
		        'pt > 25 & abs(eta) < 2.4 &'
		        'isGlobalMuon & isPFMuon &'
		        'globalTrack.normalizedChi2 < 10 &'
		        'track.hitPattern.trackerLayersWithMeasurement > 5 &'
		        'globalTrack.hitPattern.numberOfValidMuonHits > 0 &'
		        'innerTrack.hitPattern.numberOfValidPixelHits > 0 &'
		        'abs(dB) < 0.2 &'
		        'numberOfMatchedStations > 1 &'
		        '(pfIsolationR04().sumChargedHadronPt + max(pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - 0.5*pfIsolationR04().sumPUPt,0.0))/pt < 0.2 &'
		        'triggerObjectMatches.size > 0'
		)
)

process.zmuMatchedmuMatched = cms.EDProducer('CandViewShallowCloneCombiner',
                               decay = cms.string('matchedMuons@+ matchedMuons@-'),
                               cut   = cms.string('mass > 50.0 & mass < 200.0'),
                               name  = cms.string('Zmumatchedmumatched'),
                               roles = cms.vstring('matched1', 'matched2')
)

############## Making Z to mumu - em version

process.matchedMuons0EM = cms.EDProducer("MuScleFitPATMuonCorrector",
                         src = cms.InputTag("selectedPatMuonsTriggerMatchEM"),
                         debug = cms.bool(False),
                         identifier = cms.string("Data2012_53X_ReReco"),
                         applySmearing = cms.bool(False),
                         fakeSmearing = cms.bool(False)
)

process.matchedMuonsEM = selectedPatMuons.clone(
             src = cms.InputTag('matchedMuons0EM'),
             cut = cms.string(
                     'pt > 25 & abs(eta) < 2.4 &'
                     'isGlobalMuon & isPFMuon &'
                     'globalTrack.normalizedChi2 < 10 &'
                     'track.hitPattern.trackerLayersWithMeasurement > 5 &'
                     'globalTrack.hitPattern.numberOfValidMuonHits > 0 &'
                     'innerTrack.hitPattern.numberOfValidPixelHits > 0 &'
                     'abs(dB) < 0.2 &'
                     'numberOfMatchedStations > 1 &'
                     '(pfIsolationR04().sumChargedHadronPt + max(pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - 0.5*pfIsolationR04().sumPUPt,0.0))/pt < 0.2 &'
                     'triggerObjectMatches.size > 0'
	     )
)

############## e Trigger Matching

pathTriggerEle = 'path("HLT_Ele27_WP80*")'

process.selectedTriggeredPatElectrons = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                                   src     = cms.InputTag('selectedPatElectrons'+postfix),
                                                   matched = cms.InputTag('patTrigger'),    # selections of trigger objects
                                                   matchedCuts = cms.string(pathTriggerEle),  # selection of matches
                                                   maxDPtRel   = cms.double(0.5), 
                                                   maxDeltaR   = cms.double(0.3),
                                                   resolveAmbiguities    = cms.bool(True),
                                                   resolveByMatchQuality = cms.bool(True)
)

process.selectedPatElectronsTriggerMatch = cms.EDProducer("PATTriggerMatchElectronEmbedder",
		     src     = cms.InputTag('selectedPatElectrons'+postfix),
	             matches = cms.VInputTag('selectedTriggeredPatElectrons')
)

############## e Trigger Matching - em version

process.selectedTriggeredPatElectronsEM = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                                   src     = cms.InputTag('selectedPatElectrons'+postfix),
                                                   matched = cms.InputTag('patTrigger'),    # selections of trigger objects
                                                   matchedCuts = cms.string(pathTriggerEleMu),  # selection of matches
                                                   maxDPtRel   = cms.double(0.5),
                                                   maxDeltaR   = cms.double(0.3),
                                                   resolveAmbiguities    = cms.bool(True),
                                                   resolveByMatchQuality = cms.bool(True)
)

process.selectedPatElectronsTriggerMatchEM = cms.EDProducer("PATTriggerMatchElectronEmbedder",
                   src     = cms.InputTag('selectedPatElectrons'+postfix),
                   matches = cms.VInputTag('selectedTriggeredPatElectronsEM')
)

##############

switchOnTriggerMatching(process,triggerMatchers = ['selectedTriggeredPatMuons','selectedTriggeredPatElectrons','selectedTriggeredPatMuonsEM','selectedTriggeredPatElectronsEM'],sequence ='patDefaultSequence',hltProcess = '*')

removeCleaningFromTriggerMatching(process)

############## Making Z to ee

process.matchedElectrons = selectedPatElectrons.clone(
		     src = cms.InputTag('selectedPatElectronsTriggerMatch'),
		     cut = cms.string(
			'pt > 10 & abs(eta) < 2.4 &'
			'(('
			 'abs(superCluster.eta) < 1.442 &'
			 'abs(deltaEtaSuperClusterTrackAtVtx) < 0.004 &'
			 'abs(deltaPhiSuperClusterTrackAtVtx) < 0.06 &'
			 'sigmaIetaIeta < 0.01 &'
			 'hadronicOverEm < 0.12'
			')|('
			 'abs(superCluster.eta) > 1.566 & abs(superCluster.eta) < 2.5 &'
			 'abs(deltaEtaSuperClusterTrackAtVtx) < 0.007 &'
			 'abs(deltaPhiSuperClusterTrackAtVtx) < 0.03 &'
			 'sigmaIetaIeta < 0.03 &'
			 'hadronicOverEm < 0.10'
			')) &'
			'abs(dB) < 0.02 &'
			'abs(1./ecalEnergy - eSuperClusterOverP/ecalEnergy) < 0.05 &'
			'(chargedHadronIso + max((neutralHadronIso + photonIso - 0.5*puChargedHadronIso),0.0))/et < 0.15 &'
			'passConversionVeto &'
			'gsfTrack.trackerExpectedHitsInner.numberOfHits <= 1 &'
			'triggerObjectMatches.size > 0'
		     )
)

process.zeleMatchedeleMatched = cms.EDProducer('CandViewShallowCloneCombiner',
			            decay = cms.string('matchedElectrons@+ matchedElectrons@-'),
				    cut   = cms.string('mass > 50.0 & mass < 200.0'),
				    name  = cms.string('Zelematchedelematched'),
				    roles = cms.vstring('matched1', 'matched2')
)

############## Making Z to ee - em version

process.matchedElectronsEM = selectedPatElectrons.clone(
                   src = cms.InputTag('selectedPatElectronsTriggerMatchEM'),
                   cut = cms.string(
                      'pt > 10 & abs(eta) < 2.4 &'
                      '(('
                       'abs(superCluster.eta) < 1.442 &'
                       'abs(deltaEtaSuperClusterTrackAtVtx) < 0.004 &'
                       'abs(deltaPhiSuperClusterTrackAtVtx) < 0.06 &'
                       'sigmaIetaIeta < 0.01 &'
                       'hadronicOverEm < 0.12'
                      ')|('
                       'abs(superCluster.eta) > 1.566 & abs(superCluster.eta) < 2.5 &'
                       'abs(deltaEtaSuperClusterTrackAtVtx) < 0.007 &'
                       'abs(deltaPhiSuperClusterTrackAtVtx) < 0.03 &'
                       'sigmaIetaIeta < 0.03 &'
                       'hadronicOverEm < 0.10'
                      ')) &'
                      'abs(dB) < 0.02 &'
                      'abs(1./ecalEnergy - eSuperClusterOverP/ecalEnergy) < 0.05 &'
                      '(chargedHadronIso + max((neutralHadronIso + photonIso - 0.5*puChargedHadronIso),0.0))/et < 0.15 &'
                      'passConversionVeto &'
                      'gsfTrack.trackerExpectedHitsInner.numberOfHits <= 1 &'
                      'triggerObjectMatches.size > 0'
                   )
)

############## Making Z to em

process.zeleMatchedmuMatched = cms.EDProducer('CandViewShallowCloneCombiner',
                                  decay = cms.string('matchedElectronsEM@+ matchedMuonsEM@-'),
                                  cut   = cms.string('mass > 50.0 & mass < 200.0'),
                                  name  = cms.string('ZelematchedmumatchedEM'),
                                  roles = cms.vstring('matched1', 'matched2')
)

##############

process.GlobalTag.globaltag = 'FT53_V21A_AN6::All'
process.source = cms.Source("PoolSource",
	#fileNames = cms.untracked.vstring('/store/data/Run2012A/DoubleElectron/AOD/13Jul2012-v1/00000/FA1B4710-F3D9-E111-858E-0024E876636C.root')
	fileNames = cms.untracked.vstring('file:FA1B4710-F3D9-E111-858E-0024E876636C.root')
	#fileNames = cms.untracked.vstring('/store/data/Run2012A/MuEG/AOD/22Jan2013-v1/20000/00F0AA8F-D566-E211-9A55-BCAEC50971F9.root')
	#fileNames = cms.untracked.vstring('file:00F0AA8F-D566-E211-9A55-BCAEC50971F9.root')
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.options.wantSummary = True
process.maxEvents.input = 2000

getattr(process,"pfElectronsFromVertex"+postfix).d0Cut = 0.02
getattr(process,"pfElectronsFromVertex"+postfix).dzCut = 0.1
getattr(process,"pfMuonsFromVertex"+postfix).d0Cut = 0.2
getattr(process,"pfMuonsFromVertex"+postfix).dzCut = 0.5

getattr(process,"pfIsolatedElectrons"+postfix).isolationValueMapsCharged = cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFId'+postfix))
getattr(process,"pfIsolatedElectrons"+postfix).deltaBetaIsolationValueMap = cms.InputTag('elPFIsoValuePU03PFId'+postfix)
getattr(process,"pfIsolatedElectrons"+postfix).isolationValueMapsNeutral = cms.VInputTag(cms.InputTag('elPFIsoValueNeutral03PFId'+postfix), cms.InputTag('elPFIsoValueGamma03PFId'+postfix))
getattr(process,"patElectrons"+postfix).isolationValues = cms.PSet(
        pfChargedHadrons = cms.InputTag('elPFIsoValueCharged03PFId'+postfix),
        pfChargedAll = cms.InputTag('elPFIsoValueChargedAll03PFId'+postfix),
        pfPUChargedHadrons = cms.InputTag('elPFIsoValuePU03PFId'+postfix),
        pfNeutralHadrons = cms.InputTag('elPFIsoValueNeutral03PFId'+postfix),
        pfPhotons = cms.InputTag('elPFIsoValueGamma03PFId'+postfix)
        )

getattr(process,"pfIsolatedElectrons"+postfix).doDeltaBetaCorrection = True
getattr(process,"pfIsolatedElectrons"+postfix).isolationCut = 0.15
getattr(process,"pfIsolatedMuons"+postfix).doDeltaBetaCorrection     = True
getattr(process,"pfIsolatedMuons"+postfix).isolationCut = 0.2

getattr(process,"patJets"+postfix).addTagInfos = True

########### top projections in PF2PAT:

getattr(process,"pfNoPileUp"+postfix).enable = usePFnoPU
getattr(process,"pfNoMuon"+postfix).enable = useNoMuon
getattr(process,"pfNoJet"+postfix).enable = useNoJet
getattr(process,"pfNoTau"+postfix).enable = useNoTau
getattr(process,"pfNoElectron"+postfix).enable = useNoElectron

process.dump = cms.EDAnalyzer("EventContentAnalyzer")
process.MyProcess = cms.EDFilter('WbFilter')

process.p = cms.Path(
   process.goodOfflinePrimaryVertices *
   getattr(process,"patPF2PATSequence"+postfix) *
   process.recoTauClassicHPSSequence *
   process.goodJets *
   process.pfMEtSysShiftCorrSequence *
   process.producePFMETCorrections *
   process.patTrigger *
   process.selectedTriggeredPatMuons *
   process.selectedPatMuonsTriggerMatch *
   process.matchedMuons0 *
   process.matchedMuons *
   process.zmuMatchedmuMatched *
   process.selectedTriggeredPatElectrons *
   process.selectedPatElectronsTriggerMatch *
   process.matchedElectrons *
   process.zeleMatchedeleMatched *
   process.selectedTriggeredPatMuonsEM *
   process.selectedPatMuonsTriggerMatchEM *
   process.matchedMuons0EM *
   process.matchedMuonsEM *
   process.selectedTriggeredPatElectronsEM *
   process.selectedPatElectronsTriggerMatchEM *
   process.matchedElectronsEM *
   process.zeleMatchedmuMatched *
   process.MyProcess
   #process.dump
)

process.out.fileName = 'patTuple.root'

from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
from PhysicsTools.PatAlgos.patEventContent_cff import patExtraAodEventContent
from PhysicsTools.PatAlgos.patEventContent_cff import patTriggerEventContent
process.out.outputCommands = patEventContent
process.out.outputCommands += patExtraAodEventContent
process.out.outputCommands += patTriggerEventContent
process.out.outputCommands += [
	'keep *_addPileupInfo_*_*',
	'keep *_matchedElectrons_*_*',
	'keep *_matchedMuons_*_*',
	'keep *_matchedElectronsEM_*_*',
	'keep *_matchedMuonsEM_*_*',
	'keep *_goodJets_*_*',
	'keep *_zeleMatchedeleMatched_*_*',
	'keep *_zmuMatchedmuMatched_*_*',
	'keep *_zeleMatchedmuMatched_*_*',
	'keep *_goodOfflinePrimaryVertices_*_*'
]

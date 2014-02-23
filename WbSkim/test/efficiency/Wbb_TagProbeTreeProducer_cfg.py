import FWCore.ParameterSet.Config as cms

##                      _              _       
##   ___ ___  _ __  ___| |_ __ _ _ __ | |_ ___ 
##  / __/ _ \| '_ \/ __| __/ _` | '_ \| __/ __|
## | (_| (_) | | | \__ \ || (_| | | | | |_\__ \
##  \___\___/|_| |_|___/\__\__,_|_| |_|\__|___/
##                                              
################################################

MC_flag = True
GLOBAL_TAG = 'FT53_V21A_AN6::All'
if MC_flag:
    GLOBAL_TAG = 'START53_V27::All'
    
OUTPUT_FILE_NAME = "testNewWrite.root"

########################

##    ___            _           _      
##   |_ _|_ __   ___| |_   _  __| | ___ 
##    | || '_ \ / __| | | | |/ _` |/ _ \
##    | || | | | (__| | |_| | (_| |  __/
##   |___|_| |_|\___|_|\__,_|\__,_|\___|
##
process = cms.Process("TagProbe")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.GlobalTag.globaltag = GLOBAL_TAG
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
#process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

##   ____             _ ____                           
##  |  _ \ ___   ___ | / ___|  ___  _   _ _ __ ___ ___ 
##  | |_) / _ \ / _ \| \___ \ / _ \| | | | '__/ __/ _ \
##  |  __/ (_) | (_) | |___) | (_) | |_| | | | (_|  __/
##  |_|   \___/ \___/|_|____/ \___/ \__,_|_|  \___\___|
##  
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
        'root://xrootd.ba.infn.it//store/user/dellaric/grid/v04/DYJetsToLL/patTuple_424_3_hqf.root',
    )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )    
process.source.inputCommands = cms.untracked.vstring("keep *","drop *_MEtoEDMConverter_*_*")

##    _____ ___   ____    ____       _          
##   |_   _( _ ) |  _ \  |  _ \ __ _(_)_ __ ___ 
##     | | / _ \/\ |_) | | |_) / _` | | '__/ __|
##     | || (_>  <  __/  |  __/ (_| | | |  \__ \
##     |_| \___/\/_|     |_|   \__,_|_|_|  |___/
##                                              
##   
process.PassingWP80 = cms.EDFilter("PATElectronRefSelector",
    src = cms.InputTag("matchedElectrons"),
    cut = cms.string("")
)

process.PassingHLT = process.PassingWP80.clone()
process.PassingHLT.cut = cms.string(
    #process.PassingWP80.cut.value() +
    "triggerObjectMatches.size > 0"
    )
#  Tag & probe selection ######
process.tagHLT = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("PassingWP80 PassingWP80"), # charge coniugate states are implied
    checkCharge = cms.bool(False),                           
    cut   = cms.string("40 < mass < 1000"),
)

process.allTagsAndProbes = cms.Sequence(
    process.PassingWP80 +
    process.PassingHLT +
    process.tagHLT
)


##    __  __  ____   __  __       _       _               
##   |  \/  |/ ___| |  \/  | __ _| |_ ___| |__   ___  ___ 
##   | |\/| | |     | |\/| |/ _` | __/ __| '_ \ / _ \/ __|
##   | |  | | |___  | |  | | (_| | || (__| | | |  __/\__ \
##   |_|  |_|\____| |_|  |_|\__,_|\__\___|_| |_|\___||___/
##                                                        
process.McMatchWP80 = cms.EDProducer("MCTruthDeltaRMatcherNew",
    matchPDGId = cms.vint32(11),
    src = cms.InputTag("PassingWP80"),
    distMin = cms.double(0.3),
    matched = cms.InputTag("genParticles"),
    checkCharge = cms.bool(True)
)
process.McMatchHLT = cms.EDProducer("MCTruthDeltaRMatcherNew",
    matchPDGId = cms.vint32(11),
    src = cms.InputTag("PassingHLT"),
    distMin = cms.double(0.3),
    matched = cms.InputTag("genParticles"),
    checkCharge = cms.bool(True)
)
process.mc_sequence = cms.Sequence(
   process.McMatchWP80
)

##    ___    _       __    _   _ _   _____ 
##   |_ _|__| |      \ \  | | | | | |_   _|
##    | |/ _` |  _____\ \ | |_| | |   | |  
##    | | (_| | |_____/ / |  _  | |___| |  
##   |___\__,_|      /_/  |_| |_|_____|_|
##
##  offline selection --> HLT. First specify which quantities to store in the TP tree. 
if MC_flag:
    HLTmcTruthCommonStuff = cms.PSet(
        isMC = cms.bool(MC_flag),
        tagMatches = cms.InputTag("McMatchWP80"),
        motherPdgId = cms.vint32(22,23),
        makeMCUnbiasTree = cms.bool(MC_flag),
        checkMotherInUnbiasEff = cms.bool(MC_flag),
        mcVariables = cms.PSet(
          probe_eta = cms.string("eta"),
          probe_phi  = cms.string("phi"),
          probe_et  = cms.string("et"),
          probe_charge = cms.string("charge"),
        ),
        mcFlags     =  cms.PSet(
          probe_flag = cms.string("pt>0")
        ),      
        )
else:
     HLTmcTruthCommonStuff = cms.PSet(
         isMC = cms.bool(False)
         )
##  WP80 --> HLT
process.WP80ToHLT = cms.EDAnalyzer("TagProbeFitTreeProducer",
    HLTmcTruthCommonStuff,                                
    variables = cms.PSet(
      probe_patEle_eta = cms.string("eta"),
      probe_patEle_phi  = cms.string("phi"),
      probe_patEle_et  = cms.string("et"),
      probe_patEle_charge = cms.string("charge"),
      #probe_sc_et    = cms.string("superCluster.energy*sin(superClusterPosition.theta)"),    
      #probe_sc_eta    = cms.string("superCluster.eta"), 
      #probe_sc_phi    = cms.string("superCluster.phi"),
      probe_patEle_isEB           = cms.string("isEB"),
      probe_patEle_isEE           = cms.string("isEE"),
      probe_patEle_isGap          = cms.string("isGap"),
    ),
    ignoreExceptions =  cms.bool (False),
    addRunLumiInfo   =  cms.bool (False),
    addEventVariablesInfo   =  cms.bool (False),                                                        
    tagProbePairs = cms.InputTag("tagHLT"),
    arbitration   = cms.string("Random2"),
    flags = cms.PSet( 
        probe_passingHLT = cms.InputTag("PassingHLT")        
    ),
    probeMatches  = cms.InputTag("McMatchWP80"),
    allProbes     = cms.InputTag("PassingWP80")
)

process.tree_sequence = cms.Sequence(
    process.WP80ToHLT
)    

##    ____       _   _     
##   |  _ \ __ _| |_| |__  
##   | |_) / _` | __| '_ \ 
##   |  __/ (_| | |_| | | |
##   |_|   \__,_|\__|_| |_|
##
if MC_flag:
    process.tagAndProbe = cms.Path(
        #process.sc_sequence + process.eIDSequence + process.ele_sequence + 
        #process.ext_ToNearestJet_sequence + 
        process.allTagsAndProbes +
        process.mc_sequence + 
        process.tree_sequence
        )
else:
    process.tagAndProbe = cms.Path(
        #process.sc_sequence + process.eIDSequence + process.ele_sequence + 
        #process.ext_ToNearestJet_sequence + 
        process.allTagsAndProbes +
        process.tree_sequence
        )
    
process.TFileService = cms.Service(
    "TFileService", fileName = cms.string( OUTPUT_FILE_NAME )
    )

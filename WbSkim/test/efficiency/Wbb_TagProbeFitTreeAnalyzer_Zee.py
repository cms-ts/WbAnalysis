import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


##                      _              _       
##   ___ ___  _ __  ___| |_ __ _ _ __ | |_ ___ 
##  / __/ _ \| '_ \/ __| __/ _` | '_ \| __/ __|
## | (_| (_) | | | \__ \ || (_| | | | | |_\__ \
##  \___\___/|_| |_|___/\__\__,_|_| |_|\__|___/
##                                              
################################################


isMC = True
InputFileName = "SingleElectron_2012A_22Jan13.root"
OutputFilePrefix = "efficiency-data-"



################################################
HLTDef = "probe_passingHLT"
PDFName = "pdfSignalPlusBackground"

if isMC:
    InputFileName = "DYJetsToLL.root"
    #PDFName = ""
    OutputFilePrefix = "efficiency-mc-"
################################################

#### For data: except for HLT step
EfficiencyBinningSpecification = cms.PSet(
    #specifies what unbinned variables to include in the dataset, the mass is needed for the fit
    UnbinnedVariables = cms.vstring("mass"),
    #specifies the binning of parameters
    BinnedVariables = cms.PSet(
    probe_patEle_et = cms.vdouble( 10, 15, 20, 30, 40, 50, 1000 ),
    probe_patEle_eta = cms.vdouble( -2.5, -2.0, -1.566, -1.4442, -0.8, 0.0, 0.8, 1.4442, 1.566, 2.0, 2.5 )
    ),
    #first string is the default followed by binRegExp - PDFname pairs
    BinToPDFmap = cms.vstring(PDFName) # leave empty if you want to cut&count
)

#### For MC truth: do truth matching
EfficiencyBinningSpecificationMC = cms.PSet(
    UnbinnedVariables = cms.vstring("mass"),
    BinnedVariables = cms.PSet(
        probe_et = cms.vdouble( 10, 15, 20, 30, 40, 50, 1000 ),
        probe_eta = cms.vdouble( -2.5, -2.0, -1.566, -1.4442, -0.8, 0.0, 0.8, 1.4442, 1.566, 2.0, 2.5 ),
        mcTrue = cms.vstring("true")
        ),
    BinToPDFmap = cms.vstring()  
    )

###########################################################################################
#############################################################################################
#if isMC:
#    mcTruthModules = cms.PSet(
#        MCtruth_HLT_Ele27_WP80 = cms.PSet(
#            EfficiencyBinningSpecificationMC,
#            EfficiencyCategoryAndState = cms.vstring("probe_passingHLT","pass"),
#            )
#        )
#else:
#    mcTruthModules = cms.PSet()
###########################################################################################
###########################################################################################
    
    
    
    
############################################################################################
####### Efficiency  measurement
############################################################################################
    
process.WP80toHLT = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
                                         # IO parameters:
                                             InputFileNames = cms.vstring(InputFileName),
                                         InputDirectoryName = cms.string("WP80ToHLT"),
                                         InputTreeName = cms.string("fitter_tree"),
                                         OutputFileName = cms.string(OutputFilePrefix+"WP80toHLT.root"),
                                         #numbrer of CPUs to use for fitting
                                         NumCPU = cms.uint32(1),
                                         # specifies wether to save the RooWorkspace containing the data for each bin and
                                         # the pdf object with the initial and final state snapshots
                                         SaveWorkspace = cms.bool(True),
                                         floatShapeParameters = cms.bool(True),
                                         #fixVars = cms.vstring("mean"),
                                         
                                         # defines all the real variables of the probes available in the input tree and intended for use in the efficiencies
                                         Variables = cms.PSet(
        mass = cms.vstring("Tag-Probe Mass", "60.0", "120.0", "GeV/c^{2}"),
        probe_patEle_et = cms.vstring("Probe E_{T}", "0", "1000", "GeV/c"),
        probe_patEle_eta = cms.vstring("Probe #eta", "-2.5", "2.5", ""),                
        ),

                                         # defines all the discrete variables of the probes available in the input tree and intended for use in the efficiency calculations
                                         Categories = cms.PSet(
        mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
        probe_passingHLT = cms.vstring("probe_passingHLT", "dummy[pass=1,fail=0]"), 
        ),
                                         # defines all the PDFs that will be available for the efficiency calculations; uses RooFit's "factory" syntax;
                                         # each pdf needs to define "signal", "backgroundPass", "backgroundFail" pdfs, "efficiency[0.9,0,1]" and "signalFractionInPassing[0.9]" are used for initial values  
                                         PDFs = cms.PSet(
        pdfSignalPlusBackground = cms.vstring(
            ##     "CBExGaussShape::signalRes(mass, mean[2.0946e-01], sigma[8.5695e-04],alpha[3.8296e-04], n[6.7489e+00], sigma_2[2.5849e+00], frac[6.5704e-01])",  ### the signal function goes here
            "CBExGaussShape::signalResPass(mass, meanP[0.], sigmaP[8.5695e-04, 0., 3.],alphaP[3.8296e-04], nP[6.7489e+00], sigmaP_2[2.5849e+00], fracP[6.5704e-01])",  ### signal resolution for "pass" sample
            "CBExGaussShape::signalResFail(mass, meanF[2.0946e-01, -5., 5.], sigmaF[8.5695e-04, 0., 5.],alphaF[3.8296e-04], nF[6.7489e+00], sigmaF_2[2.5849e+00], fracF[6.5704e-01])",  ### signal resolution for "fail" sample     
            "ZGeneratorLineShape::signalPhy(mass)", ### NLO line shape
            "RooCMSShape::backgroundPass(mass, alphaPass[60.,50.,70.], betaPass[0.001, 0.,0.1], betaPass, peakPass[90.0])",
            "RooCMSShape::backgroundFail(mass, alphaFail[60.,50.,70.], betaFail[0.001, 0.,0.1], betaFail, peakFail[90.0])",
            "FCONV::signalPass(mass, signalPhy, signalResPass)",
            "FCONV::signalFail(mass, signalPhy, signalResFail)",     
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[1.0]"     
            #"Gaussian::signal(mass, mean[91.2, 89.0, 93.0], sigma[2.3, 0.5, 10.0])",
            #"RooExponential::backgroundPass(mass, cPass[-0.02,-5,0])",
            #"RooExponential::backgroundFail(mass, cFail[-0.02,-5,0])",
            #"efficiency[0.9,0,1]",
            #"signalFractionInPassing[0.9]"
            ),
        ),

                                         # defines a set of efficiency calculations, what PDF to use for fitting and how to bin the data;
                                         # there will be a separate output directory for each calculation that includes a simultaneous fit, side band subtraction and counting. 
                                         Efficiencies = cms.PSet(
        #mcTruthModules,
        HLT_Ele27_WP80 = cms.PSet(
            EfficiencyBinningSpecification,
            EfficiencyCategoryAndState = cms.vstring("probe_passingHLT","pass"),
            )
        )
                                         )

process.fit = cms.Path(
    process.WP80toHLT
    )

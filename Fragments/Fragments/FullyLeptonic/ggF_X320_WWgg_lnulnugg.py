
import FWCore.ParameterSet.Config as cms

externalLHEProducer = cms.EDProducer("ExternalLHEProducer",
    args = cms.vstring('/cvmfs/cms.cern.ch/phys_generator/gridpacks/2017/13TeV/madgraph/V5_2.4.2/GluGluToRadionToHH_M320/v1/GluGluToRadionToHH_M320_slc6_amd64_gcc481_CMSSW_7_1_30_tarball.tar.xz'),
    nEvents = cms.untracked.uint32(5000),
    numberOfParameters = cms.uint32(1),
    outputFile = cms.string('cmsgrid_final.lhe'),
    scriptName = cms.FileInPath('GeneratorInterface/LHEInterface/data/run_generic_tarball_cvmfs.sh')
)

#Link to datacards:
#https://github.com/cms-sw/genproductions/tree/master/bin/MadGraph5_aMCatNLO/cards/production/2017/13TeV/exo_diboson/Spin-0


from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunes2017.PythiaCP5Settings_cfi import *

generator = cms.EDFilter("Pythia8HadronizerFilter",
                         maxEventsToPrint = cms.untracked.int32(1),
                         pythiaPylistVerbosity = cms.untracked.int32(1),
                         filterEfficiency = cms.untracked.double(1.0),
                         pythiaHepMCVerbosity = cms.untracked.bool(False),
                         comEnergy = cms.double(13000.),
                         PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CP5SettingsBlock,
        processParameters = cms.vstring(
            '24:mMin = 0.05', # the lower limit of the allowed mass range generated by the Breit-Wigner (in GeV)
            '24:onMode = off', # Turn off all W decays 
            '24:onIfAny = 11 13', # Add W->enu, W->munu, (no tau decay)
            '25:m0 = 125.09', # Higgs mass 
            '25:onMode = off', # Turn off all H decays 
            '25:onIfMatch = 22 22', # Add H->gg decay
            '25:onIfMatch = 24 -24', # Add H->WW decay 
            'ResonanceDecayFilter:filter = on',
            'ResonanceDecayFilter:exclusive = on', #off: require at least the specified number of daughters, on: require exactly the specified number of daughters
            'ResonanceDecayFilter:eMuAsEquivalent = on', #on: treat electrons and muons as equivalent
            'ResonanceDecayFilter:allNuAsEquivalent  = on', #on: treat all three neutrino flavours as equivalent
            'ResonanceDecayFilter:udscAsEquivalent   = on', #on: treat udsc quarks as equivalent
            'ResonanceDecayFilter:mothers = 25,24', #list of mothers not specified -> count all particles in hard process+resonance decays (better to avoid specifying mothers when including leptons from the lhe in counting, since intermediate resonances are not gauranteed to appear in general
            'ResonanceDecayFilter:daughters = 11,12,11,12,22,22', # lnu,lnu,gg. Technically tau neutrinos are allowed in the final state, but there is no way to produce them here
          ),
        parameterSets = cms.vstring('pythia8CommonSettings',
                                    'pythia8CP5Settings',
                                    'processParameters'
                                    )
        )
                         )


ProductionFilterSequence = cms.Sequence(generator)
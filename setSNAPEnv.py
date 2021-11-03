#Set up useful environment variables that are common to all runs and write to file
#Typically, these will be constant for a given IPTS and a given instrument configuration
#So the recommendation is to use a file name that is just IPTS-12345 but, since everything
#is mobile, also provided the possiblity of an optional Filename tag to track changes during an IPTS
# M. Guthrie 10 Sept 2021
from mantid.simpleapi import *
from mantidqt.utils.qt.qappthreadcall import QAppThreadCall
import matplotlib.pyplot as plt
import sys
sys.path.append('/SNS/SNAP/shared/Malcolm/devel/MantidPython')
import numpy as np
import shutil
import os
import json
from datetime import datetime
import mgutilities20211012 as mg
import importlib
importlib.reload(mg)
input = QAppThreadCall(mg.workbench_input_fn)
#import datetime

#input experiment-specific variables here:
#SIPTS = 26636 #ipts for sample data
#ConfigID first digit labels detector positions, second is for wavelength setting, 
#third digit is for guide in or guide out, last digit == 0 for d-spacing, == 1 for Q-spacing
RunCycle = '2021B'
ConfigID = 1010  # 1 = "standard" DAC setting: WEST=65 deg EAST=105 deg
                # 1 = West= 90 deg East= 90 deg
nHst = 6 #set to process 6 SNAP columns as separate histograms
        # config definition must have nHst values for histogram-dependent parameters
#FileNameTag = 'config0_inQ' #optional, use when details below change within a given ipts
#
# Positional calibration for detectors
#
detCalName = 'SNAP51984.detcal'#'changelogs'  # detcal filename, expected to be in shared folder
                     # if you don't want to use a detcal, replace with 'none' (all lower case)
                     # to apply changelogs, type 'changelogs' here
#
# Vanadium run information
#

Vrun = 51940
VBrun =51941

vanHeight = 0.085#0.075 #illuminated height of vanadium pin (cm) STANDARDISE?
vanRadius = 0.035#0.035 #radius of vanadium pin (cm) STANDARISE?

#0.5044,0.5191,0.5350,0.5526,0.5936,0.6178,0.6453,0.6768,0.7134,0.7566,0.8089,0.8737,0.9571,1.0701,1.2356,1.5133,2.1401
VPeaks = '0.8089,0.8737,0.9571,1.0701,1.2356,1.5133,2.1401' #string of comma separated d-spacings vanadium peaks
VPeaks = VPeaks + ',2.064,2.220' #If necessary, can add non-V peaks here if necessary (e.g. V2O3)  START STRING WITH COMMA!


InQ = 0 #0 for d-spacing set equal to 1 for

#= 1 to check vanadium Vanadium Peak stripping settings
#= 2 to check vanadium smoothing and d-limits
#= 3 to check vanadium attenuation correction

#
#below are (mostly) fixed variables and auto-generated variable
#usually won't need to edit these
#
#this string is obsolete (delete at some point)
detLogVal = '-65.713,0.044,105.414,0.040'#det_arc1,det_lin1,det_arc2,det_lin2


instrumentTag = 'SNAP'
extn = 'nxs.h5' #changed from 'nxs' to 'nxs.h5' in later runs
datafolder = 'nexus' #changed from 'data' to 'nexus' at some point
rootDir = '/SNS/'

#sharedDir = rootDir+str(instrumentTag)+'/IPTS-'+str(SIPTS)+'/shared/' #equivalent to "shared" directory in SNS file structure  
#dataDir = rootDir+str(instrumentTag)+'/IPTS-'+str(SIPTS)+'/'+datafolder+'/' #location of sample nexus files, typically '../data' at SNS
VIPTS = GetIPTS(RunNumber=Vrun,Instrument=instrumentTag)
VBIPTS = GetIPTS(RunNumber=VBrun,Instrument=instrumentTag)
VDir = VIPTS + 'nexus/'#rootDir+str(instrumentTag)+'/IPTS-'+str(VIPTS)+'/'+datafolder+'/' #location of vanadium nexus files
VBDir = VBIPTS + 'nexus/'#rootDir+str(instrumentTag)+'/IPTS-'+str(VBIPTS)+'/'+datafolder+'/' #location of vanadium background nexus files
instFileLoc = rootDir+'/SNAP/shared/Malcolm/dataFiles/SNAP_Definition.xml'
#
# Config specific setting (my attempt to build framework for standardising)
# Actual values are not standardised yet
#
if ConfigID == 1010:
    Author = 'Malcolm'
    SetUpDate = datetime.now()
    SetUpComment = 'd-spacing: DAC setting used for ice run Aug 2021'
    CenterWavelength = 2.1
    GuideIn = 1 #
    EastAng = 105
    WestAng = 65
    TBinning = '1950,10,14750' 
    MonTBinning = '3520,30,15800' #tof_min, tof_bin, tof_max for monitor spectra !THIS ALSO NEEDS TO BE STANDARDISED!
    dmins = '0.65,0.65,0.65,0.65,0.65,0.906'
    dmaxs = '2.135,2.323,2.632,2.849,3.642,4.919'
    dBinSize ='-0.001'
    QBinSize ='0.01'
    vanPeakFWHM = '3,3,8,8,8,11' #choose value for each spectrum.
    vanPeakTol = '0.01,0.01,0.01,0.01,0.02,0.01' ##value for each column (from 1 to 6)
    vanSmoothing = '30,30,30,30,30,30' #value for each column (from 1 to 6) 
if ConfigID == 1111:
    Author = 'Malcolm'
    SetUpDate = datetime.now()
    SetUpComment = 'Q: spacing. Biancas PDF DAC 2.1 - retains high Q component'
    CenterWavelength = 1.4
    GuideIn = 1 #
    TBinning = '1950,10,14750' 
    MonTBinning = '3520,30,15800' #tof_min, tof_bin, tof_max for monitor spectra !THIS ALSO NEEDS TO BE STANDARDISED!
    dmins = '0.33,0.34,0.386,0.4,0.46,0.59'
    dmaxs = '2.19,2.45,2.84,3.0,3.83,5.18'
    dBinSize ='-0.001'
    QBinSize ='0.01'
    vanPeakFWHM = '7,3,8,8,9,11' #choose value for each spectrum.
    vanPeakTol = '0.01,0.01,0.01,0.01,0.01,0.01' ##value for each column (from 1 to 6)
    vanSmoothing = '30,30,25,20,20,10' #value for each column (from 1 to 6) 
if ConfigID == 1013:
    Author = 'Malcolm and Bianca'
    SetUpDate = datetime.now()
    SetUpComment = 'd-spacing: set up for Graphite '
    CenterWavelength = 2.1
    GuideIn = 1 #
    EastAng = 105
    WestAng = 65
    TBinning = '1950,10,14750' 
    MonTBinning = '3520,30,15800' #tof_min, tof_bin, tof_max for monitor spectra !THIS ALSO NEEDS TO BE STANDARDISED!
    dmins = '0.65,0.65,0.65,0.65,0.65,0.906'
    dmaxs = '2.135,2.323,2.632,2.849,3.642,4.919'
    dBinSize ='-0.001'
    QBinSize ='0.01'
    vanPeakFWHM = '3,3,8,8,8,11' #choose value for each spectrum.
    vanPeakTol = '0.01,0.01,0.01,0.01,0.02,0.01' ##value for each column (from 1 to 6)
    vanSmoothing = '30,30,30,30,30,30' #value for each column (from 1 to 6) 
###############################################################################################
"""
Everything below here is fixed and shouldn't be edited
"""
###############################################################################################

#generate dictionary with all the useful information above

ConfigDir = '/SNS/SNAP/shared/Calibration/' + str(ConfigID).zfill(4) +'/'

configDict = {\
"ConfigVersion":1.0,\
"Author":Author,\
"SetUpDate":SetUpDate.strftime("%d-%b-%Y (%H:%M:%S.%f)"),\
"SetUpComment":SetUpComment,\
"CenterWavelength":CenterWavelength,\
"GuideIn":GuideIn,\
"EastAng":EastAng,\
"WestAng":WestAng,\
"RunCycle":RunCycle,\
"ConfigID":ConfigID,\
"detCalName":ConfigDir + detCalName,\
"detLogVal":detLogVal,\
"VIPTS":VIPTS,\
"VBIPTS":VBIPTS,\
"Vrun":Vrun,\
"VBrun":VBrun,\
"vanPeakFWHM":vanPeakFWHM,\
"vanPeakTol":vanPeakTol,\
"vanSmoothing":vanSmoothing,\
"vanHeight":vanHeight,\
"vanRadius":vanRadius,\
"VPeaks":VPeaks,\
"instrumentTag":instrumentTag,\
"extn":extn,\
"datafolder":datafolder,\
"VDir":VDir,\
"VBDir":VBDir,\
"instFileLoc":instFileLoc,\
"TBinning":TBinning,\
"InQ":InQ,\
"dBinSize":dBinSize,\
"QBinSize":QBinSize,\
"MonTBinning":MonTBinning,\
"dmins":dmins,\
"dmaxs":dmaxs\
}

##############################################################################
# if requested, generate a vanadium correction to check that settings are good
##############################################################################
testVanConfig = 1 # used to inspect vanadium configuration

if testVanConfig ==1:
    if plt.fignum_exists('Vanadium Setup'): 
        plt.close('Vanadium Setup')
    mg.generateVCorr('',0,configDict,InQ,2)
if testVanConfig ==2:
    if plt.fignum_exists('Vanadium Setup'): 
        plt.close('Vanadium Setup')
    mg.generateVCorr('',0,configDict,InQ,3)
if testVanConfig ==3:
    if plt.fignum_exists('Vanadium Setup'): 
        plt.close('Vanadium Setup')
    mg.generateVCorr('',0,configDict,InQ,4)

##############################################################################
# Save configuration to file using vanadium run start time. So testVanConfig
# now hardwired to be == 1 
##############################################################################

ws = mtd['SNAP' + str(Vrun)]
timestr = str(ws.getRun().startTime().to_datetime64())

#ConfigFileName = sharedDir + 'IPTS' + str(SIPTS) + '_cfg' + str(ConfigID).zfill(4) + '.json'
#timestr = vanRunTime.strftime("_%d-%b-%Y-%H%M%S")

#Check if directory already exists and create if not
if not os.path.exists(ConfigDir):
    os.makedirs(ConfigDir)

ConfigFileName = ConfigDir +'config' + str(ConfigID).zfill(4) + timestr+'.json'
print(ConfigFileName)

if os.path.exists(ConfigFileName):
    tog = input('WARNING!','Config file already exists, do you want to overwrite? (y/n)','str')
    if tog.lower() == 'y':
        print('writing config file')
        with open(ConfigFileName, "w") as outfile:
            json.dump(configDict, outfile)
    else:
        pass #Default is always to not overwrite...this could be costly.
else:
    print('writing config file')
    with open(ConfigFileName, "w") as outfile:
        json.dump(configDict, outfile)

#test, read it back in

with open(ConfigFileName, "r") as json_file:
    dictIn = json.load(json_file)



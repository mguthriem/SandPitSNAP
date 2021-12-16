# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
from mantidqt.widgets.instrumentview.api import get_instrumentview
import matplotlib.pyplot as plt
import numpy as np
import sys
#import time from datetime 
# import time 
from datetime import datetime, date
from os.path import exists
import mgutilities as mg

#run = 51984
run = sys.argv[1]
print(sys.argv[2].lower)
calibrant = str(sys.argv[2]).strip().lower()
print(calibrant)
if calibrant =='diamond':
    print('Expecting a diamond powder data set')
    peaksToFit = '1.26116,1.07552'
    #Diamond peak positions from Shikata et al J. Appl. Phys 57, 111301 (2018)
    #First 3 peak are: 2.05947,1.26116,1.07552
else:
    print('Only works for diamond!')
    sys.exit()

IPTSLoc = GetIPTS(RunNumber=run,Instrument='SNAP')
inputFile = IPTSLoc + '/nexus/SNAP_%s'%(run) + '.nxs.h5'
# check file exists
if not exists(inputFile):
    print('error! input nexus file does not exist')
    sys.exit() #terminate
# load nexus file and determine instrument state
LoadEventNexus(Filename=inputFile, OutputWorkspace='calFile')
stateID = mg.getSNAPPars('calFile')

# Fit DIFCs and Zeros
SumNeighbours(InputWorkspace='calFile', OutputWorkspace='calFile_8x8', SumX=8, SumY=8)

PDCalibration(InputWorkspace='calFile_8x8', TofBinning='2500,10,14500',\
    PeakFunction='Gaussian',BackgroundType='Quadratic',\
    PeakPositions=peaksToFit, CalibrationParameters='DIFC+TZERO', OutputCalibrationTable='_cal',\
    DiagnosticWorkspaces='_diag')

# Save Calibration File to correct State folder
stateFolder = '/SNS/SNAP/shared/Calibration/%s/'%(stateID)
now = date.today()
calibFilename = stateFolder + 'SNAP_calibrate_d%s_'%(run)+ now.strftime("%Y%m%d") + '.h5'

print(calibFilename)
SaveDiffCal(CalibrationWorkspace='_cal',Filename=calibFilename)

# Check calibration by using SNAPReduce to create focused detector columns
SNAPReduce(RunNumbers=str(run),\
    Calibration='Calibration File', CalibrationFilename=calibFilename,\
    Binning='0.5,-0.002,4', GroupDetectorsBy='Column')

#Lastly clean up unused workspaces
DeleteWorkspace('_cal')
DeleteWorkspace('_cal_mask')
#DeleteWorkspace('calFile')
#DeleteWorkspace('calFile_8x8')
DeleteWorkspace('Column')
myiv = get_instrumentview('calFile')
myiv.show_view()


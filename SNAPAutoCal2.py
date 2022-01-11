# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
#from mantidqt.widgets.instrumentview.api import get_instrumentview
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from mantid.plots.utility import MantidAxType
from mantid.api import AnalysisDataService as ADS
import numpy as np
import sys
from datetime import datetime, date
from os.path import exists
import mgutilities as mg
import importlib
importlib.reload(mg)

#M. Guthrie 16 Dec 2021 added built in graphics to allow user to confirm input data quality and resultant fit

#run = 51984
run = sys.argv[1]
try: 
    print(sys.argv[2].lower)
    calibrant = str(sys.argv[2]).strip().lower()
except:
    print('WARNING: need to specify calibrant after run number')
    calibrant = input('Enter now: (e.g. diamond):')
    calibrant = str(calibrant).strip().lower()

#print('Calibrant is: ',calibrant)
if calibrant =='diamond':
    print('Expecting a diamond powder data set')
    peaksToFit = '2.05947,1.26116,1.07552'
    #Diamond peak positions from Shikata et al J. Appl. Phys 57, 111301 (2018)
    #First 3 peak are: 2.05947,1.26116,1.07552
elif calibrant == '':
    print('Must enter type of calibrant sample after run number')
    sys.exit()
else:
    print('Currently only works for diamond!')
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

#Show user plot of input data to allow inspection (and sanity check)

mg.gridPlot(['calFile_8x8'],[],[[1584,2786,7567],[16851,13808,10802]],[],[],[],'Sample spectra')
#plt.show(block=False)
tog = input('Do input data look OK to continue ([y],n):?')
tog = str(tog).strip().lower()
if tog =='n':
    sys.exit()

#print('starting to fit spectra...')

PDCalibration(InputWorkspace='calFile_8x8', TofBinning='2500,10,15000',\
    PeakFunction='Gaussian',BackgroundType='Quadratic',\
    PeakPositions=peaksToFit, CalibrationParameters='DIFC+TZERO', OutputCalibrationTable='_cal',\
    DiagnosticWorkspaces='_diag')


mg.gridPlot(['calFile_8x8','_diag_fitted'],[],[[1584,2786,7567],[16851,13808,10802]],[],[],[],'Sample spectra')

# plot chi2 of fits for all peaks
ws = mtd['_diag_fitparam']
tableData = ws.toDict()
wsindex = tableData['wsindex']
chi2 = tableData['chi2']

fig, axes = plt.subplots(num='Chi2 for all peaks fitted', subplot_kw={'projection': 'mantid'})
axes.plot(wsindex,chi2)
axes.tick_params(axis='x', which='major', **{'gridOn': False, 'tick1On': True, 'tick2On': False, 'label1On': True, 'label2On': False})
axes.tick_params(axis='y', which='major', **{'gridOn': False, 'tick1On': True, 'tick2On': False, 'label1On': True, 'label2On': False})
axes.set_xlabel('wsindex')
axes.set_ylabel('chi2')
axes.set_ylim([0.0, 2.0])
#axes.set_xlim([0,18431])
#legend = axes.legend().draggable().legend
plt.show()

# Save Calibration File to correct State folder
stateFolder = '/SNS/SNAP/shared/Calibration/%s/'%(stateID)
now = date.today()
calibFilename = stateFolder + 'temp.h5' #SNAP_calibrate_d%s_'%(run)+ now.strftime("%Y%m%d") + '.h5'

#tog = input('Output calibration file ([y],n):?')
#tog = str(tog).strip().lower()
#if tog =='n':
#    sys.exit()
#print(calibFilename)
SaveDiffCal(CalibrationWorkspace='_cal',Filename=calibFilename)

# Check calibration by using SNAPReduce to create focused detector columns
SNAPReduce(RunNumbers=str(run),\
    Calibration='Calibration File', CalibrationFilename=calibFilename,\
    Binning='0.5,-0.002,4', GroupDetectorsBy='Column')

##SCRIPT TO PLOT OUTPUT SPECTRA FROM SNAPREDUCE
SNAPRedWSName = 'SNAP_' + str(run) + '_Column_red'
SNAPredOut = ADS.retrieve(SNAPRedWSName)

fig, axes = plt.subplots(num='SNAPreduce Output', subplot_kw={'projection': 'mantid'})
axes.plot(SNAPredOut, color='#1f77b4', label='SNAPredOut: spec 1', wkspIndex=0)
axes.plot(SNAPredOut, color='#ff7f0e', label='SNAPredOut: spec 2', wkspIndex=1)
axes.plot(SNAPredOut, color='#2ca02c', label='SNAPredOut: spec 3', wkspIndex=2)
axes.plot(SNAPredOut, color='#d62728', label='SNAPredOut: spec 4', wkspIndex=3)
axes.plot(SNAPredOut, color='#9467bd', label='SNAPredOut: spec 5', wkspIndex=4)
axes.plot(SNAPredOut, color='#8c564b', label='SNAPredOut: spec 6', wkspIndex=5)
axes.tick_params(axis='x', which='major', **{'gridOn': False, 'tick1On': True, 'tick2On': False, 'label1On': True, 'label2On': False, 'size': 6, 'tickdir': 'out', 'width': 1})
axes.tick_params(axis='y', which='major', **{'gridOn': False, 'tick1On': True, 'tick2On': False, 'label1On': True, 'label2On': False, 'size': 6, 'tickdir': 'out', 'width': 1})
axes.set_title('SNAPredOut')
axes.set_xlabel('d-Spacing ($\\AA$)')
axes.set_ylabel('Counts (microAmp.hour $\\AA$)$^{-1}$')
legend = axes.legend(fontsize=8.0).draggable().legend

plt.show()


#check if all is OK then output final calibration...
calibFilename = stateFolder + 'SNAP_calibrate_d%s_'%(run)+ now.strftime("%Y%m%d") + '.h5'
tog = input('Happy with calibration? Save file ([y],n):?')
tog = str(tog).strip().lower()
if tog =='n':
    sys.exit()
else:
    print(calibFilename)
    SaveDiffCal(CalibrationWorkspace='_cal',Filename=calibFilename)

#Lastly clean up unused workspaces
DeleteWorkspace('_cal')
DeleteWorkspace('_cal_mask')
#DeleteWorkspace('calFile')
#DeleteWorkspace('calFile_8x8')
DeleteWorkspace('Column')
#myiv = get_instrumentview('calFile')
#myiv.show_view()

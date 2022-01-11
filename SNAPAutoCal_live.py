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

#M. Guthrie 16 Dec 2021: added built in graphics to allow user to confirm input data quality 
#and resultant fit

class StopExecution(Exception):
    def _render_traceback_(self):
        pass



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
    print('Expecting a diamond powder dataset')
    peaksToFit = '2.0595,1.2612,1.0755, 0.8918, 0.8183'
    #peaksToFit = '1.2612'
    #Diamond peak positions from Shikata et al J. Appl. Phys 57, 111301 (2018).
    #First 3 peak are: 2.05947,1.26116,1.07552
elif calibrant =='lab6':
    print('Expecting a LaB6 powder dataset')
    peaksToFit = '4.1561,2.9388,2.3995,2.0781,1.8587,1.6967,1.3854,1.3143,1.0390,1.0080,0.9796'
    #Peak positions are taken from Booth et al PRB 63 2243021-2243028 (2001).
elif calibrant == '':
    print('Must enter type of calibrant sample after run number')
    sys.exit()
else:
    print('Currently only works for diamond and LaB6!')
    sys.exit()

pks = peaksToFit.split(',')
nPks = len(pks)

calFile = 'snapLive4'

#IPTSLoc = GetIPTS(RunNumber=run,Instrument='SNAP')
#inputFile = IPTSLoc + '/nexus/SNAP_%s'%(run) + '.nxs.h5'
# check file exists
#if not exists(inputFile):
#    print('error! input nexus file does not exist')
#    raise StopExecution
    #sys.exit() #terminate
# load nexus file and determine instrument state
#LoadEventNexus(Filename=inputFile, OutputWorkspace='calFile')
stateID = mg.getSNAPPars(calFile)

# Fit DIFCs and Zeros
SumNeighbours(InputWorkspace=calFile, OutputWorkspace='calFile_8x8', SumX=8, SumY=8)

#Show user plot of input data to allow inspection (and sanity check)

mg.gridPlot(['calFile_8x8'],[],[[1584,2786,7567],[16851,13808,10802]],[],[],[],'Sample spectra')
#plt.show(block=False)
tog = input('Do input data look OK to continue ([y],n):?')
tog = str(tog).strip().lower()
#tog = 'y'
if tog =='n':
    print('OK. STOPPING!')
    sys.exit()
    #raise StopExecution

#print('starting to fit spectra...')

PDCalibration(InputWorkspace='calFile_8x8', TofBinning='1500,10,16000',\
    PeakFunction='Gaussian',BackgroundType='Linear',\
    PeakPositions=peaksToFit, CalibrationParameters='DIFC', OutputCalibrationTable='_cal',\
    DiagnosticWorkspaces='_diag')


mg.gridPlot(['calFile_8x8','_diag_fitted'],[],[[1584,2786,7567],[16851,13808,10802]],[],[],[],'Sample spectra')



# Save temp copy of calibration File to correct State folder.
# option to inspect and make permanent copy later
stateFolder = '/SNS/SNAP/shared/Calibration/%s/'%(stateID)
now = date.today()
calibFilename = stateFolder + 'temp.h5' #SNAP_calibrate_d%s_'%(run)+ now.strftime("%Y%m%d") + '.h5'
SaveDiffCal(CalibrationWorkspace='_cal',Filename=calibFilename)

# Check calibration by using SNAPReduce to create focused detector columns

tog = input('Run SNAPReduce to inspect columns (y,[n]):?')
tog = str(tog).strip().lower()
#tog = 'n'
if tog =='y':
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

#Option to inspect fitting data versus scattering angle

tog = input('Inspect fitted peak widths (y,[n]):?')
tog = str(tog).strip().lower()
#tog = 'y'
if tog =='y':
    #get 2theta info from detector table
    CreateDetectorTable(InputWorkspace='calFile',DetectorTableWorkspace='calFile_detectorInfo')
    #get 2theta info from detector table for SumNeighbours output
    CreateDetectorTable(InputWorkspace='calFile_8x8',DetectorTableWorkspace='calFile_detectorInfo_8x8')
    # get ttheta values for full detector and SumNeighbours (SN) instrument 
    ws = mtd['calFile_detectorInfo_8x8']
    tableData = ws.toDict()
    ttSN = np.array(tableData['Theta']) #why isn't it called ttheta!!?
    ksortSN = np.argsort(ttSN) #index to sort values in order of increasing ttheta

    ws = mtd['calFile_detectorInfo']
    tableData = ws.toDict()
    tt = np.array(tableData['Theta']) #why isn't it called ttheta!!?
    ksort = np.argsort(tt) #index to sort values in order of increasing ttheta   

    # get chi2 and store as function of ordered tt in workspace
    ws = mtd['_diag_fitparam']
    tableData = ws.toDict()
    chi2 = np.array(tableData['chi2'])# np array of chi2 sorted with increasing ttheta
    CreateWorkspace(OutputWorkspace='Chi2AngularInfo',DataX=ttSN[ksortSN],DataY=chi2[ksortSN])
    #plot chi


    Chi2AngularInfo = ADS.retrieve('Chi2AngularInfo')

    fig, axes = plt.subplots(num='Chi2AngularInfo-11', subplot_kw={'projection': 'mantid'})
    axes.plot(Chi2AngularInfo, color='#1f77b4', label='Chi2AngularInfo: spec 1', wkspIndex=0)
    axes.tick_params(axis='x', which='major', **{'gridOn': False, 'tick1On': True, 'tick2On': False, 'label1On': True, 'label2On': False, 'size': 6, 'tickdir': 'out', 'width': 1})
    axes.tick_params(axis='y', which='major', **{'gridOn': False, 'tick1On': True, 'tick2On': False, 'label1On': True, 'label2On': False, 'size': 6, 'tickdir': 'out', 'width': 1})
    axes.set_title('Chi2AngularInfo')
    axes.set_xlabel('ttheta (deg)')
    axes.set_ylabel('chi2')
    axes.set_ylim([0.0, 5.0])
    #legend = axes.legend(fontsize=8.0).draggable().legend
    plt.show()

    # get widths for all fitted peaks
    ws = mtd['_diag_width']
    tableData = ws.toDict()
    xData = [] 
    yData = []
    fitted_pks = []
    interpx = np.linspace(40.0,120.0,6)
    nSpec=0
    for i in range(nPks):
        
        dictVal = '@' + pks[i].replace(" ","")    
        try:
            tmpY = np.array(tableData[dictVal])
            tmpY = tmpY[ksort]
            fitted_pks.append(float(pks[i])) #only d-spacings for peaks that have been successfully fitted
            print('Got data for peak: ',i,' name: ',dictVal)
            interpy = np.interp(interpx,tt[ksort],tmpY)
            xData.append(np.ndarray.tolist(interpx))
            yData.append(np.ndarray.tolist(interpy))
        except:
            print('Couldn\'t find data for peak', i,dictVal)
        
    print('size of xData=',len(xData))
    nSpec = len(fitted_pks)
    print('number of peaks with fitted data: ',nSpec)
    xData = np.array(xData)
    yData = np.array(yData)
    CreateWorkspace(OutputWorkspace='GaussWidths_vs_tt',DataX=xData,DataY=yData,NSpec=nSpec)
    nfitted,ntt=xData.shape
    # parse peak widths into array as a function of d-spacing for each ttheta bin
    #xData = np.array(pks).astype(np.float) #new array with d-spacings in it
    yData = np.transpose(yData)
    xData = []
    for i in range(ntt):
        xData.append(fitted_pks) #x values are now all d-spacings for peaks at this ttheta
    xData=np.array(xData)    
    CreateWorkspace(OutputWorkspace='GaussWidths_vs_d',DataX=xData,DataY=yData,NSpec=ntt,UnitX='dspacing',\
        YUnitLabel='sigma')    

    GaussWidths_vs_d = ADS.retrieve('GaussWidths_vs_d')

    fig, axes = plt.subplots(num='GaussWidths_vs_d-12', subplot_kw={'projection': 'mantid'})
    axes.plot(GaussWidths_vs_d, color='#1f77b4', label='GaussWidths_vs_d: spec 1', marker='o', wkspIndex=0)
    axes.plot(GaussWidths_vs_d, color='#ff7f0e', label='GaussWidths_vs_d: spec 2', marker='o', wkspIndex=1)
    axes.plot(GaussWidths_vs_d, color='#2ca02c', label='GaussWidths_vs_d: spec 3', marker='o', wkspIndex=2)
    axes.plot(GaussWidths_vs_d, color='#d62728', label='GaussWidths_vs_d: spec 4', marker='o', wkspIndex=3)
    axes.plot(GaussWidths_vs_d, color='#9467bd', label='GaussWidths_vs_d: spec 5', marker='o', wkspIndex=4)
    axes.plot(GaussWidths_vs_d, color='#8c564b', label='GaussWidths_vs_d: spec 6', marker='o', wkspIndex=5)
    axes.tick_params(axis='x', which='major', **{'gridOn': False, 'tick1On': True, 'tick2On': False, 'label1On': True, 'label2On': False, 'size': 6, 'tickdir': 'out', 'width': 1})
    axes.tick_params(axis='y', which='major', **{'gridOn': False, 'tick1On': True, 'tick2On': False, 'label1On': True, 'label2On': False, 'size': 6, 'tickdir': 'out', 'width': 1})
    axes.set_title('GaussWidths_vs_d')
    axes.set_xlabel('d-spacing (Ang)')
    axes.set_ylabel('Gauss sigma')
    #legend = axes.legend(fontsize=8.0).set_draggable().legend

    plt.show()
# Scripting Plots in Mantid:
# https://docs.mantidproject.org/tutorials/python_in_mantid/plotting/02_scripting_plots.html


#check if all is OK then output final calibration...
calibFilename = stateFolder + 'SNAP_calibrate_d%s_'%(run)+ now.strftime("%Y%m%d") + '.h5'
tog = input('Happy with calibration? Save file ([y],n):?')
tog = str(tog).strip().lower()
if tog =='n':
    print('Finished')
    sys.exit()
else:
    print(calibFilename)
    SaveDiffCal(CalibrationWorkspace='_cal',Filename=calibFilename)

#Lastly clean up unused workspaces
##DeleteWorkspace('_cal')
#DeleteWorkspace('_cal_mask')
#DeleteWorkspace('calFile')
#DeleteWorkspace('calFile_8x8')
#DeleteWorkspace('Column')
#myiv = get_instrumentview('calFile')
#myiv.show_view()

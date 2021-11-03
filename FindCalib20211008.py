# Description: least squares optimisation of the SNAP detector angles and L2's by fitting a single reference
# Bragg peak from a known crystalline calibrant
# M. Guthrie 5 OCT 2021
#
# 6 OCT 2021 M. Guthrie modified to work with fits of peak across up to 18 banks


from mantid.simpleapi import *
import matplotlib.pyplot as plt
import sys
sys.path.append('/SNS/SNAP/shared/Malcolm/release/')
from mantidqt.utils.qt.qappthreadcall import QAppThreadCall
from mantid.plots.utility import MantidAxType
from mantid.api import AnalysisDataService as ADS
import numpy as np
import mgutilities20211001 as mg
import importlib
from lmfit import minimize, Parameters, fit_report
import lmfit
print('lmfit version: ',lmfit.__version__)
importlib.reload(mg)
input = QAppThreadCall(mg.workbench_input_fn)
global alliter,allchi2

######################################################################
#Instructions:
#1) update the settings between here and "Don't edit below".
#2) First run the script with RefineAngles = False. This will fix the arcs
#   and refine the L2. when it's done, copy the refined lims below, set
#   RefineAngles to True and then re-run script to refine arcs 
#3) to run script, go to iPython tab in mantidWorkbench and type:
#   %run /SNS/SNAP/shared/Malcolm/release/FindCalib.py
#   then enter
#
#tips: the algorithm needs to the calibration data to have a peak close
#to the specified d-spacing target. Make sure that you choose a peak that
#is visible in all parts of the detector! The diamond 220 works for this,
#
######################################################################

SIPTS = 27111 #ipts for calibration data
CalibRun = 51984 #run number for a powder calibration run

#Set refinement flags
RefineAngles = True #true to refine angles, False to refine L2
detGrouping = 'Column' #e.g. 'Column' or 'bank' use correct capitalisation (I'm not responsible!)
#Initial values
arc1 = -65.641  #degrees
lin1 = 0.04353  #meters
arc2 = 105.3559  #degrees
lin2 = 0.04190  #meters
#Provide correct d-spacing for target peak (MUST BE PRESENT IN ALL BANKS
xTarget = 1.2611585#1.9941
#Useful numbers:
#1.2611585 -220 for diamond from Shikata et al 2018
maxIt = 15 #maximum number of iterations to run

#####################################################################
#Don't edit below here
#####################################################################
instrumentTag = 'SNAP'
extn = 'nxs.h5' #changed from 'nxs' to 'nxs.h5' in later runs
datafolder = 'nexus' #changed from 'data' to 'nexus' at some point
rootDir = '/SNS/'
sharedDir = rootDir+str(instrumentTag)+'/IPTS-'+str(SIPTS)+'/shared/' #equivalent to "shared" directory in SNS file structure  
dataDir = rootDir+str(instrumentTag)+'/IPTS-'+str(SIPTS)+'/'+datafolder+'/' #location of sample nexus files, typically '../data' at SNS

if detGrouping=='Column':
    nHst = 6
elif detGrouping=='bank':
    nHst=18

if mtd.doesExist('%s%s'%(instrumentTag,CalibRun)):
    print('Nexus data already loaded')
else:
    LoadEventNexus(Filename=r'%s%s_%s.%s'%(dataDir,instrumentTag,CalibRun,extn),OutputWorkspace='%s%s'%(instrumentTag,CalibRun),Precount='1', LoadMonitors=False)
    NormaliseByCurrent(InputWorkspace='%s%s'%(instrumentTag,CalibRun),OutputWorkspace='%s%s'%(instrumentTag,CalibRun))  
    CompressEvents(InputWorkspace='%s%s'%(instrumentTag,CalibRun),OutputWorkspace='%s%s'%(instrumentTag,CalibRun))

WSName = '%s%s'%(instrumentTag,CalibRun)
alliter = np.empty(0)
allchi2 = np.empty(0)
#Define fitting parameters
fit_params = Parameters()
fit_params.add('det_arc1', value=arc1)#, min =arc1-arcShift, max = arc1+arcShift)
fit_params.add('det_lin1', value=lin1)#, min =lin1-linShift, max = lin1+linShift)
fit_params.add('det_arc2', value=arc2)#, min =arc2-arcShift, max = arc2+arcShift)
fit_params.add('det_lin2', value=lin2)#, min =lin2-linShift, max = lin2+linShift)
xFirst = xTarget -0.08#minimum of fitting range
xLast = xTarget +0.08 #maximum of fitting range


if plt.fignum_exists('chi2Plot'): 
    plt.close('chi2Plot') #get rid of window if it already exists
if plt.fignum_exists('Calib'): 
    plt.close('Calib') #get rid of window if it already exists

      # Set up column grouping workspace if it doesn't yet exist
CreateGroupingWorkspace(InstrumentFilename='/SNS/SNAP/shared/Malcolm/dataFiles/SNAP_Definition.xml', \
GroupDetectorsBy=detGrouping, OutputWorkspace='SNAPGpWS')

if RefineAngles:
    fit_params['det_arc1'].vary = True
    fit_params['det_lin1'].vary = False
    fit_params['det_arc2'].vary = True
    fit_params['det_lin2'].vary = False
else:
    fit_params['det_arc1'].vary = False
    fit_params['det_lin1'].vary = True
    fit_params['det_arc2'].vary = False
    fit_params['det_lin2'].vary = True


def do_CHI2_plot(w):

    chi2Plot = ADS.retrieve('chi2valsCurrentIter')
    
    fig, axes = plt.subplots(num='chi2Plot', subplot_kw={'projection': 'mantid'})
    axes.clear()
    axes.plot(chi2Plot, color='#1f77b4', label='running Chi2', wkspIndex=0)

    plt.show()


def iteration_output(pars,iter,resid,*args,**kws):
    global alliter,allchi2
    chi2 = np.sum(np.square(resid))/resid.size
    alliter = np.append(alliter,iter+2)
    allchi2 = np.append(allchi2,chi2)
    print('Iteration:',iter+2,'chi2:',chi2)
    chi2valsCurrentIter = CreateWorkspace(DataX=alliter,DataY=allchi2,NSpec=1,UnitX='Iteration',YUnitLabel='Chi2')
    chiplot2 = do_CHI2_plot([chi2valsCurrentIter])
    return

def xDistrib(fit_params,WSName,xFirst,xLast,xTarget,nHst):

    #Defined function to minimise and returns a residual


    arc1 = str(fit_params['det_arc1'].value)
    lin1 = str(fit_params['det_lin1'].value)
    arc2 = str(fit_params['det_arc2'].value)
    lin2 = str(fit_params['det_lin2'].value)

    #GpWSName = 
    AddSampleLog(Workspace=WSName,LogName='det_arc1',LogText=arc1,LogType='Number Series')
    AddSampleLog(Workspace=WSName,LogName='det_lin1',LogText=lin1,LogType='Number Series')
    AddSampleLog(Workspace=WSName,LogName='det_arc2',LogText=arc2,LogType='Number Series')
    AddSampleLog(Workspace=WSName,LogName='det_lin2',LogText=lin2,LogType='Number Series')
    LoadInstrument(Workspace=WSName,MonitorList='-1,1179648', RewriteSpectraMap='False',InstrumentName='SNAP')
    ConvertUnits(InputWorkspace=WSName,OutputWorkspace='%s_d'%(WSName),Target='dSpacing')
    Rebin(InputWorkspace='%s_d'%(WSName),OutputWorkspace='%s_d'%(WSName),Params='0.8,-0.001,2.5')
    DiffractionFocussing(InputWorkspace='%s_d'%(WSName), OutputWorkspace='%s_d'%(WSName)+str(nHst).zfill(2), GroupingWorkspace='SNAPGpWS')

    sumChi2 = 0
   
    resid = np.zeros(nHst)
    for i in range(nHst):
        #Function="name=Gaussian,Height=7596.55,PeakCentre=1.2616,Sigma=0.00512778;name=Polynomial,n=0,A0=0"
        Function="name=Gaussian,Height=7596.55,PeakCentre="+str(xTarget)+",Sigma=0.015;name=Polynomial,n=0,A0=0"
        InputWorkspace= '%s_d'%(WSName)+str(nHst).zfill(2)
        WorkspaceIndex=i
        Output="res"
        StartX=xFirst
        EndX=xLast
        Normalise=True
        OutputCompositeMembers=True
        fit_output = Fit(EndX=EndX,Function=Function,InputWorkspace=InputWorkspace,Normalise=Normalise,Output=Output,OutputCompositeMembers=OutputCompositeMembers,StartX=StartX,WorkspaceIndex=WorkspaceIndex)
        paramTable = fit_output.OutputParameters
        xFit = paramTable.column(1)[1]
        sig = paramTable.column(1)[2]#Gauss peak width, normalise to this

        #print('d-spacing is:',xFit,sig)
        resid[i]= (xFit-xTarget)/sig
        #print('sumChi2:',sumChi2)
    #print(resid)    
    return resid

#res = minimize(xDistrib, fit_params,args=(WSName,xFirst,xLast,xTarget),iter_cb=iteration_output,epsfcn=1e-6)
if RefineAngles:
    res = minimize(xDistrib, fit_params,args=(WSName,xFirst,xLast,xTarget,nHst),iter_cb=iteration_output,maxfev=maxIt-2)
else:
    res = minimize(xDistrib, fit_params,args=(WSName,xFirst,xLast,xTarget,nHst),iter_cb=iteration_output,epsfcn=1e-6,maxfev=maxIt-2)

print(fit_report(res))


FindCalib_BestFit = ADS.retrieve('%s_d'%(WSName)+str(nHst).zfill(2))
if plt.fignum_exists('Calib'): 
    plt.close('Calib') #get rid of window if it already exists

colours = ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f','#bcbd22',\
    '#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f','#bcbd22']


fig, axes = plt.subplots(num='Calib', subplot_kw={'projection': 'mantid'})
for i in range(nHst):
    axes.plot(FindCalib_BestFit, color=colours[i], label='spec '+str(i+1), wkspIndex=i)
PlotYLims = axes.get_ylim()
#show target d-spacing
axes.plot([xTarget,xTarget],[PlotYLims[0],PlotYLims[1]],color='#000000',linestyle='--', linewidth=0.5)
axes.tick_params(axis='x', which='major', **{'gridOn': False, 'tick1On': True, 'tick2On': False, 'label1On': True, 'label2On': False, 'size': 6, 'tickdir': 'out', 'width': 1})
axes.tick_params(axis='y', which='major', **{'gridOn': False, 'tick1On': True, 'tick2On': False, 'label1On': True, 'label2On': False, 'size': 6, 'tickdir': 'out', 'width': 1})
axes.set_title('Best fit')
axes.set_xlabel('d-Spacing ($\\AA$)')
axes.set_ylabel('Counts ($\\AA$)$^{-1}$')
axes.set_xlim([xFirst, xLast])
#axes.set_ylim([-332330.0, 48774000.0])
legend = axes.legend(fontsize=8.0).draggable().legend
plt.show()

ISAWFileName=sharedDir + '%s%s'%(instrumentTag,CalibRun) +'.detcal'
tog = input('Dialogue','Do you want to save and ISAW detcal file? (y/n)','str')
if tog.lower() == 'y':
    SaveIsawDetCal(InputWorkspace=WSName,Filename=ISAWFileName)
    print('wrote:',ISAWFileName)
from mantid.kernel import *
from mantid.api import *
from mantid.simpleapi import *
import sys
sys.path.append('/SNS/users/66j/SNAPSoftware/SANDPitSNAP/')
import mgutilities as mg
import importlib
importlib.reload(mg)

run = sys.argv[1]

IPTSLoc = GetIPTS(run,'SNAP')

try:
    NexusLoc = GetIPTS(run,'SNAP') + '/nexus/SNAP_' +str(run) + '.nxs.h5'
except:
    print('ERROR: Nexus file not found')

#Load Nexus Log
CreateWorkspace(OutputWorkspace='Dummy2',DataX=0,DataY=0)
LoadNexusLogs(Workspace='Dummy2',Filename=NexusLoc)
ws = mtd['Dummy2']
logRun = ws.getRun()

# get log data from nexus file
print('\n/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_')
print('Log Values:')
fail = False
try:
    det_arc1 = logRun.getLogData('det_arc1').value[0]
    print('det_arc1 is:',det_arc1)
except:
    print('ERROR: Nexus file doesn\'t contain value for det_arc1')
    fail = True

try:    
    det_arc2 = logRun.getLogData('det_arc2').value[0]
    print('det_arc2 is:',det_arc2)
except:
    print('ERROR: Nexus file doesn\'t contain value for det_arc2')
    fail = True

try:
    wav = logRun.getLogData('BL3:Chop:Skf1:WavelengthUserReq').value[0]
    print('wav Skf1 wavelengthUserReq is:',wav, 'Ang.')
except:
    print('ERROR: Nexus file doesn\'t contain value for central wavelength')
    fail = True

try:
    freq = logRun.getLogData('BL3:Det:TH:BL:Frequency').value[0]
    print('frequency setting is:',freq, 'Hz')
except:
    print('ERROR: Nexus file doesn\'t contain value for central wavelength')
    fail = True

try:
    GuideIn = logRun.getLogData('BL3:Mot:OpticsPos:Pos').value[0]
    print('guide status is:',GuideIn)
except:
    print('ERROR: Nexus file doesn\'t contain guide status')
    fail = True
print('/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_')


if not fail:
    stateID = mg.checkSNAPState([det_arc1,det_arc2,wav,freq,0.0],[GuideIn,0])
else:
    print('Insufficient log data, can\'t determine state')
DeleteWorkspace(Workspace='Dummy2')
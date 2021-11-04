from mantid.kernel import *
from mantid.api import *
from mantid.simpleapi import *
import sys
sys.path.append('/SNS/users/66j/SNAPSoftware/SANDPitSNAP/')
import mgutilities as mg
import importlib
importlib.reload(mg)

class SNAPConfigurator(PythonAlgorithm):

    def PyInit(self):
        # Declare properties
        #validator = IntArrayBoundedValidator(lower=0)
        self.declareProperty("RunNumber", 12345, direction=Direction.Input)
        self.declareProperty("Manually choose config file", False,
                             "manually select DRC file")
        self.declareProperty(FileProperty(name="DRCFile", defaultValue="/SNS/SNAP/shared/Calibration/",
                                    extensions=['.json'],
                                    direction=Direction.Input,
                                    action=FileAction.OptionalLoad),
                    doc="The calibration file to convert to d_spacing.")

    def PyExec(self):
        # Run the algorithm
        if self.getProperty("Manually choose config file").value:
            print('I will read from file')
        else:
            run = self.getProperty("RunNumber").value
            #print('I will get info from a run')
            print('Run is:',run)
            IPTSLoc = GetIPTS(run,'SNAP')
            NexusLoc = GetIPTS(run,'SNAP') + '/nexus/SNAP_' +str(run) + '.nxs.h5'
            CreateWorkspace(OutputWorkspace='Dummy2',DataX=0,DataY=0)
            LoadNexusLogs(Workspace='Dummy2',Filename=NexusLoc)
            ws = mtd['Dummy2']
            logRun = ws.getRun()
            det_arc1 = logRun.getLogData('det_arc1').value[0]
            print('det_arc1 is:',det_arc1)
            det_arc2 = logRun.getLogData('det_arc2').value[0]
            print('det_arc2 is:',det_arc2)
            wav = logRun.getLogData('BL3:Chop:Skf1:WavelengthUserReq').value[0]
            print('wav Skf1 wavelengthUserReq is:',wav)
            GuideIn = logRun.getLogData('BL3:Mot:OpticsPos:Pos').value[0]
            #
            # ADD CODE HERE TO CHECK IF GUIDE IS NEW OR OLD
            #
            print('guide status is:',GuideIn)
            DeleteWorkspace(Workspace='Dummy2')
            stateID = mg.getSNAPStateID(det_arc1,det_arc2,wav,[1.0,1.0,0.1],GuideIn,0)
            if stateID=='0000' #no match was found
                mg.createSNAPStateID('short name',det_arc1,det_arc2,wav,GuideIn)


            

# Register algorithm with Mantid
AlgorithmFactory.subscribe(SNAPConfigurator)

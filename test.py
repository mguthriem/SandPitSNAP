#test 
import sys
sys.path.append('/SNS/users/66j/SNAPSoftware/SANDPitSNAP/')
import mgutilities as mg
import importlib
importlib.reload(mg)

#state = mg.getSNAPStateID(-65.371,104.944,2.1,[1.0,1.0,0.1],1,0)
#print(state)
stateID = mg.checkSNAPState([-65.,105.0,2.1,0.0],[2,0])

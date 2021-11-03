import glob
import os
import os.path, time
from datetime import datetime

pattern = '/SNS/SNAP/shared/Calibration/**/1010/config1010*.json'
FindMostRecent = False

refDate = datetime.now().timestamp()

for fname in glob.glob(pattern, recursive=True):
    ShortestTimeDifference = 10000000000 # a large number of seconds
    if os.path.isfile(fname):
        #rint(fname)
        #print("Created: %s" % time.ctime(os.path.getctime(fname)))
        #print('epoch:',os.path.getctime(fname))
        #print('refdate epoch:',refDate)
        delta = refDate - os.path.getctime(fname)
        #print('difference:',delta)
        if delta <= ShortestTimeDifference:
            MostRecentFile = fname
            ShortestTimeDifference = delta
if ShortestTimeDifference == 10000000000:
    print('no matching file found')
else:
    print('Most recent matching file:',fname)
    print('Created: %s'% time.ctime(os.path.getctime(fname)))

        #print(refDate-vvos.path.getctime(fname))
        #timestr = SetUpDate.strftime("_%d-%b-%Y-%H%M%S")
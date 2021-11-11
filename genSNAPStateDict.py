import json

SNAPStateDict = {\
"westArc":(-50.0,-65.3,-76.1,-90.0,-105.0,-115.0),\
"eastArc":[-50.2,-66.8,-75.7,-90.0,104.9,115.0],\
"waveSet":[1.5,2.1,2.4],\
"tolerance":[1.0,1.0,0.1,0.1],\
"spareFloat":[0.0],\
"floatParameterOrder":['westArc','eastArc','waveSet','spareFloat'],\
"guideStatus":[0,1,2],\
"spareInt":[0],\
"intParameterOrder":['guideStatus','spareInt'],
}
fname = '/SNS/SNAP/shared/Calibration/SNAPStateDict_20211105.json'
with open(fname, "w") as outfile:
    json.dump(SNAPStateDict, outfile)

import json

SNAPStateDict = {\
"westArc":[-50.0,-65.0,-90.0,-105.0],\
"eastArc":[-50.0,-65.0,-90.0,-105.0],\
"waveSet":[1.5,2.1,2.4],\
"tolerance":[1.0,1.0,0.1],\
"spareFloat":[0.0],\
"floatParameterOrder":['westArc','eastArc','waveSet','spareFloat'],\
"guideStatus":[0,1,2],\
"spareInt":[0]\
"intParameterOrder":['guideStatus','spareInt'],\
}
fname = '/SNS/SNAP/shared/Calibration/SNAPStateDict_20211104.json'
with open(fname, "w") as outfile:
    json.dump(SNAPStateDict, outfile)

    # with open(ConfigFileName, "w") as outfile:
    #     json.dump(configDict, outfile)
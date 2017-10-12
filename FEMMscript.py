#relies on sign function in Lua code

def nodesGen(map): #produces addnode Lua code from list of points
	nodeCode = []
	counter = 0
	for node in map:
		nodeCode.append("ei_addnode(" + node[0] + "," + node[1] + ") -- " + str(counter))
		counter += 1
	return nodeCode
	
def segsGen(map, arcs): #produces addsegment Lua code from points; ignores connections that should be arcs
	segsCode = []
	counter = 0
	for i in range(len(map)):
		if (str(counter)+"-"+str((counter+1)%(len(map))) not in arcs.keys() and str((counter+1)%(len(map)))+"-"+str(counter) not in arcs.keys()):
			segsCode.append("ei_addsegment(" + map[i][0] + "," + map[i][1] + "," + map[(i+1)%len(map)][0] + "," + map[(i+1)%len(map)][1] + ") -- " + str(counter) + "-" + str((counter+1)%(len(map))))
		counter += 1
	return segsCode
	
def arcsGen(map, arcs): #produces addarcsegment Lua code from points and list of nodes connected by arcs
	arcsCode = []
	counter = 0
	for i in range(len(map)):
		if (str(counter)+"-"+str((counter+1)%(len(map))) in arcs.keys()):
			arcsCode.append("ei_addarc(" + map[i][0] + "," + map[i][1] + "," + map[(i+1)%len(map)][0] + "," + map[(i+1)%len(map)][1] + "," + arcs[str(counter)+"-"+str(counter+1)][0] + "*180/pi," +  arcs[str(counter)+"-"+str(counter+1)][1] + ") -- " + str(counter) + "-" + str((counter+1)%(len(map))))
		elif (str((counter+1)%(len(map)))+"-"+str(counter) in arcs.keys()):
			arcsCode.append("ei_addarc(" + map[(i+1)%len(map)][0] + "," + map[(i+1)%len(map)][1] + "," + map[i][0] + "," + map[i][1] + "," + arcs[str(counter+1)+"-"+str(counter)][0] + "*180/pi," +  arcs[str(counter+1)+"-"+str(counter)][1] + ") -- " + str((counter+1)%(len(map))) + "-" + str(counter))
		counter += 1
	return arcsCode

def selectSegs(map, arcs): #produces selectsegment Lua code, ignoring arc connections
	selectSegsCode = []
	counter = 0
	for i in range(len(map)):
		if (str(counter)+"-"+str((counter+1)%(len(map))) not in arcs.keys() and str((counter+1)%(len(map)))+"-"+str(counter) not in arcs.keys() and str(counter)+"-"+str((counter+1)%len(map)) not in presets["selectexcept"].split(",")):
			x0, y0 = map[i][0], map[i][1]
			x1, y1 = map[(i+1)%len(map)][0], map[(i+1)%len(map)][1]
			signx = "sign(" + x1 + "-" + x0 + ")"
			signy = "sign(" + y1 + "-" + y0 + ")"
			selectSegsCode.append("ei_selectsegment(" + x0 + " + 0.1*" + signx + "," + y0 + " + 0.1*" + signy + ") -- " + str(counter) + "-" + str((counter+1)%(len(map)))) 
		counter += 1
	return selectSegsCode

def selectArcs(map, arcs): #produces selectarcsegment Lua code, using list of arc connections
	selectArcsCode = []
	counter = 0
	for i in range(len(map)):
		if (str(counter)+"-"+str((counter+1)%(len(map))) in arcs.keys() or str((counter+1)%(len(map)))+"-"+str(counter) in arcs.keys()):
			x0, y0 = map[i][0], map[i][1]
			x1, y1 = map[(i+1)%len(map)][0], map[(i+1)%len(map)][1]
			signx = "sign(" + x1 + "-" + x0 + ")"
			signy = "sign(" + y1 + "-" + y0 + ")"
			selectArcsCode.append("ei_selectarcsegment(" + x0 + " + 0.1*" + signx + "," + y0 + "+ 0.1*" + signy + ") -- " + str(counter) + "-" + str((counter+1)%(len(map))))
		counter += 1
	return selectArcsCode

def segmentProps(map, presets): #adds setsegmentprop code using presets
	segPropCode = []
	segPropCode.append("ei_setsegmentprop(" + presets["segmentpropname"] + "," + presets["segmentelementsize"] + "," + presets["segmentautomesh"] + "," + presets["segmenthide"] + "," + presets["segmentgroup"] + "," + presets["segmentinconductor"]+ ")")
	segPropCode.append("ei_clearselected()")
	return segPropCode
	
def arcProps(map, presets): #adds setarcsegmentprop code using presets
	arcPropCode = []
	arcPropCode.append("ei_setarcsegmentprop(" + presets["arcmaxsegdeg"] + "," + presets["arcpropname"] + "," + presets["archide"] + "," + presets["arcgroup"] + "," + presets["arcinconductor"]  + ")")
	arcPropCode.append("ei_clearselected()")
	return arcPropCode

def blockProps(presets): #adds block properties using presets
	blockPropCode = []
	x, y = split(presets["blockxy"])
	blockPropCode.append("ei_addblocklabel(" + x + "," + y + ")")
	blockPropCode.append("ei_selectlabel(" + x + "," + y + ")")
	blockPropCode.append("ei_setblockprop(" + presets["blockname"] + "," + presets["blockautomesh"] + "," + presets["meshsize"] + "," + presets["blockgroup"] + ")")
	blockPropCode.append("ei_clearselected()")
	return blockPropCode

def split(str): #splits (x,y) string it [x,y] list
	counter = 0
	for i in range(len(str)):
		if str[i] == "(":
			counter += 1
		elif str[i] == ")":
			counter -= 1
		elif (counter == 0 and str[i] == ","):
			return str[:i], str[i+1:]

def process(file):
	file = file.read().split("nodes:")
	presetsList = file[0].split("\n")
	presetsList.pop(len(presetsList)-1) #remove orphaned newline
	arcsList = presetsList.pop(len(presetsList)-1)[5:].split(",") #[5:] removes "arcs:"
	map = file[1].split("\n")[1:] #remove orphaned newline
	presets = {}
	arcs = {}
	for preset in presetsList:
		p = preset.split(":")
		presets[p[0]] = p[1]
	for i in range(0, len(arcsList)-1, 3):
		arcs[arcsList[i]] = [arcsList[i+1], arcsList[i+2]]
	for i in range(len(map)):
		map[i] = split(map[i]) #splits (x,y) string into [x, y] list
	return presets, arcs, map

#files should be placed in the same folder as this file and their paths within this folder should be put in the array files as strings
files = ["vacuumChamber.txt", "topElectrodeHalf.txt", "topElectrodeFull.txt", "middleElectrode.txt", "bottomElectrodeHalf.txt", "bottomElectrodeFull.txt", "mainArea.txt"]
for file in files:
	print "-- " + file.split(".")[0] #makes Lua comment labelling section of code
	print "if " + file.split(".")[0] + " == 1 then" #relies on file name being variable to toggle this code
	nodeFile = open(file, 'r')
	presets, arcs, map = process(nodeFile)
	nodeCode = nodesGen(map)
	
	if nodeCode: #if list isn't empty
		for node in nodeCode:
			print node
		print "\n"	
		
		segsCode = segsGen(map, arcs)
		if segsCode: #if list isn't empty
			for seg in segsCode:
				print seg
			print "\n"
			selectSegsCode = selectSegs(map, arcs)
			for seg in selectSegsCode:
				print seg
			print "\n"
			segmentPropCode = segmentProps(map, presets)
			for line in segmentPropCode:
				print line
			print "\n"
		arcsCode = arcsGen(map, arcs)
		if arcsCode: #if list isn't empty
			for arc in arcsCode:
				print arc
			print "\n"
			selectArcsCode = selectArcs(map, arcs)
			for arc in selectArcsCode:
				print arc
			print "\n"
			arcPropCode = arcProps(map, presets)
			for line in arcPropCode:
				print line
			print "\n"
	blockPropCode = blockProps(presets)
	for line in blockPropCode:
		print line
	print "end\n"
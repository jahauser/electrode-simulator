-- EDM cell electrode modelling, written by Jake Hauser during his 2016 summer fellowship
--############################### Modelling Functions ###############################

 -- Returns 1 for positive x, -1 for negative x, and 0 for zero
function sign(x)
	if x > 0 then
		return 1
	elseif x < 0 then
		return -1
	end
	return 0
end

-- Generates x-coordinate values (corresponding to r-values in FEMM) for a variety of parametric electrode profile equations
-- t is the parameter for these equations
-- uses global values: curveType, scale, pi, k0, k1, k2
function x(t)
	if curveType == "ernstB" then
		return scale*(t-.25*(1-sqrt(1-k0^2))*sinh(2*t)) -- Ernst B equation; source: UNIFORM-FIELD ELECTRODES WITH MINIMUM WIDTH - Gerard J. ERNST, 1984
	elseif curveType == "ernstC" then
		return scale*(t-k1*sinh(2*t)) -- Ernst C equation; source: UNIFORM-FIELD ELECTRODES WITH MINIMUM WIDTH - Gerard J. ERNST, 1984
	elseif curveType == "rogowski" then
		return t -- The Rogowski profile is not usually a parametric equation so x = t; source: https://de.wikipedia.org/wiki/Rogowski-Profil
	elseif curveType == "flat" then
		return t -- Horizontal line for flat electrode
	end
end

-- Generates y-coordinate values (corresponding to z-values in FEMM) for a variety of electrode profile equations
-- t is the parameter for these equations
-- uses global values: curveType, scale, pi, k0, k1, k2
function y(t)
	if curveType == "ernstB" then
		return scale*(pi/2 + k0*cosh(t)) -- Ernst B equation; source: UNIFORM-FIELD ELECTRODES WITH MINIMUM WIDTH - Gerard J. ERNST, 1984
	elseif curveType == "ernstC" then
		return scale*(pi/2 + k0*cosh(t) - k2*cosh(3*t)) -- Ernst C equation; source: UNIFORM-FIELD ELECTRODES WITH MINIMUM WIDTH - Gerard J. ERNST, 1984
	elseif curveType == "rogowski" then
		return electrodeSeparation/pi * (pi/2 + exp(t*pi/electrodeSeparation)) -- Rogowski equation; source: https://de.wikipedia.org/wiki/Rogowski-Profil
	elseif curveType == "flat" then
		return 0 -- Horizontal line for flat electrode
	end
end

-- Approximates the angle of a curve's slope at point x(t), y(t) using the slope angle between the nearest points on either side (points at t+step and t-step)
-- Returns angle not slope to avoid attempting atan(nil)
-- uses functions: x(), y()
function angle(t, step)
	dy = y(t+step)-y(t-step)
	dx = x(t+step)-x(t-step)
	if dx <= 0 then -- When dx = 0, dy/dx is undefined; when dx < 0 the curve has gone too far (it should only go to vertical) -- in both cases, nil is returned
		return nil
	end
	return atan(dy/dx)
end

-- Produces a 2D array of points and slope angles in a curve, called a map
-- Returns this map along with the radius (difference between final and inital x), thickness (difference between final and inital y), slope angle at the final point, and number of total points
-- Can generate a map either over a desired t-range (indicated by tf = a value) or until the curve becomes vertical (indicated by tf = nil)
-- step is the amount t increments by; vshift is a translation applied to y-coordinates; vflip is a multiplier for y-coordinates
-- uses functions: angle(), x(), y()
function mapGen(step, vshift, vflip, t0, tf)
	map = {} -- Empty array
	i = 1 -- i steps through the integers to indicate position in the array
	t = t0 -- t steps by a given parameter for use in x(t) and y(t) functions
	theta = angle(t, step) -- sets a starting theta so "while theta ~= nil do" takes place at least once
	
	-- Values used for radius and thickness calculations, not put in map
	x0 = x(t)
	y0 = y(t)
	
	if tf == nil then -- Generates points until curve becomes vertical
		while theta ~= nil do
			t = t0 + (i-1)*step -- (i-1) not i since i starts at 1 but t should start at 0
			map[i] = {x(t), y(t), angle(t, step)} -- adds point
			
			theta = angle(t0 + i*step, step) -- if angle at next t-value == nil (curve ending is vertical) then the loop will break
			i = i + 1 -- increments i
		end
		
		-- Generates final point separately to get final values for radius and thickness calculations
		t = (i-1)*step
		map[i] = {x(t), y(t), pi/2}	
		
		xf = x(t)
		yf = y(t)
		thetaf = pi/2 -- Even if the final angle isn't exactly pi/2 radians, it should be very close and assuming such improves accuracy of connection between two curves to make an electrode
	else -- Generates points over certain t-range
		for t=t0, tf-step, step do
			map[i] = {x(t), y(t), angle(t, step)} -- adds point
			i = i + 1 -- increments i
		end
		
		-- Generates final point separately to get final values for radius and thickness calculations
		t = tf
		map[i] = {x(t), y(t), angle(t, step)}
		
		xf = x(t)
		yf = y(t)
		thetaf = angle(t, step) -- The final angle is unlikely to be pi/2 radians in this case (curve ending isn't vertical)
	end
	
	radius = xf - x0
	thickness = yf - y0
	
	-- Adjusting coordinates for proper drawing
	for j=1, i do
		map[j][1] = map[j][1] - x0 -- In cases where the initial x-value isn't zero (eg. Rogowski profiles) this shifts the curve to start at r = 0
		map[j][2] = (map[j][2] - y0 - thickness) * vflip + vshift -- Inside the parentheses, the curve is shifted so that it starts at y = -thickness and ends at y = 0; afterwards it is flipped about the horizontal if necessary and shifted to make room for a rounding connection or to be a top/bottom electrode
	end
	
	imax = i -- Highest i, corresponds to number of points
	
	return map, radius, thickness, thetaf, imax
end

-- Uses mapGen() to actually draw points, arcs, and segments in FEMM
-- Draws the first node outside a loop. In the first iteration of the loop it draws the second and the connection between the first and second. In future iterations it draws the next point and the connection to the previous one until it has drawn all points.
-- Currently unimplemented is the ability to pause drawing over an x-interval to make gaps for elements like gas intake or O-ring trenches
-- step, vshift, vflip, t0, and tf are all communicated to mapGen(); maxseg, meshsize, and inconductor are boundary FEMM settings for arcs and segments
-- use globals: pi; functions: mapGen(), sign()
function drawCurve(step, maxseg, meshsize, inconductor, vshift, vflip, t0, tf)
	t0 = t0 or 0 -- if no t0 is provided it is set to 0
	tf = tf or nil -- if no tf is provided it is set to nil
	map, radius, thickness, thetaf, imax = mapGen(step, vshift, vflip, t0, tf)
	ei_addnode(map[1][1], map[1][2]) -- draws first node

	for i = 2, imax, 1 do
		ei_addnode(map[i][1], map[i][2]) -- draws next node 
		
		-- determines the difference between the new and previous points' slope angles (in degrees, which FEMM takes) 
		-- if this difference is less than one we have to draw a segment not an arc since FEMM refuses to make arcs with an angle less than one 
		dtheta = (map[i][3]-map[i-1][3])/pi*180
		
		-- distinction made between vflip*dtheta >= 1 and <= -1 because FEMM accepts only angles > 1 degree and builds curves counter-clockwise, so the order of points is flipped and the angle multiplied by -1 for negative vflip*dtheta
		if vflip*dtheta >= 1 then
			ei_addarc(map[i-1][1], map[i-1][2], map[i][1], map[i][2], vflip*dtheta, maxseg)
			ei_selectarcsegment(map[i-1][1] + 0.1*sign(map[i][1]-map[i-1][1]), map[i-1][2] + 0.1*sign(map[i][2]-map[i-1][2]))
			ei_setarcsegmentprop(meshsize, "", 0, 1, inconductor)
		elseif vflip*dtheta <= -1 then
			ei_addarc(map[i][1], map[i][2], map[i-1][1], map[i-1][2], -vflip*dtheta, maxseg)
			ei_selectarcsegment(map[i-1][1] + 0.1*sign(map[i][1]-map[i-1][1]), map[i-1][2] + 0.1*sign(map[i][2]-map[i-1][2]))
			ei_setarcsegmentprop(meshsize, "", 0, 1, inconductor)
		else
			ei_addsegment(map[i-1][1], map[i-1][2], map[i][1], map[i][2])
			ei_selectsegment(map[i-1][1] + 0.1*sign(map[i][1]-map[i-1][1]), map[i-1][2] + 0.1*sign(map[i][2]-map[i-1][2]))
			ei_setsegmentprop("", meshsize, 0, 0, 1, inconductor)
		end
	end
end

--###################################################################################
--################################ Model Build Info #################################
--###################################################################################

--################################## Key Globals ####################################
pi = 3.14159265359

-- "electrode" and "curveType" are sort of confusing variables - sorry. To properly compare different electrode shapes, it's useful for them to have similar thicknesses and radii.
-- Therefore, while "electrode" stays the same and defines the desired electrode shape, "curveType" can change to calculate one profile's thickness and radius in order to properly shift the desired electrode curve.
-- "flat", "rogowski", "ernstB", or "ernstC"
electrode = "ernstC" -- Electrode shape; unchanged through document. Change this to change what electrode is drawn.
curveType = "ernstC" -- Curve style; changed in document. The value listed here is used to calculate a certain ideal radius and thickness. Later, "curveType" is changed to match "electrode"

electrodeSeparation = 180 -- Distance between near edges of a pair of electrodes at r = 0
t0, tf = 0, nil

-- k values, used for Ernst profiles and other electrode shapes made with similar dimensions to Ernst profiles
k0ByRadius = 0 -- 0: user selects desired k0; 1: k0 approximated from desired electrode radius
if k0ByRadius == 0 then
	k0 = 0.03 -- Desired k0 value
elseif k0ByRadius == 1 then
	radius =  225 -- Desired electrode radius
	k0 = 2*sqrt(1/cosh(2*pi*radius/electrodeSeparation+1)) -- Approximation for k0 from desired radius; 0.5% error between 200 and 300 mm, error increases as radius decreases
end
k1 = 1/8*k0*k0
k2 = 1/90*k0*k0*k0

scale = electrodeSeparation/2/(pi/2+k0) -- Multiplier for Ernst profiles -- makes them the proper size for a given separation

-- Component Modulation
topElectrode = 1
middleElectrode = 1
bottomElectrode = 1
vacuumChamber = 1
mainArea = 1
upperBlock = 1
lowerBlock = 1
fullTop = 1
fullBottom = 1

--################################### FEMM Setup ####################################
-- New Window
clearconsole() -- clears FEMM's built-in Lua console
newdocument(1) -- 1 creates an electrostatics problem

-- Voltages
highVoltage = 1000 * electrodeSeparation -- high voltage V; defined for 1 kV/mm between each electrode pair
lowVoltage = 0 -- low voltage V

-- Mesh Size Constraints
meshSizeCoarser = 100
meshSizeCoarse = 10
meshSizeMedium = 1
meshSizeFine = 0.1
meshSizeUltraFine = 0.01

-- Environment Setup
ei_probdef("millimeters","axi",1e-8,1,10) -- Creates an axisymmetric problem in millimeters with a precision of 1e-8 and a minimum mesh angle of 10 degrees
ei_addconductorprop("elecCond",highVoltage,0,1) -- Adds a conductor property named "elecCond" defined by a set voltage, highVoltage
ei_addconductorprop("zeroCond",lowVoltage,0,1) -- Adds a conductor property named "zeroCond" defined by a set voltage, lowVoltage
ei_addmaterial("Glass",4.7,4.7,0) -- Adds a material named "Glass" with relative permittivities of 4.7 in the x/r- and y/z-directions (value take from WikiPedia)
ei_addmaterial("Vacuum",1,1,0) -- Adds a material named "Vacuum" with relative permittivities of 1 in the x/r- and y/z-directions

-- Ideal radius and thickness values off which to model electrodes
-- Thickness always describes a half-electrode as they're drawn in halves
step = 0.025 -- higher values are usually required for Ernst profiles and lower for Rogowski and flat
nullMap, idealRadius, idealThickness, nullTheta = mapGen(step, 0, 1, 0, nil) -- Only the radius and thickness of this curve matter

-- Electrode Variables
innerRadius = idealRadius - idealThickness
blockRadius = innerRadius
if electrode == "ernstC" then
	curveType = "ernstC"
	step = 0.025
	t0 = 0
	tf = nil	
	nullMap, curveRadius, curveThickness, thetaf = mapGen(step, 0, 1, t0, tf)
	
	electrodeThickness = curveThickness
	roundedThickness = 0
elseif electrode == "ernstB" then
	curveType = "ernstB"
	step = 0.025
	t0 = 0
	tf = nil
	nullMap, curveRadius, curveThickness, thetaf = mapGen(step, 0, 1, t0, tf)
	
	electrodeThickness = curveThickness
	roundedThickness = 0
elseif electrode == "rogowski" then
	curveType = "rogowski"

	step = 1

	t0 = -250
	tf = -50
	nullMap, curveRadius, curveThickness, thetaf = mapGen(step, 0, 1, t0, tf)
	
	electrodeThickness = idealThickness
	roundedThickness = idealThickness-curveThickness
	roundedWidth = roundedThickness*(1/sin(pi/2-thetaf)-1/tan(pi/2-thetaf))

	
elseif electrode == "flat" then
	curveType = "flat"
	step = 1
	t0 = 0
	tf = innerRadius
	nullMap, curveRadius, curveThickness, thetaf = mapGen(step, 0, 1, t0, tf)
	
	electrodeThickness = idealThickness
	roundedThickness = idealThickness
	roundedWidth = roundedThickness
end

-- Vacuum Chamber Variables
vacuumScale = 2
vacuumHeight = 2*(2*electrodeThickness + electrodeSeparation) * vacuumScale
vacuumRadius = (innerRadius+electrodeThickness)*2 * vacuumScale
vacuumMesh = meshSizeCoarse -- Defined up here for easy manipulation of vacuumChamber size and detail

-- Main Problem Area Variables; main area is used for increased precision when the vacuum chamber is large
mainHeight = 2*(2*electrodeThickness + electrodeSeparation)
mainRadius = (innerRadius+electrodeThickness)*2
mainAreaMesh = meshSizeCoarse

--#################################### Geometry #####################################

-- vacuumChamber; built first to keep boundary conditions along whole segment, even when it's split by future nodes
-- built here from separate python script and .txt document
-- nodes are numbered to make clear which segments connect which points
if vacuumChamber == 1 then
	ei_addnode((0),(vacuumHeight/2)) -- 0
	ei_addnode((vacuumRadius),(vacuumHeight/2)) -- 1
	ei_addnode((vacuumRadius),(-vacuumHeight/2)) -- 2
	ei_addnode((0),(-vacuumHeight/2)) -- 3

	ei_addsegment((0),(vacuumHeight/2),(vacuumRadius),(vacuumHeight/2)) -- 0-1
	ei_addsegment((vacuumRadius),(vacuumHeight/2),(vacuumRadius),(-vacuumHeight/2)) -- 1-2
	ei_addsegment((vacuumRadius),(-vacuumHeight/2),(0),(-vacuumHeight/2)) -- 2-3
	ei_addsegment((0),(-vacuumHeight/2),(0),(vacuumHeight/2)) -- 3-0

	ei_selectsegment((0) + 0.1*sign((vacuumRadius)-(0)),(vacuumHeight/2) + 0.1*sign((vacuumHeight/2)-(vacuumHeight/2))) -- 0-1
	ei_selectsegment((vacuumRadius) + 0.1*sign((vacuumRadius)-(vacuumRadius)),(vacuumHeight/2) + 0.1*sign((-vacuumHeight/2)-(vacuumHeight/2))) -- 1-2
	ei_selectsegment((vacuumRadius) + 0.1*sign((0)-(vacuumRadius)),(-vacuumHeight/2) + 0.1*sign((-vacuumHeight/2)-(-vacuumHeight/2))) -- 2-3

	ei_setsegmentprop("",vacuumMesh,0,0,1,"zeroCond")
	ei_clearselected()

	ei_setarcsegmentprop(vacuumMesh,"",0,1,"zeroCond")
	ei_clearselected()

	ei_addblocklabel((vacuumRadius*7/8),(0))
	ei_selectlabel((vacuumRadius*7/8),(0))
	ei_setblockprop("Vacuum",0,vacuumMesh,0)
	ei_clearselected()
end

-- topElectrode
if topElectrode == 1 then
	
	drawCurve(step, 0.2, meshSizeFine, "zeroCond", electrodeSeparation+electrodeThickness*2-roundedThickness, 1, t0, tf)
	
	ei_addblocklabel((innerRadius/2),(electrodeThickness*1.5+electrodeSeparation))
	ei_selectlabel((innerRadius/2),(electrodeThickness*1.5+electrodeSeparation))
	ei_setblockprop("<No Mesh>",1,0,0)
	ei_clearselected()
	if tf ~= nil then
		ei_addnode(curveRadius+roundedWidth, (2*electrodeThickness + electrodeSeparation))
		ei_addarc(curveRadius, (2*electrodeThickness + electrodeSeparation)-roundedThickness, curveRadius+roundedWidth, (2*electrodeThickness + electrodeSeparation), (pi/2-thetaf)*180/pi, 0.2)
		ei_selectarcsegment(curveRadius+0.1, (2*electrodeThickness + electrodeSeparation)-roundedThickness+0.1)
		ei_setarcsegmentprop(meshSizeFine, "", 0, 1, "zeroCond")
	end
	
	if fullTop == 1 then
		drawCurve(step, 0.2, meshSizeFine, "zeroCond", electrodeSeparation+electrodeThickness*2+roundedThickness, -1, t0, tf)
	
		ei_addblocklabel((innerRadius/2),(electrodeThickness*2.5+electrodeSeparation))
		ei_selectlabel((innerRadius/2),(electrodeThickness*2.5+electrodeSeparation))
		ei_setblockprop("<No Mesh>",1,0,0)
		ei_clearselected()
		if tf ~= nil then
			ei_addarc(curveRadius+roundedWidth, (2*electrodeThickness + electrodeSeparation), curveRadius, (2*electrodeThickness + electrodeSeparation)+roundedThickness, (pi/2-thetaf)*180/pi, 0.2)
			ei_selectarcsegment(curveRadius+roundedWidth-0.1, (2*electrodeThickness + electrodeSeparation)+0.1)
			ei_setarcsegmentprop(meshSizeFine, "", 0, 1, "zeroCond")
		end
	end
end

-- middleElectrode
if middleElectrode == 1 then
	drawCurve(step, 0.2, meshSizeFine, "elecCond", -roundedThickness, 1, t0, tf)
	drawCurve(step, 0.2, meshSizeFine, "elecCond", roundedThickness, -1, t0, tf)
	
	if tf ~= nil then
		ei_addarc(curveRadius, -roundedThickness, curveRadius, roundedThickness, 2*(pi/2-thetaf)*180/pi, 0.2)
		ei_selectarcsegment(curveRadius+0.1, -roundedThickness+0.1)
		ei_setarcsegmentprop(meshSizeFine, "", 0, 1, "elecCond")
	end
	
	ei_addblocklabel(innerRadius/2,0)
	ei_selectlabel(innerRadius/2,0)
	ei_setblockprop("<No Mesh>",1,0,0)
	ei_clearselected()
end

-- bottomElectrode
if bottomElectrode == 1 then
	drawCurve(step, 0.2, meshSizeFine, "zeroCond", -electrodeSeparation-electrodeThickness*2+roundedThickness, -1, t0, tf)
	
	ei_addblocklabel((innerRadius/2),(-electrodeThickness*1.5-electrodeSeparation))
	ei_selectlabel((innerRadius/2),(-electrodeThickness*1.5-electrodeSeparation))
	ei_setblockprop("<No Mesh>",1,0,0)
	ei_clearselected()
	if tf ~= nil then
		ei_addnode(curveRadius+roundedWidth, -(2*electrodeThickness + electrodeSeparation))
		ei_addarc(curveRadius+roundedWidth, -(2*electrodeThickness + electrodeSeparation), curveRadius, -(2*electrodeThickness + electrodeSeparation)+roundedThickness, (pi/2-thetaf)*180/pi, 0.2)
		ei_selectarcsegment(curveRadius+roundedWidth + 0.1, -(2*electrodeThickness + electrodeSeparation) + 0.1)
		ei_setarcsegmentprop(meshSizeFine, "", 0, 1, "zeroCond")
	end
	
	if fullBottom == 1 then
		drawCurve(step, 0.2, meshSizeFine, "zeroCond", -electrodeSeparation-electrodeThickness*2-roundedThickness, 1, t0, tf)
	
		ei_addblocklabel((innerRadius/2),(-electrodeThickness*2.5-electrodeSeparation))
		ei_selectlabel((innerRadius/2),(-electrodeThickness*2.5-electrodeSeparation))
		ei_setblockprop("<No Mesh>",1,0,0)
		ei_clearselected()
		if tf ~= nil then
			ei_addarc(curveRadius, -(2*electrodeThickness + electrodeSeparation)-roundedThickness, curveRadius+roundedWidth, -(2*electrodeThickness + electrodeSeparation), (pi/2-thetaf)*180/pi, 0.2)
			ei_selectarcsegment(curveRadius+0.1, -(2*electrodeThickness + electrodeSeparation)-roundedThickness+0.1)
			ei_setarcsegmentprop(meshSizeFine, "", 0, 1, "zeroCond")
		end
	end
end

-- mainArea; used when vacuumChamber is large to increase precision in a smaller area around the electrodes
-- built here from separate python script and .txt document
-- nodes are numbered to make clear which segments connect which points
if mainArea == 1 then
	ei_addnode((0),(mainHeight/2)) -- 0
	ei_addnode((mainRadius),(mainHeight/2)) -- 1
	ei_addnode((mainRadius),(-mainHeight/2)) -- 2
	ei_addnode((0),(-mainHeight/2)) -- 3

	ei_addsegment((0),(mainHeight/2),(mainRadius),(mainHeight/2)) -- 0-1
	ei_addsegment((mainRadius),(mainHeight/2),(mainRadius),(-mainHeight/2)) -- 1-2
	ei_addsegment((mainRadius),(-mainHeight/2),(0),(-mainHeight/2)) -- 2-3
	ei_addsegment((0),(-mainHeight/2),(0),(mainHeight/2)) -- 3-0

	ei_selectsegment((0) + 0.1*sign((mainRadius)-(0)),(mainHeight/2) + 0.1*sign((mainHeight/2)-(mainHeight/2))) -- 0-1
	ei_selectsegment((mainRadius) + 0.1*sign((mainRadius)-(mainRadius)),(mainHeight/2) + 0.1*sign((-mainHeight/2)-(mainHeight/2))) -- 1-2
	ei_selectsegment((mainRadius) + 0.1*sign((0)-(mainRadius)),(-mainHeight/2) + 0.1*sign((-mainHeight/2)-(-mainHeight/2))) -- 2-3
	ei_selectsegment((0) + 0.1*sign((0)-(0)),(-mainHeight/2) + 0.1*sign((mainHeight/2)-(-mainHeight/2))) -- 3-0

	ei_setsegmentprop("",mainAreaMesh,0,0,1,"<None>")
	ei_clearselected()

	ei_setarcsegmentprop(mainAreaMesh,"",0,1,"<None>")
	ei_clearselected()

	ei_addblocklabel((mainRadius/2),(mainHeight/4))
	ei_selectlabel((mainRadius/2),(mainHeight/4))
	ei_setblockprop("Vacuum",0,mainAreaMesh,0)
	ei_clearselected()
end

-- upperBlock; creates a rectangle between the upper and middle electrodes for block integral calculations
if upperBlock == 1 then
	-- because the segment is drawn from the surface y at r = 0 not at r = blockRadius, the segment is shifted a couple mm from blockRadius
	-- it is still vertical because the electrodes are mirror images of each other
	-- the problem could be fixed by drawing a segment from the top of the EDM cell to the bottom, but this requires adding additional block labels to the right of the segment as it passes through each electrode
	ei_addsegment(blockRadius, electrodeThickness, blockRadius, electrodeThickness+electrodeSeparation)
	
	ei_addblocklabel((blockRadius/2),(mainHeight/4))
	ei_selectlabel((blockRadius/2),(mainHeight/4))
	ei_setblockprop("Vacuum",0,mainAreaMesh,0)
	ei_clearselected()
end

-- lowerBlock; creates a rectangle between the lower and middle electrodes for block integral calculations
if lowerBlock == 1 then
	-- because the segment is drawn from the surface y at r = 0 not at r = blockRadius, the segment is shifted a couple mm from blockRadius
	-- it is still vertical because the electrodes are mirror images of each other
	-- the problem could be fixed by drawing a segment from the top of the EDM cell to the bottom, but this requires adding additional block labels to the right of the segment as it passes through each electrode
	ei_addsegment(blockRadius, -electrodeThickness, blockRadius, -electrodeThickness-electrodeSeparation)
	
	ei_addblocklabel((blockRadius/2),(-mainHeight/4))
	ei_selectlabel((blockRadius/2),(-mainHeight/4))
	ei_setblockprop("Vacuum",0,mainAreaMesh,0)
	ei_clearselected()
end

--###################################################################################
--############################# End of Model Build Info #############################
--###################################################################################

--############################ Postprocessing Functions #############################

-- records voltage and electric field information over area
function gather(xmin, xmax, xstep, ymin, ymax, ystep)
	map = {}

	for i = 1, (xmax-xmin)/xstep, 1 do
		map[i] = {}
		x = xmin + xstep*(i-1)
		
		for j = 1, (ymax-ymin)/ystep do
			y = ymin + ystep*(j-1)
			
			V, Dx, Dy, Ex, Ey, ex, ey, nrg = eo_getpointvalues(x, y)
			V, Ex, Ey = V or 0, Ex or 0, Ey or 0 -- 0 if no value exists
			
			map[i][j] = {x, y, V, Ex, Ey}
		end
	end
	return map
end

-- prepares a text file with information from gather() in a specific format used by PENTRACK simulations
-- files work with PENTRACK, but few tests have been done using these files
function pentrack(map, handle)
	write(handle, getn(map), " 1 ", getn(map[1]), " 2\n1 X*100 [METRE]\n2 Z*100 [METRE]\n3 RBX*1E4 [TESLA]\n4 RBZ*1E4 [TESLA]\n5 RV [JOULE]\n 0\n")
	for i = 1, getn(map) do
		for j = 1, getn(map[i]) do
			write(handle, map[i][j][1]/1000,"\t",map[i][j][2]/1000,"\t0\t0\t",sqrt(map[i][j][4]^2+map[i][j][5]^2),"\n")
		end
	end
end

-- gathers field information along a horizontal contour
function horizontalContour(xmax, xstep, y, handle)
	write(handle, "k,vH,vR,x,y,E,Emin,Emax,(Emax-Emin)/Emax\n")
	Emin = 1e9
	Emax = -1e9
	for x = 0, xmax, xstep do
		V, Dx, Dy, Ex, Ey, ex, ey, nrg = eo_getpointvalues(x, y)
		Ex = Ex or 0
		Ey = Ey or 0
		E = sqrt(Ex^2+Ey^2)
		Emin = min(Emin, E)
		Emax = max(Emax, E)
		write(handle, k0, ",", x, ",", y, ",", E, ",", Emin, ",", Emax, ",", (Emax-Emin)/Emax,"\n")
	end
end

-- for a range in y, finds the x values at which (Emax-Emin)/Emax exceeds a certain maximum value
function deviationContour(maxDeviation, xstep, ymin, ymax, ystep, handle)
	write(handle, "k0,x,y,dev\n")
	for y = ymin, ymax, ystep do
	
		dev = 0
		x = 0
		
		Emin = 1e9
		Emax = -1e9
		while dev < maxDeviation do
			V, Dx, Dy, Ex, Ey, ex, ey, nrg = eo_getpointvalues(x, y)
			Ex = Ex or 0
			Ey = Ey or 0
			E = sqrt(Ex^2+Ey^2)
			Emin = min(Emin, E)
			Emax = max(Emax, E)
			dev = (Emax-Emin)/Emax
			--messagebox(Ex..", "..Ey..", "..E..", "..Emin..", "..Emax..", "..dev)
			x = x + xstep
		end
		write(handle, k0..","..x..","..y..","..dev.."\n")
	end
end

-- gathers er and ez information in main blocks
function post(handle)
	eo_selectblock(blockRadius/2, mainHeight/4)
	er, ez = eo_blockintegral(4)
	write(handle, k0..","..electrodeSeparation..","..vacuumHeight..","..vacuumRadius..","..er..","..ez..","..er/ez.."\n")
end

--################################# Postprocessing ##################################
-- Postprocessing Variables
docName = electrode..".fee"
pentrackFile = "pentrack"..electrode..".txt"
fullMapFile = "fullMap"..electrode..".csv"
hContourFile1 = "hContourFile"..electrode.."High.csv"
hContourFile2 = "hContourFile"..electrode.."Med.csv"
hContourFile3 = "hContourFile"..electrode.."Low.csv"
postFile = "post"..electrode..".csv"
deviationFile = "deviation"..electrode..".csv"

densityMax = 3e6 -- for electric field gradient
xmin, xmax, xstep, ymin, ymax, ystep = 0, 180, 3, -210, 210, 3 -- for PENTRACK gathering

-- Inspect Model
ei_zoom(0, -mainHeight/2, mainRadius, mainHeight/2)
ei_saveas(docName)
ei_analyze()
ei_loadsolution()

-- FEMM Postprocessing
eo_showdensityplot(1,0,2,densityMax,0)
eo_zoom(0, -mainHeight/2, mainRadius, mainHeight/2)

map = gather(xmin, xmax, xstep, ymin, ymax, ystep)

-- Pentrack
handle=openfile(pentrackFile,'w')
pentrack(map, handle)
closefile(handle)

-- Deviation Contour
handle = openfile(deviationFile, 'w')
deviationContour(0.05, 0.05, electrodeThickness, electrodeThickness+electrodeSeparation, 0.05, handle)
closefile(handle)

-- Horizontal Contours
-- Top Contour
handle = openfile(hContourFile1, 'w')
horizontalContour(180, 1, 3*mainHeight/8, handle)
closefile(handle)

-- Middle Contour
handle = openfile(hContourFile2, 'w')
horizontalContour(180, 1, mainHeight/4, handle)
closefile(handle)

-- Bottom Contour
handle = openfile(hContourFile3, 'w')
horizontalContour(180, 1, mainHeight/8, handle)
closefile(handle)

-- Post Info
handle = openfile(postFile, 'a')
post(handle)
closefile(handle)
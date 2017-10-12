--Electrode FEA v2

--New features:	1. Added study mode option, one can choose to do the o-ring trench study or just a flat electrode.
--				2. Added gas filling port on the top electrode. One can choose to include this feature or not.

--Updates: 1. fixed some geometry errors, now the code can correctly generate customized parameters combination
--         2. fixed the boundary property errors.
--		   3. changed the material of insulator from porcelain to glass
--         4. optimized the some of the variables and formular being used
--         5. added more comments

--v1.0 was coded by Rex Chen


-------------------Programming Structure--------------------
-- 0. User Input
-- 1. Specify the value of each parameter including customizable parameters and fixed parameters
-- 2. Define the problem and setup the environment
-- 3. Build geometry
--		3a. Build the vacuum chamber
--		3b. Build the lower electrode
--		3c. Build the insulator
--		3d. BUild the upper electrode
-- 4. Add block labels and assign materials
-- 5. Specify boundary conditions and mesh size
-- 6. Start meshing and solving
-- 7. Postprocessor (add contour, save plot...)
---------------------------------------------------------------

clearconsole()
newdocument(1)

--###################################################################################
--################################ User Input Area ##################################
--###################################################################################

--################# Study Mode ######################################################
boolOringTrench = 1 -- 1: O-ring trench optimization study. 
					-- 0: NO o-ring trench, square oring
					
boolGasFilling = 1  -- 1: include the gas filling port on top electrode. 
					-- 0: NO gas filling port
--###################################################################################

--############### Input Variable ####################################################
-- Refer to Martin's thesis for the meaning of each parameter
-- s should be non-zero
r = 5 --electrode radius around insulator groove
s = 1 --distance between side of trench and insulator
d = 2 --height of electrode wall around insulator, NOT INCLUDING radius(r)
a = 0 --electrode wall to insulator angle
--###################################################################################

--################ Problem Setting ##################################################
highVoltage = 100000 -- high voltage V
lowVoltage = 0 -- low voltage V
electrodeSeparation = 65 --distance between upper and lower electrode, default: 65
insulatorInnerRadius = 50.8 --insulator inner radius. Default: 50.8
squareOringHeight = 10 --only used when in "NO o-ring trench" study mode
--###################################################################################

--################### FEA Setting ###################################################
meshSizeCoarse = 10
meshSizeMedium = 1
meshSizeFine = 0.1
meshSizeUltraFine = 0.01
contourPtNum = 1000 --specify the number of points taken on the contour line in postprocessor
--###################################################################################

--###################################################################################
--############################# End of User Input ###################################
--###################################################################################


-- 1. Specify the value of each non customizable fixed parameters
--NON-CUSTOMIZABLE parameters
--these parameters are either have standard values or related to the geometry of the HV device setup 
vacuumWidth = 374.65 --vacuum tank width
vacuumHeight = 838.2 --vacuum tank height
oRingGrooveWidth = 3.8100 --o-ring groove s
oRingGrooveHeight = 2.8956 --o-ring groove height
aRad = a * 0.0174532925 --convert electrode wall to insulator angle form degree to radius
wallShift = d*tan(aRad) --distance bottom of electrode wall is shifted by
insulatorThickness = 6.35 --insulator thickness (outer radius - inner radius)
oRingCenter =(insulatorInnerRadius+insulatorThickness+insulatorInnerRadius)/2 --center of o-ring groove
bottom = 294.95 -- height from the bottom of vacuum chamber to bottom of lower electrode
electrodeThickness = 20 --thicknes of each electrode (discounting portruding curves and profiles)

--2. Define the problem and setup the environment
ei_probdef("millimeters","axi",1e-8,1,10)
ei_addconductorprop("elecCond",highVoltage,0,1)
ei_addconductorprop("zeroCond",lowVoltage,0,1)

ei_addmaterial("Glass",4.7,4.7,0) --from WikiPedia
ei_addmaterial("Vacuum",1,1,0)

-- 3. Build geometry
-- 3a. Build the vacuum chamber
ei_addnode(0,0)
ei_addnode(0,vacuumHeight)
ei_addnode(vacuumWidth,0)
ei_addnode(vacuumWidth,vacuumHeight)
ei_clearselected()

-- lines
ei_addsegment(0,0,0,vacuumHeight)
ei_addsegment(0,vacuumHeight,vacuumWidth,vacuumHeight)
ei_addsegment(0,0,vacuumWidth,0)
ei_addsegment(vacuumWidth,0,vacuumWidth,vacuumHeight)
ei_clearselected()

-- 3b. build lower electrode
--fixed nodes
ei_addnode(0,bottom + electrodeThickness)
ei_addnode(0,bottom +13 + 2.7)
ei_addnode(5,bottom +13)
ei_addnode(5,bottom)
ei_addnode(32,bottom)
ei_addnode(37,bottom - 5)
ei_addnode(37,bottom - 5 - 1)
ei_addnode(42,bottom - 11)
ei_addnode(80.7906 + (insulatorInnerRadius - 50.8),bottom - 11)
ei_addnode(84.6325 + (insulatorInnerRadius - 50.8),bottom - 11 - 1.8)
ei_addnode(100 + (insulatorInnerRadius - 50.8),bottom - 20)
ei_addnode(120 + (insulatorInnerRadius - 50.8),bottom)
ei_addnode(100 + (insulatorInnerRadius - 50.8),bottom + electrodeThickness)
ei_clearselected()

--fixed lines
ei_addsegment(0,bottom + electrodeThickness,0,bottom +13 + 2.7)
ei_addsegment(0,bottom +13 + 2.7,5,bottom +13)
ei_addsegment(5,bottom +13,5,bottom)
ei_addsegment(5,bottom,32,bottom)
ei_addarc(37,bottom - 5,32,bottom,90,0.2)
ei_addsegment(37,bottom - 5,37,bottom - 5 - 1)
ei_addarc(37,bottom - 5 - 1,42,bottom - 11,90,0.2)
ei_addsegment(42,bottom - 11,80.7906+(insulatorInnerRadius-50.8),bottom - 11)
ei_addarc(84.6325+(insulatorInnerRadius-50.8),bottom - 11 - 1.8,80.7906+(insulatorInnerRadius-50.8),bottom - 11,50.2082,0.2)
ei_addarc(84.6325+(insulatorInnerRadius-50.8),bottom - 11 - 1.8,100+(insulatorInnerRadius-50.8),bottom - 20,50.2082,0.2)
ei_addarc(100+(insulatorInnerRadius-50.8),bottom - 20,120+(insulatorInnerRadius-50.8),bottom,90,0.2)
ei_addarc(120+(insulatorInnerRadius-50.8),bottom,100+(insulatorInnerRadius-50.8),bottom + electrodeThickness,90,0.2)
ei_clearselected()

--o-ring trench area
if boolOringTrench == 1 then
	ei_addnode(insulatorInnerRadius+insulatorThickness+s+wallShift+r*cos(aRad), bottom + electrodeThickness)
	ei_addnode(insulatorInnerRadius+insulatorThickness+s+wallShift,bottom + electrodeThickness -(r - r*sin(aRad)))
	ei_addnode(insulatorInnerRadius+insulatorThickness+s,bottom + electrodeThickness -(r - r*sin(aRad))-d)

	ei_addnode(insulatorInnerRadius-s-wallShift-r*cos(aRad), bottom + electrodeThickness)
	ei_addnode(insulatorInnerRadius-s-wallShift,bottom + electrodeThickness -(r - r*sin(aRad)))
	ei_addnode(insulatorInnerRadius-s,bottom + electrodeThickness -(r - r*sin(aRad))-d)

	ei_addnode(oRingCenter+(oRingGrooveWidth/2),bottom + electrodeThickness -(r - r*sin(aRad))-d)
	ei_addnode(oRingCenter-(oRingGrooveWidth/2),bottom + electrodeThickness -(r - r*sin(aRad))-d)
	ei_addnode(oRingCenter+(oRingGrooveWidth/2),bottom + electrodeThickness -(r - r*sin(aRad))-d-oRingGrooveHeight)
	ei_addnode(oRingCenter-(oRingGrooveWidth/2),bottom + electrodeThickness -(r - r*sin(aRad))-d-oRingGrooveHeight)
	ei_clearselected()

	ei_addsegment(100+(insulatorInnerRadius-50.8),bottom + electrodeThickness,insulatorInnerRadius+insulatorThickness+s+wallShift+r*cos(aRad), bottom + electrodeThickness)
	ei_addsegment(0,bottom + electrodeThickness, insulatorInnerRadius-s-wallShift-r*cos(aRad), bottom + electrodeThickness)

	ei_addarc(insulatorInnerRadius-s-wallShift,bottom + electrodeThickness -(r - r*sin(aRad)),insulatorInnerRadius-s-wallShift-r*cos(aRad), bottom + electrodeThickness, 90-a,0.2)
	ei_addarc(insulatorInnerRadius+insulatorThickness+s+wallShift+r*cos(aRad), bottom + electrodeThickness,insulatorInnerRadius+insulatorThickness+s+wallShift,bottom + electrodeThickness -(r - r*sin(aRad)), 90-a,0.2)

	ei_addsegment(insulatorInnerRadius+insulatorThickness+s+wallShift,bottom + electrodeThickness -(r - r*sin(aRad)), insulatorInnerRadius+insulatorThickness+s,bottom + electrodeThickness -(r - r*sin(aRad))-d)
	ei_addsegment(insulatorInnerRadius-s-wallShift,bottom + electrodeThickness -(r - r*sin(aRad)),insulatorInnerRadius-s,bottom + electrodeThickness -(r - r*sin(aRad))-d)

	ei_addsegment(insulatorInnerRadius+insulatorThickness+s,bottom + electrodeThickness -(r - r*sin(aRad))-d,oRingCenter+(oRingGrooveWidth/2),bottom + electrodeThickness -(r - r*sin(aRad))-d)
	ei_addsegment(insulatorInnerRadius-s,bottom + electrodeThickness -(r - r*sin(aRad))-d,oRingCenter-(oRingGrooveWidth/2),bottom + electrodeThickness -(r - r*sin(aRad))-d)

	ei_addsegment(oRingCenter+(oRingGrooveWidth/2),bottom + electrodeThickness -(r - r*sin(aRad))-d,oRingCenter+(oRingGrooveWidth/2),bottom + electrodeThickness -(r - r*sin(aRad))-d-oRingGrooveHeight)
	ei_addsegment(oRingCenter-(oRingGrooveWidth/2),bottom + electrodeThickness -(r - r*sin(aRad))-d,oRingCenter-(oRingGrooveWidth/2),bottom + electrodeThickness -(r - r*sin(aRad))-d-oRingGrooveHeight)
	ei_addsegment(oRingCenter-(oRingGrooveWidth/2),bottom + electrodeThickness -(r - r*sin(aRad))-d-oRingGrooveHeight,oRingCenter+(oRingGrooveWidth/2),bottom + electrodeThickness -(r - r*sin(aRad))-d-oRingGrooveHeight)

	ei_clearselected()	
	
else
	--no o-ring trench
	ei_addsegment(0,bottom + electrodeThickness,100 + (insulatorInnerRadius - 50.8),bottom + electrodeThickness)
	ei_clearselected()
end

-- -- 3c. build the  insulator
if boolOringTrench == 1 then
--nodes
	ei_addnode(insulatorInnerRadius+insulatorThickness,bottom + electrodeThickness + electrodeSeparation+(r - r*sin(aRad))+d)
	ei_addnode(insulatorInnerRadius,bottom + electrodeThickness + electrodeSeparation+(r - r*sin(aRad))+d)
	ei_addnode(insulatorInnerRadius+insulatorThickness,bottom + electrodeThickness - (r - r*sin(aRad)) - d)
	ei_addnode(insulatorInnerRadius,bottom + electrodeThickness -(r - r*sin(aRad)) - d)
	ei_clearselected()
--lines
	ei_addsegment(insulatorInnerRadius+insulatorThickness,bottom + electrodeThickness + electrodeSeparation+(r - r*sin(aRad))+d,insulatorInnerRadius,bottom + electrodeThickness + electrodeSeparation+(r - r*sin(aRad))+d)
	ei_addsegment(insulatorInnerRadius,bottom + electrodeThickness + electrodeSeparation+(r - r*sin(aRad))+d,insulatorInnerRadius,bottom + electrodeThickness -(r - r*sin(aRad))-d)
	ei_addsegment(insulatorInnerRadius+insulatorThickness,bottom + electrodeThickness -(r - r*sin(aRad))-d,insulatorInnerRadius,bottom + electrodeThickness -(r - r*sin(aRad))-d)
	ei_addsegment(insulatorInnerRadius+insulatorThickness,bottom + electrodeThickness +electrodeSeparation+(r - r*sin(aRad))+d,insulatorInnerRadius+insulatorThickness,bottom + electrodeThickness -(r - r*sin(aRad))-d)
	ei_clearselected()
else
--insulator contact points on electrode surface
	ei_addnode(insulatorInnerRadius+insulatorThickness,bottom + electrodeThickness)
	ei_addnode(insulatorInnerRadius,bottom + electrodeThickness)
	ei_addnode(insulatorInnerRadius+insulatorThickness,bottom + electrodeThickness + electrodeSeparation)
	ei_addnode(insulatorInnerRadius,bottom + electrodeThickness + electrodeSeparation)
--oring
	ei_addnode(insulatorInnerRadius+insulatorThickness,bottom + electrodeThickness + squareOringHeight)	
	ei_addnode(insulatorInnerRadius,bottom + electrodeThickness + squareOringHeight)
	ei_addnode(insulatorInnerRadius+insulatorThickness,bottom + electrodeThickness + electrodeSeparation - squareOringHeight)
	ei_addnode(insulatorInnerRadius,bottom + electrodeThickness + electrodeSeparation - squareOringHeight)
	ei_clearselected()

	--lines
	ei_addsegment(insulatorInnerRadius+insulatorThickness,bottom + electrodeThickness,insulatorInnerRadius+insulatorThickness,bottom + electrodeThickness + electrodeSeparation)
	ei_addsegment(insulatorInnerRadius,bottom + electrodeThickness,insulatorInnerRadius,bottom + electrodeThickness + electrodeSeparation)
	ei_addsegment(insulatorInnerRadius+insulatorThickness,bottom + electrodeThickness + squareOringHeight,insulatorInnerRadius,bottom + electrodeThickness + squareOringHeight)
	ei_addsegment(insulatorInnerRadius+insulatorThickness,bottom + electrodeThickness + electrodeSeparation - squareOringHeight,insulatorInnerRadius,bottom + electrodeThickness + electrodeSeparation - squareOringHeight)
end

-- 3d. build the upper electrode
-- fixed nodes
ei_addnode(0,bottom + 2*electrodeThickness + electrodeSeparation)
ei_addnode(0,bottom + electrodeThickness + electrodeSeparation)
ei_addnode(78.0911+(insulatorInnerRadius-50.8),bottom + 2*electrodeThickness + electrodeSeparation)
ei_addnode(80.0828+(insulatorInnerRadius-50.8),bottom + electrodeThickness + electrodeSeparation + 21.818)
ei_addnode(100+(insulatorInnerRadius-50.8),bottom + electrodeThickness + electrodeSeparation + 40)
ei_addnode(120+(insulatorInnerRadius-50.8),bottom + electrodeThickness + electrodeSeparation + 20)
ei_addnode(100+(insulatorInnerRadius-50.8),bottom + electrodeThickness + electrodeSeparation)
ei_clearselected()

--fixed lines
ei_addsegment(0,bottom + 2*electrodeThickness + electrodeSeparation,78.0911+(insulatorInnerRadius-50.8),bottom + 2* electrodeThickness + electrodeSeparation)
ei_addarc(78.0911+(insulatorInnerRadius-50.8),bottom + 2* electrodeThickness + electrodeSeparation,80.0828+(insulatorInnerRadius-50.8),bottom + electrodeThickness + electrodeSeparation + 21.818,84.7841,0.2)
ei_addarc(100+(insulatorInnerRadius-50.8),bottom + electrodeThickness + electrodeSeparation + 40,80.0828+(insulatorInnerRadius-50.8),bottom + electrodeThickness + electrodeSeparation + 21.818,84.7841,0.2)
ei_addarc(120+(insulatorInnerRadius-50.8),bottom + 2*electrodeThickness + electrodeSeparation,100+(insulatorInnerRadius-50.8),bottom + electrodeThickness + electrodeSeparation + 40,90,0.2)
ei_addarc(100+(insulatorInnerRadius-50.8),bottom + electrodeThickness + electrodeSeparation,120+(insulatorInnerRadius-50.8),bottom + electrodeThickness + electrodeSeparation + 20,90,0.2)
ei_clearselected()

-- o-ring trench area
if boolOringTrench == 1 then
	ei_addnode(insulatorInnerRadius+insulatorThickness+s+wallShift+r*cos(aRad), bottom + electrodeThickness + electrodeSeparation)
	ei_addnode(insulatorInnerRadius+insulatorThickness+s+wallShift,bottom + electrodeThickness + electrodeSeparation+(r - r*sin(aRad)))
	ei_addnode(insulatorInnerRadius+insulatorThickness+s,bottom + electrodeThickness + electrodeSeparation+(r - r*sin(aRad))+d)

	ei_addnode(insulatorInnerRadius-s-wallShift-r*cos(aRad), bottom + electrodeThickness + electrodeSeparation)
	ei_addnode(insulatorInnerRadius-s-wallShift,bottom + electrodeThickness + electrodeSeparation+(r - r*sin(aRad)))
	ei_addnode(insulatorInnerRadius-s,bottom + electrodeThickness + electrodeSeparation+(r - r*sin(aRad))+d)

	ei_addnode(oRingCenter+(oRingGrooveWidth/2),bottom + electrodeThickness + electrodeSeparation+(r - r*sin(aRad))+d)
	ei_addnode(oRingCenter-(oRingGrooveWidth/2),bottom + electrodeThickness + electrodeSeparation+(r - r*sin(aRad))+d)
	ei_addnode(oRingCenter+(oRingGrooveWidth/2),bottom + electrodeThickness + electrodeSeparation+(r - r*sin(aRad))+d+oRingGrooveHeight)
	ei_addnode(oRingCenter-(oRingGrooveWidth/2),bottom + electrodeThickness + electrodeSeparation+(r - r*sin(aRad))+d+oRingGrooveHeight)

	ei_clearselected()

	ei_addsegment(100+(insulatorInnerRadius-50.8),bottom + electrodeThickness + electrodeSeparation,insulatorInnerRadius+insulatorThickness+s+wallShift+r*cos(aRad), bottom + electrodeThickness + electrodeSeparation)
	ei_addsegment(0,bottom + electrodeThickness + electrodeSeparation, insulatorInnerRadius-s-wallShift-r*cos(aRad), bottom + electrodeThickness + electrodeSeparation)

	ei_addarc(insulatorInnerRadius-s-wallShift-r*cos(aRad), bottom + electrodeThickness + electrodeSeparation, insulatorInnerRadius-s-wallShift,bottom + electrodeThickness + electrodeSeparation+(r - r*sin(aRad)),90-a,0.2)
	ei_addarc(insulatorInnerRadius+insulatorThickness+s+wallShift,bottom + electrodeThickness + electrodeSeparation+(r - r*sin(aRad)),insulatorInnerRadius+insulatorThickness+s+wallShift+r*cos(aRad), bottom + electrodeThickness + electrodeSeparation, 90-a,0.2)

	ei_addsegment(insulatorInnerRadius+insulatorThickness+s+wallShift,bottom + electrodeThickness + electrodeSeparation+(r - r*sin(aRad)), insulatorInnerRadius+insulatorThickness+s,bottom + electrodeThickness + electrodeSeparation+(r - r*sin(aRad))+d)
	ei_addsegment(insulatorInnerRadius-s-wallShift,bottom + electrodeThickness + electrodeSeparation+(r - r*sin(aRad)),insulatorInnerRadius-s,bottom + electrodeThickness + electrodeSeparation+(r - r*sin(aRad))+d)

	ei_addsegment(insulatorInnerRadius+insulatorThickness+s,bottom + electrodeThickness + electrodeSeparation+(r - r*sin(aRad))+d,oRingCenter+(oRingGrooveWidth/2),bottom + electrodeThickness + electrodeSeparation+(r - r*sin(aRad))+d)
	ei_addsegment(insulatorInnerRadius-s,bottom + electrodeThickness + electrodeSeparation+(r - r*sin(aRad))+d,oRingCenter-(oRingGrooveWidth/2),bottom + electrodeThickness + electrodeSeparation+(r - r*sin(aRad))+d)

	ei_addsegment(oRingCenter+(oRingGrooveWidth/2),bottom + electrodeThickness + electrodeSeparation+(r - r*sin(aRad))+d,oRingCenter+(oRingGrooveWidth/2),bottom + electrodeThickness + electrodeSeparation+(r - r*sin(aRad))+d+oRingGrooveHeight)
	ei_addsegment(oRingCenter-(oRingGrooveWidth/2),bottom + electrodeThickness + electrodeSeparation+(r - r*sin(aRad))+d,oRingCenter-(oRingGrooveWidth/2),bottom + electrodeThickness + electrodeSeparation+(r - r*sin(aRad))+d+oRingGrooveHeight)
	ei_addsegment(oRingCenter-(oRingGrooveWidth/2),bottom + electrodeThickness + electrodeSeparation+(r - r*sin(aRad))+d+oRingGrooveHeight,oRingCenter+(oRingGrooveWidth/2),bottom + electrodeThickness + electrodeSeparation+(r - r*sin(aRad))+d+oRingGrooveHeight)

	ei_clearselected()	
else
	ei_addsegment(0,bottom + electrodeThickness + electrodeSeparation,100+(insulatorInnerRadius-50.8),bottom + electrodeThickness + electrodeSeparation)
	ei_clearselected()	
end

-- gas filling port
if boolGasFilling == 1 then
	ei_addnode(4.5,bottom + electrodeThickness + electrodeSeparation)
	ei_addnode(1.5,bottom + electrodeThickness + electrodeSeparation + 3)
	ei_addnode(1.5,bottom + electrodeThickness + electrodeSeparation + 3 + 0.5367)
	ei_addnode(8.9281,bottom + electrodeThickness + electrodeSeparation + 3 + 0.5367 + 4.4633)
	ei_addnode(8.9281,bottom + electrodeThickness + electrodeSeparation + 3 + 0.5367 + 4.4633 + 1)
	ei_addnode(9.4797,bottom + electrodeThickness + electrodeSeparation + 3 + 0.5367 + 4.4633 + 1)
	ei_addnode(9.8091,bottom + electrodeThickness + electrodeSeparation + 3 + 0.5367 + 4.4633 + 1 + 11.0049)
	ei_clearselected()
	
	ei_addarc(1.5,bottom + electrodeThickness + electrodeSeparation + 3,4.5,bottom + electrodeThickness + electrodeSeparation, 90,0.2)
	ei_addsegment(1.5,bottom + electrodeThickness + electrodeSeparation + 3,1.5,bottom + electrodeThickness + electrodeSeparation + 3 + 0.5367)
	ei_addsegment(1.5,bottom + electrodeThickness + electrodeSeparation + 3 + 0.5367, 8.9281,bottom + electrodeThickness + electrodeSeparation + 3 + 0.5367 + 4.4633)
	ei_addsegment(8.9281,bottom + electrodeThickness + electrodeSeparation + 3 + 0.5367 + 4.4633,8.9281,bottom + electrodeThickness + electrodeSeparation + 3 + 0.5367 + 4.4633 + 1)
	ei_addsegment(8.9281,bottom + electrodeThickness + electrodeSeparation + 3 + 0.5367 + 4.4633 + 1,9.4797,bottom + electrodeThickness + electrodeSeparation + 3 + 0.5367 + 4.4633 + 1)
	ei_addsegment(9.4797,bottom + electrodeThickness + electrodeSeparation + 3 + 0.5367 + 4.4633 + 1,9.8091,bottom + electrodeThickness + electrodeSeparation + 3 + 0.5367 + 4.4633 + 1 + 11.0049)	
	ei_clearselected()
	
	ei_selectsegment(0.1, bottom + electrodeThickness + electrodeSeparation)
	ei_deleteselectedsegments()
	ei_clearselected()	
end
	
-- 4. Add block labels and assign materials
-- block labels
ei_addblocklabel(100+(insulatorInnerRadius-50.8),bottom + electrodeThickness + electrodeSeparation + 20) --upper electrode
ei_addblocklabel(13.1,bottom+electrodeThickness+electrodeSeparation/2) --gas cell
ei_addblocklabel(140.5+(insulatorInnerRadius-50.8),bottom+electrodeThickness+electrodeSeparation/2) --safety box
ei_addblocklabel(100+(insulatorInnerRadius-50.8),bottom) --lower electrode
ei_addblocklabel((insulatorInnerRadius+insulatorInnerRadius+insulatorThickness)/2,bottom+electrodeThickness+electrodeSeparation/2) --insulator
if boolOringTrench == 1 then
	-- oring groove
	ei_addblocklabel((insulatorInnerRadius+insulatorInnerRadius+insulatorThickness)/2,bottom + electrodeThickness + electrodeSeparation+(r - r*sin(aRad))+d+(oRingGrooveHeight/2)) --upper o-ring
	ei_addblocklabel((insulatorInnerRadius+insulatorInnerRadius+insulatorThickness)/2,bottom + electrodeThickness -(r - r*sin(aRad))-d-(oRingGrooveHeight/2)) --lower o-ring
else
	-- square o-ring
	ei_addblocklabel(insulatorInnerRadius + 0.5*insulatorThickness,bottom + electrodeThickness + 0.5*squareOringHeight)
	ei_addblocklabel(insulatorInnerRadius + 0.5*insulatorThickness,bottom + electrodeThickness + electrodeSeparation - 0.5*squareOringHeight)
end
ei_clearselected()

-- assign materials
ei_selectlabel(100+(insulatorInnerRadius-50.8),bottom + electrodeThickness + electrodeSeparation + 20) --upper electrode
ei_selectlabel(100+(insulatorInnerRadius-50.8),bottom) --lower electrode
ei_setblockprop("<No Mesh>",1,0,0)
ei_clearselected()

ei_selectlabel(13.1,bottom+electrodeThickness+electrodeSeparation/2) --gas cell
ei_setblockprop("Vacuum",0,meshSizeMedium,0)
ei_clearselected()	

ei_selectlabel(140.5+(insulatorInnerRadius-50.8),bottom+electrodeThickness+electrodeSeparation/2) --safety box
ei_setblockprop("Vacuum",0,meshSizeCoarse,0)
ei_clearselected()

ei_selectlabel((insulatorInnerRadius+insulatorInnerRadius+insulatorThickness)/2,bottom+electrodeThickness+electrodeSeparation/2) --insulator
ei_setblockprop("Glass",0,meshSizeMedium,0)
ei_clearselected()

if boolOringTrench == 1 then
	-- oring groove
	ei_selectlabel((insulatorInnerRadius+insulatorInnerRadius+insulatorThickness)/2,bottom + electrodeThickness + electrodeSeparation+(r - r*sin(aRad))+d+(oRingGrooveHeight/2)) --upper o-ring
	ei_selectlabel((insulatorInnerRadius+insulatorInnerRadius+insulatorThickness)/2,bottom + electrodeThickness -(r - r*sin(aRad))-d-(oRingGrooveHeight/2)) --lower o-ring
	ei_setblockprop("Vacuum",0,meshSizeFine,0)
else
	-- square o-ring
	ei_getmaterial("Rubber")
	ei_selectlabel(insulatorInnerRadius + 0.5*insulatorThickness,bottom + electrodeThickness + 0.5*squareOringHeight)
	ei_selectlabel(insulatorInnerRadius + 0.5*insulatorThickness,bottom + electrodeThickness + electrodeSeparation - 0.5*squareOringHeight)
	ei_setblockprop("Rubber",0,meshSizeFine,0)
end
ei_clearselected()

-- 5. Specify boundary conditions
-- To specify boundary condition, one must firstly select all the line/arc
-- segments the body contained. To select the proper line/arc segment, 
-- choose a node and shift by a tiny number (like 0.01) to either x or y direction

-- vacuum chamber
ei_selectsegment(0+0.1,0)
ei_selectsegment(0+0.1,vacuumHeight)
ei_selectsegment(vacuumWidth,0+0.1)
ei_setsegmentprop("",meshSizeCoarse,0,0,1,"zeroCond")
ei_clearselected()

--lower electrode
--fixed lines and curves
ei_selectsegment(0+0.01,bottom +13 + 2.7-0.01)
ei_selectsegment(5,bottom +13-0.1)
ei_selectsegment(5+0.01,bottom)
ei_selectsegment(37,bottom - 5-0.1)
ei_selectsegment(42+0.01,bottom - 11)
ei_setsegmentprop("",meshSizeFine,0,0,1,"elecCond")
ei_clearselected()

ei_selectarcsegment(32,bottom)
ei_selectarcsegment(37,bottom - 11)
ei_selectarcsegment(80 - 0.1+(insulatorInnerRadius-50.8),bottom - 11)
ei_selectarcsegment(100 - 0.1+(insulatorInnerRadius-50.8),bottom - 20)
ei_selectarcsegment(100 + 0.1+(insulatorInnerRadius-50.8),bottom - 20)
ei_selectarcsegment(100 + 0.1+(insulatorInnerRadius-50.8),bottom + electrodeThickness)
ei_setarcsegmentprop(meshSizeFine,"",0,1,"elecCond")
ei_clearselected()

-- O-ring trench area
if boolOringTrench == 1 then
	-- two flat surfaces on the electrode
	ei_selectsegment(insulatorInnerRadius+insulatorThickness+s+wallShift+r*cos(aRad) + 0.01, bottom + electrodeThickness)	
	ei_selectsegment(insulatorInnerRadius-s-wallShift-r*cos(aRad)-0.01, bottom + electrodeThickness)
	
	-- O-ring groove
	ei_selectsegment(oRingCenter+(oRingGrooveWidth/2),bottom + electrodeThickness -(r - r*sin(aRad))-d-0.01)
	ei_selectsegment(oRingCenter-(oRingGrooveWidth/2),bottom + electrodeThickness -(r - r*sin(aRad))-d-0.01)
	ei_selectsegment(oRingCenter+(oRingGrooveWidth/2)-0.01,bottom + electrodeThickness -(r - r*sin(aRad))-d-oRingGrooveHeight)
	ei_setsegmentprop("",meshSizeFine,0,0,1,"elecCond")
	ei_clearselected()
	
	-- line of contact of insulator and the trench
	ei_selectsegment(insulatorInnerRadius+insulatorThickness - 0.01,bottom + electrodeThickness -(r - r*sin(aRad))-d)
	ei_selectsegment(insulatorInnerRadius + 0.01,bottom + electrodeThickness -(r - r*sin(aRad)) - d)
	ei_setsegmentprop("",meshSizeUltraFine,0,0,1,"elecCond")
	ei_clearselected()

	-- trench wall, d
	if d ~= 0 then
		ei_selectsegment(insulatorInnerRadius+insulatorThickness+s+wallShift+0.01,bottom + electrodeThickness -(r - r*sin(aRad))-0.01)
		ei_selectsegment(insulatorInnerRadius-s-wallShift-0.01,bottom + electrodeThickness - (r - r*sin(aRad)) -0.01)	
		ei_setsegmentprop("",meshSizeUltraFine,0,0,1,"elecCond")
		ei_clearselected()	
	end
	
	-- seperation, s
	if s ~= 0 then
		ei_selectsegment(insulatorInnerRadius+insulatorThickness+0.01,bottom + electrodeThickness -(r - r*sin(aRad))-d)
		ei_selectsegment(insulatorInnerRadius-0.01,bottom + electrodeThickness -(r - r*sin(aRad))-d)
		ei_setsegmentprop("",meshSizeUltraFine,0,0,1,"elecCond")
		ei_clearselected()
	end
	
	-- Trench Curvature
	ei_selectarcsegment(insulatorInnerRadius+insulatorThickness+s+wallShift+r*cos(aRad), bottom + electrodeThickness)
	ei_selectarcsegment(insulatorInnerRadius-s-wallShift-r*cos(aRad), bottom + electrodeThickness)
	ei_setarcsegmentprop(meshSizeFine,"",0,1,"elecCond")
	ei_clearselected()	
	
else
	-- no o-ring trench 
	ei_selectsegment(0.01,bottom + electrodeThickness)
	ei_selectsegment(insulatorInnerRadius + 0.01,bottom + electrodeThickness)
	ei_selectsegment(insulatorInnerRadius + insulatorThickness + 0.01,bottom + electrodeThickness)
	ei_setsegmentprop("",meshSizeFine,0,0,1,"elecCond")
	ei_clearselected()
end

--upper electrode
--fixed lines curves
ei_selectsegment(78.0911+(insulatorInnerRadius-50.8)-0.01,bottom + 2*electrodeThickness + electrodeSeparation)
ei_setsegmentprop("",meshSizeFine,0,0,1,"zeroCond")
ei_clearselected()

ei_selectarcsegment(80+(insulatorInnerRadius-50.8),bottom + electrodeThickness + electrodeSeparation + 20)
ei_selectarcsegment(80+(insulatorInnerRadius-50.8),bottom + electrodeThickness + electrodeSeparation + 30)
ei_selectarcsegment(120+(insulatorInnerRadius-50.8),bottom + electrodeThickness + electrodeSeparation + 30)
ei_selectarcsegment(120+(insulatorInnerRadius-50.8),bottom + electrodeThickness + electrodeSeparation + 10)
ei_setarcsegmentprop(meshSizeFine,"",0,1,"zeroCond")
ei_clearselected()

if boolOringTrench == 1 then
	-- two flat surfaces on the electrode
	ei_selectsegment(0.01+insulatorInnerRadius+insulatorThickness+s+wallShift+r*cos(aRad), bottom + electrodeThickness + electrodeSeparation)
	ei_selectsegment(insulatorInnerRadius-s-wallShift-r*cos(aRad)-0.01, bottom + electrodeThickness + electrodeSeparation)

	-- O-ring groove
	ei_selectsegment(oRingCenter+(oRingGrooveWidth/2),bottom + electrodeThickness + electrodeSeparation+(r - r*sin(aRad))+d+0.01)
	ei_selectsegment(oRingCenter-(oRingGrooveWidth/2),bottom + electrodeThickness + electrodeSeparation+(r - r*sin(aRad))+d+0.01)
	ei_selectsegment(oRingCenter-(oRingGrooveWidth/2)+0.01,bottom + electrodeThickness + electrodeSeparation+(r - r*sin(aRad))+d+oRingGrooveHeight)
	ei_setsegmentprop("",meshSizeFine,0,0,1,"zeroCond")
	ei_clearselected()
	
	-- line of contact of insulator and the trench
	ei_selectsegment(insulatorInnerRadius + 0.01,bottom + electrodeThickness + electrodeSeparation+(r - r*sin(aRad))+d)
	ei_selectsegment(insulatorInnerRadius + insulatorThickness - 0.01,bottom + electrodeThickness + electrodeSeparation+(r - r*sin(aRad))+d)
	ei_setsegmentprop("",meshSizeUltraFine,0,0,1,"zeroCond")
	ei_clearselected()

	-- trench wall, d
	if d ~= 0 then
		ei_selectsegment(insulatorInnerRadius+insulatorThickness+s+wallShift+0.01,bottom + electrodeThickness + electrodeSeparation + (r - r*sin(aRad)) + 0.01)
		ei_selectsegment(insulatorInnerRadius-s-wallShift-0.01,bottom + electrodeThickness + electrodeSeparation + (r - r*sin(aRad)) + 0.01)	
		ei_setsegmentprop("",meshSizeUltraFine,0,0,1,"zeroCond")
		ei_clearselected()
	end
	
	-- seperation, s
	if s ~= 0 then
		ei_selectsegment(insulatorInnerRadius+insulatorThickness+0.01,bottom + electrodeThickness + electrodeSeparation + (r - r*sin(aRad)) + d)
		ei_selectsegment(insulatorInnerRadius-0.01,bottom + electrodeThickness + electrodeSeparation + (r - r*sin(aRad)) + d)
		ei_setsegmentprop("",meshSizeUltraFine,0,0,1,"zeroCond")
		ei_clearselected()
	end
	
	-- Trench Curvature
	ei_selectarcsegment(insulatorInnerRadius+insulatorThickness+s+wallShift+r*cos(aRad), bottom + electrodeThickness + electrodeSeparation)
	ei_selectarcsegment(insulatorInnerRadius-s-wallShift-r*cos(aRad), bottom + electrodeThickness + electrodeSeparation)
	ei_setarcsegmentprop(meshSizeFine,"",0,1,"zeroCond")
	ei_clearselected()
else
	-- no o-ring trench 
	ei_selectsegment(insulatorInnerRadius - 0.01,bottom + electrodeThickness + electrodeSeparation)
	ei_selectsegment(insulatorInnerRadius + 0.01,bottom + electrodeThickness + electrodeSeparation)
	ei_selectsegment(insulatorInnerRadius + insulatorThickness + 0.01,bottom + electrodeThickness + electrodeSeparation)
	ei_setsegmentprop("",meshSizeFine,0,0,1,"zeroCond")
	ei_clearselected()
end

if boolGasFilling == 1 then
	ei_selectsegment(1.5,bottom + electrodeThickness + electrodeSeparation + 3 + 0.1)
	ei_selectsegment(1.5+0.01,bottom + electrodeThickness + electrodeSeparation + 3 + 0.5367 + 0.1)
	ei_selectsegment(8.9281,bottom + electrodeThickness + electrodeSeparation + 3 + 0.5367 + 4.4633 + 0.1)
	ei_selectsegment(8.9281 + 0.1,bottom + electrodeThickness + electrodeSeparation + 3 + 0.5367 + 4.4633 + 1)
	ei_selectsegment(9.8091,bottom + electrodeThickness + electrodeSeparation + 3 + 0.5367 + 4.4633 + 1 + 11.0049-0.1)
	ei_setsegmentprop("",meshSizeFine,0,0,1,"zeroCond")
	ei_clearselected()
	
	ei_selectarcsegment(4.5,bottom + electrodeThickness + electrodeSeparation)
	ei_setarcsegmentprop(meshSizeFine,"",0,1,"zeroCond")
	ei_clearselected()	
end

-- Localize Er
-- This section, we create two boxes located at the upper and lower
-- trench area. Please refer to the report to see detailed explaination.
if boolOringTrench == 1 then
	ei_addnode(0.75*insulatorInnerRadius,bottom + electrodeThickness)
	ei_addnode(0.75*insulatorInnerRadius,bottom + electrodeThickness + electrodeSeparation / 6)
	ei_addnode(insulatorInnerRadius,bottom + electrodeThickness + electrodeSeparation / 6)
	ei_addsegment(0.75*insulatorInnerRadius,bottom + electrodeThickness,0.75*insulatorInnerRadius,bottom + electrodeThickness + electrodeSeparation / 6)
	ei_addsegment(insulatorInnerRadius,bottom + electrodeThickness + electrodeSeparation / 6,0.75*insulatorInnerRadius,bottom + electrodeThickness + electrodeSeparation / 6)
	ei_addblocklabel((insulatorInnerRadius + 0.75*insulatorInnerRadius) / 2, (bottom + electrodeThickness + bottom + electrodeThickness + electrodeSeparation / 6) / 2)

	ei_addnode(0.75*insulatorInnerRadius,bottom + electrodeThickness + electrodeSeparation)
	ei_addnode(0.75*insulatorInnerRadius,bottom + electrodeThickness + electrodeSeparation *5 / 6)
	ei_addnode(insulatorInnerRadius,bottom + electrodeThickness + electrodeSeparation*5 / 6)
	ei_addsegment(0.75*insulatorInnerRadius,bottom + electrodeThickness + electrodeSeparation,0.75*insulatorInnerRadius,bottom + electrodeThickness + electrodeSeparation *5 / 6)
	ei_addsegment(insulatorInnerRadius,bottom + electrodeThickness + electrodeSeparation *5 / 6,0.75*insulatorInnerRadius,bottom + electrodeThickness + electrodeSeparation *5 / 6)
	ei_addblocklabel((insulatorInnerRadius + 0.75*insulatorInnerRadius) / 2, (bottom + electrodeThickness + electrodeSeparation + bottom + electrodeThickness + electrodeSeparation*5 / 6) / 2)

	ei_selectlabel((insulatorInnerRadius + 0.75*insulatorInnerRadius) / 2, (bottom + electrodeThickness + bottom + electrodeThickness + electrodeSeparation / 6) / 2)
	ei_selectlabel((insulatorInnerRadius + 0.75*insulatorInnerRadius) / 2, (bottom + electrodeThickness + electrodeSeparation + bottom + electrodeThickness + electrodeSeparation*5 / 6) / 2)
	ei_setblockprop("Vacuum",0,meshSizeMedium,0)
	ei_clearselected()
end

-- 6. Start meshing and solving
ei_zoomnatural() --zoom to natural view
ei_saveas("temp.fee")
ei_analyze()
ei_loadsolution()

-- -- 7. postprocessor
eo_shownames(1) --show block label
eo_showdensityplot(1,0,2,3.5e6,0)
eo_showvectorplot(2,1)

-- drawing the contour line. start from the center of the bottom electrode, follow the surface profile
-- and end at the triple point
if s ~= 0 then
	eo_addcontour(0,bottom+electrodeThickness+0.01)
	eo_addcontour(insulatorInnerRadius-s-wallShift-r*cos(aRad), bottom + electrodeThickness+0.01)
	eo_addcontour(insulatorInnerRadius-s-wallShift+0.01,bottom + electrodeThickness -(r - r*sin(aRad))+0.01)
	
	if r~= 0 then
		eo_bendcontour(-90+a,1)
	end

	if d~= 0 then
		eo_addcontour(insulatorInnerRadius-s+0.01,bottom + electrodeThickness -(r - r*sin(aRad))-d+0.01)
	end

	eo_addcontour(insulatorInnerRadius-0.01,bottom + electrodeThickness -(r - r*sin(aRad)) - d+0.01)

	-- fileName = "C:\\\FEMMdata\\Ertest\\".."s"..s.."r"..r.."d"..d..".txt"
	-- eo_makeplot(4,contourPtNum,fileName,1) 
	-- eo_clearcontour()
end

-- evaluate the localized Er and Ez
eo_selectblock((insulatorInnerRadius + 0.75*insulatorInnerRadius) / 2, (bottom + electrodeThickness + bottom + electrodeThickness + electrodeSeparation / 6) / 2)
Erlower, Ezlower = eo_blockintegral(4)
eo_clearblock()

eo_selectblock((insulatorInnerRadius + 0.75*insulatorInnerRadius) / 2, (bottom + electrodeThickness + electrodeSeparation + bottom + electrodeThickness + electrodeSeparation*5 / 6) / 2)
Erupper, Ezupper = eo_blockintegral(4)
eo_clearblock()

-- file = openfile("average E.txt","a")
-- appendto ("average E.txt")
-- write("s","\t",s,"\t","r", "\t", r,"\t", "d","\t", d,"\t","Erlower","\t", Erlower,"\t","Ezlower","\t",Ezlower,"\t","Erupper","\t", Erupper,"\t","Ezupper","\t",Ezupper,"\n")
-- closefile(file)
-- eo_close()

-- messagebox("Done")

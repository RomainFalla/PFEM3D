Problem = {
    id = "Conduction",
	simulationTime = 1,
	verboseOutput = false,
	
	Mesh = {
		hchar = 0.05,
		alpha = 1.2,
		omega = 0.7,
		gamma = 0.7,
		boundingBox = {-1, -0.25, 5, 1.25},
		mshFile = "examples/2D/conduction/geometry.msh"
	},
	
	Extractors = {
		{
			kind = "GMSH",
			outputFile = "results.msh",
			timeBetweenWriting = 0.01,
			whatToWrite = {"T", "normals"},
			writeAs = "NodesElements" 
		}
	},
	
	Material = {
		rho = 1,
		k = 237,
		cv = 1000
	},
	
	IC = {
		UpFixed = true,
		DownFixed = true,
		RightFixed = true,
		LeftFixed = true,
	},
	
	Solver = {
	    id = "PSPG",
		adaptDT = true,
		coeffDTincrease = 1.5,
		coeffDTDecrease = 2,
		maxDT = 0.1,
		initialDT = 0.025,
		
		HeatEq = {
			minRes = 1e-6,
			maxIter = 10,
			BC = {

			}
		}
	}
}

function Problem.IC:initStates(pos)
	return {10*math.exp(-((pos[1]-0.5)^2 + (pos[2]-0.5)^2)/0.05)}
end

function Problem.Solver.HeatEq.BC:DownQ(pos, initPos, states, t)
	return {0, 0}
end

function Problem.Solver.HeatEq.BC:UpQ(pos, initPos, states, t)
	return {0, 0}
end

function Problem.Solver.HeatEq.BC:LeftT(pos, initPos, states, t)
	return {0}
end

function Problem.Solver.HeatEq.BC:RightT(pos, initPos, states, t)
	return {0}
end

-- function Problem.Solver.HeatEq.BC:DownQ(pos, initPos, states, t)
	-- return {0}
-- end

-- function Problem.Solver.HeatEq.BC:RightQ(pos, initPos, states, t)
	-- return {0}
-- end
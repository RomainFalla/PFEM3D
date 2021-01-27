Problem = {
    id = "Boussinesq",
	simulationTime = 50,
	verboseOutput = false,
	
	Mesh = {
		hchar = 0.02,
		alpha = 1.2,
		omega = 0.7,
		gamma = 0.7,
		boundingBox = {-0.05, -0.05, 4.05, 100},
		mshFile = "examples/2D/rayleigh/geometry.msh"
	},
	
	Extractors = {
		{
			kind = "GMSH",
			outputFile = "results.msh",
			timeBetweenWriting = 0.4,
			whatToWrite = {"T", "ke"},
			writeAs = "NodesElements" 
		}
	},
	
	Material = {
		mu = 1e-3,
		rho = 1000,
		k = 0.6,
		alpha = 69e-6,
		Tr = 650,
		cv = 1,
		gamma = 0
	},
	
	IC = {
		TopFixed = true,
		BottomFixed = true,
		LeftFixed = true,
		RightFixed = true
	},
	
	Solver = {
	    id = "PSPG",
		adaptDT = true,
		coeffDTincrease = 1.5,
		coeffDTDecrease = 2,
		maxDT = 0.1,
		initialDT = 0.1,
		solveHeatFirst = true,
		
		MomContEq = {
			minRes = 1e-6,
			maxIter = 10,
			bodyForce = {0, -9.81},
			BC = {

			}
		},
		
		HeatEq = {
			minRes = 1e-6,
			maxIter = 10,
			BC = {

			}
		}
	}
}

function Problem.IC:initStates(pos)
	return {0, 0, 0, 650}
end

function Problem.IC:initTopStates(pos)
	return {0, 0, 0, 300}
end

function Problem.IC:initBottomStates(pos)
	return {0, 0, 0, 1000}
end

function Problem.Solver.HeatEq.BC:TopT(pos, initPos, states, t) 
	return {300}
end

function Problem.Solver.MomContEq.BC:TopV(pos, initPos, states, t) 
	return {0, 0}
end

function Problem.Solver.HeatEq.BC:BottomT(pos, initPos, states, t) 
	return {1000}
end

function Problem.Solver.MomContEq.BC:BottomV(pos, initPos, states, t) 
	return {0, 0}
end

function Problem.Solver.HeatEq.BC:LeftQ(pos, initPos, states, t) 
	return {0}
end

function Problem.Solver.MomContEq.BC:LeftV(pos, initPos, states, t) 
	return {0, 0}
end

function Problem.Solver.HeatEq.BC:RightQ(pos, initPos, states, t) 
	return {0}
end

function Problem.Solver.MomContEq.BC:RightV(pos, initPos, states, t) 
	return {0, 0}
end

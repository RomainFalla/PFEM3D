Problem = {
    id = "Boussinesq",
	simulationTime = 6,
	verboseOutput = false,
	
	Mesh = {
		hchar = 0.02,
		alpha = 1.2,
		omega = 0.7,
		gamma = 0.7,
		boundingBox = {-0.05, -0.05, 4.05, 100},
		mshFile = "examples/2D/coolingBallInFluid/geometry.msh"
	},
	
	Extractors = {
		{
			kind = "GMSH",
			outputFile = "results.msh",
			timeBetweenWriting = 0.05,
			whatToWrite = {"T", "ke", "p"},
			writeAs = "NodesElements" 
		}
	},
	
	Material = {
		mu = 1e-3,
		rho = 1000,
		k = 0.6,
		alpha = 69e-6,
		Tr = 650,
		cv = 4186,
		gamma = 0
	},
	
	IC = {
		BoundaryFixed = true
	},
	
	Solver = {
	    id = "PSPG",
		adaptDT = true,
		coeffDTincrease = 1.5,
		coeffDTDecrease = 2,
		maxDT = 0.01,
		initialDT = 0.01,
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
	local T = 300
	if(pos[2] > 1 and pos[1] > 0 and pos[1] < 4) then
        T = 1000
    end
        
    return {0, 0, 0, T};
end

function Problem.Solver.HeatEq.BC:BoundaryQ(pos, initPos, states, t) 
	return {0}
end

function Problem.Solver.MomContEq.BC:BoundaryV(pos, initPos, states, t) 
	return {0, 0}
end

Problem = {
    id = "IncompNewtonNoT",
	simulationTime = 4,
	verboseOutput = false,
	
	Mesh = {
		hchar = 0.02,
		alpha = 1.2,
		omega = 0.7,
		gamma = 0.7,
		boundingBox = {-2, -1, 12, 100},
		mshFile = "examples/2D/sloshing/geometry.msh"
	},
	
	Extractors = {
		{
			kind = "MinMax",
			outputFile = "topPositionSloshing.txt",
			timeBetweenWriting = 0.01,
			minMax = "max",
			coordinate = 1 
		},
		{
			kind = "GMSH",
			outputFile = "results.msh",
			timeBetweenWriting = 0.05,
			whatToWrite = {"p", "ke"},
			writeAs = "NodesElements" 
		}
	},
	
	Material = {
		mu = 0.01,
		rho = 1,
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
		
		MomContEq = {
			minRes = 1e-6,
			maxIter = 10,
			bodyForce = {0, -9.81},
			BC = {

			}
		}
	}
}

function Problem.IC:initStates(pos)
	return {0, 0, 0}
end

function Problem.Solver.MomContEq.BC:BoundaryV(pos, initPos, states, t)
	return {0, 0}
end
Problem = {
    id = "IncompNewtonNoT",
	simulationTime = 4,
	verboseOutput = false,
	
	Mesh = {
		hchar = 0.025,
		alpha = 1.2,
		omega = 0.7,
		gamma = 0.7,
		boundingBox = {-1, -1, -1, 1, 1, 100},
		mshFile = "examples/3D/dropFallInFluid/geometry.msh"
	},
	
	Extractors = {
		{
			kind = "GMSH",
			outputFile = "results.msh",
			timeBetweenWriting = 0.01,
			whatToWrite = {"p", "ke"},
			writeAs = "NodesElements" 
		}
	},
	
	Material = {
		mu = 1e-3,
		rho = 1000,
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
		maxDT = 0.001,
		initialDT = 0.001,
		
		MomContEq = {
			minRes = 1e-6,
			maxIter = 10,
			bodyForce = {0, 0, -9.81},
			BC = {

			}
		}
	}
}

function Problem.IC:initStates(pos)
	return {0, 0, 0, 0}
end

function Problem.Solver.MomContEq.BC:BoundaryV(pos, initPos, states, t)
	return {0, 0, 0}
end
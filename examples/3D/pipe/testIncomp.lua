Problem = {
    id = "IncompNewtonNoT",
	simulationTime = 4,
	verboseOutput = true,
	
	Mesh = {
		hchar = 0.1,
		alpha = 1.2,
		omega = 0.7,
		gamma = 0.7,
		boundingBox = {-0.1, -0.1, 5, 1.1},
		mshFile = "examples/2D/pipe/geometry.msh"
	},
	
	Extractors = {
		{
			kind = "GMSH",
			outputFile = "results.msh",
			timeBetweenWriting = 0.1,
			whatToWrite = {"p", "u"},
			writeAs = "NodesElements" 
		}
	},
	
	Material = {
		mu = 200,
		rho = 1000,
		gamma = 0
	},
	
	IC = {
		BoundaryFixed = false,
		FluidInputFixed = true
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
			bodyForce = {0, 0},
			BC = {

			}
		}
	}
}

function Problem.IC:initStates(pos)
	return {0, 0, 0}
end

function Problem.IC:initFluidInputStates(pos)
	return {1, 0, 0}
end

function Problem.Solver.MomContEq.BC:BoundaryV(pos, initPos, states, t)
	return {0, 0}
end

function Problem.Solver.MomContEq.BC:FluidInputV(pos, initPos, states, t)
	return {1, 0}
end
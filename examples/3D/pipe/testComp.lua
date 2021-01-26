Problem = {
    id = "WCompNewtonNoT",
	simulationTime = 4,
	verboseOutput = true,
	
	Mesh = {
		hchar = 0.05,
		alpha = 1.2,
		omega = 0.35,
		gamma = 0.7,
		boundingBox = {-1, -1, -0.5, 5, 1, 1.5},
		mshFile = "examples/3D/pipe/geometry.msh"
	},
	
	Extractors = {
		{
			kind = "Point",
			outputFile = "line_h_0.1.txt",
			timeBetweenWriting = 0.1,
			whatToWrite = "u",
			points = {{3, -0.4, 0.5}, {3, -0.3, 0.5}, {3, -0.2, 0.5}, {3, -0.1, 0.5}, {3, 0, 0.5}, {3, 0.1, 0.5}, {3, 0.2, 0.5}, {3, 0.3, 0.5}, {3, 0.4, 0.5}}
		},
		{
			kind = "GMSH",
			outputFile = "results.msh",
			timeBetweenWriting = 0.05,
			whatToWrite = {"p", "u"},
			writeAs = "NodesElements" 
		}
	},
	
	Material = {
		mu = 200,
		gamma = 0,
		K0 = 2200000,
		K0p = 7.6,
		rhoStar = 1000
	},
	
	IC = {
		BoundaryFixed = false,
		FluidInputFixed = true
	},
	
	Solver = {
	    id = "CDS",
		adaptDT = true,
		maxDT = 0.001,
		initialDT = 1e-8,
		securityCoeff = 0.05,
		
		MomEq = {
			bodyForce = {0, 0, 0},
			BC = {
			
			}
		},
		
		ContEq = {
			strongContinuity = false,
			BC = {

			}
		}
	}
}

function Problem.IC:initStates(pos)
    local rhoStar = Problem.Material.rhoStar 
	return {1, 0, 0, 0, rhoStar, 0, 0, 0}
end

function Problem.IC:initBoundaryStates(pos)
	local rhoStar = Problem.Material.rhoStar
	return {0, 0, 0, 0, rhoStar, 0, 0, 0}
end

function Problem.Solver.MomEq.BC:BoundaryV(pos, initPos, states, t)
	return {0, 0, 0}
end

function Problem.Solver.MomEq.BC:FluidInputV(pos, initPos, states, t)
	return {0, 0, 0}
end
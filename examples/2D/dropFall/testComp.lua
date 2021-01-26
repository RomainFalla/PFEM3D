Problem = {
    id = "WCompNewtonNoT",
	simulationTime = 4,
	verboseOutput = false,
	
	Mesh = {
		hchar = 0.0292,
		alpha = 1.2,
		omega = 0.7,
		gamma = 0.7,
		boundingBox = {-2, -1, 12, 100},
		mshFile = "examples/2D/dropFall/geometry.msh"
	},
	
	Extractors = {
		{
			kind = "GMSH",
			outputFile = "results.msh",
			timeBetweenWriting = 0.2,
			whatToWrite = {"p", "ke"},
			writeAs = "NodesElements" 
		}
	},
	
	Material = {
		mu = 1e-3,
		gamma = 0,
		K0 = 2200000,
		K0p = 7.6,
		rhoStar = 1000
	},
	
	IC = {
		BoundaryFixed = true,
	},
	
	Solver = {
	    id = "CDS",
		adaptDT = true,
		maxDT = 0.001,
		initialDT = 1e-8,
		securityCoeff = 1,
		
		MomEq = {
			bodyForce = {0, -9.81},
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
	return {0, 0, 0, Problem.Material.rhoStar, 0, 0}
end

function Problem.Solver.MomEq.BC:BoundaryV(pos, initPos, states, t)
	return {0, 0}
end

Problem = {
    id = "WCompNewtonNoT",
	simulationTime = 4,
	verboseOutput = false,
	
	Mesh = {
		hchar = 0.0073,
		alpha = 1.2,
		omega = 0.7,
		gamma = 0.7,
		boundingBox = {-0.05, -0.05, 0.634, 100},
		mshFile = "examples/2D/dropFallInFluid/geometry.msh"
	},
	
	Extractors = {
		{
			kind = "GMSH",
			outputFile = "results.msh",
			timeBetweenWriting = 0.05,
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
	local K0 = Problem.Material.K0
	local K0p = Problem.Material.K0p
	local rhoStar = Problem.Material.rhoStar
	local g = -Problem.Solver.MomEq.bodyForce[2]

	if(pos[2] <= 0.146) then
		local rho, p
		rho = rhoStar*((K0p - 1)/K0*rhoStar*g*(0.146 - pos[2]) + 1)^(1/(K0p - 1))
		p = K0/K0p*((rho/rhoStar)^K0p - 1)
		
		return {0, 0, p, rho, 0, 0}
	else
		return {0, 0, 0, rhoStar, 0, 0}
	end
end

function Problem.Solver.MomEq.BC:BoundaryV(pos, initPos, states, t)
	return {0, 0}
end

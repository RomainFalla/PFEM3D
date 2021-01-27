Problem = {
    id = "WCompNewtonNoT",
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
		gamma = 0,
		K0 = 2200000,
		K0p = 7.6,
		rhoStar = 1
	},
	
	IC = {
		BoundaryFixed = true,
	},
	
	Solver = {
	    id = "CDS",
		adaptDT = true,
		maxDT = 0.1,
		initialDT = 1e-8,
		securityCoeff = 2000,
		
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
    local rhoStar = Problem.Material.rhoStar
    local K0 = Problem.Material.K0
	local K0p = Problem.Material.K0p
	local g = -Problem.Solver.MomEq.bodyForce[2]
	local z0 = 1 + 0.1*math.sin(math.pi*pos[1]/1 - math.pi/2);
	
	if(pos[2] <= z0) then
	    local rho = rhoStar*((K0p - 1)/K0*rhoStar*g*(z0 - pos[2]) + 1)^(1/(K0p - 1))
		local p = K0/K0p*((rho/rhoStar)^K0p - 1)
		return {0, 0, p, rho, 0, 0}
	else 
		return {0, 0, 0, rhoStar, 0, 0}
	end
end

function Problem.Solver.MomEq.BC:BoundaryV(pos, initPos, states, t)
	return {0, 0}
end

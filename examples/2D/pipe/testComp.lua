Problem = {
    id = "WCompNewtonNoT",
	simulationTime = 4,
	verboseOutput = false,
	
	Mesh = {
		hchar = 0.1,
		alpha = 1.1,
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
		gamma = 0,
		K0 = 2200000,
		K0p = 7.6,
		rhoStar = 1000
	},
	
	IC = {
		BoundaryFixed = true,
		FluidInputFixed = true,
	},
	
	Solver = {
	    id = "CDS",
		adaptDT = true,
		maxDT = 0.005,
		initialDT = 1e-8,
		securityCoeff = 10,
		
		MomEq = {
			bodyForce = {0, 0},
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
    --[[local K0 = Problem.Material.K0
	local K0p = Problem.Material.K0p
	local g = -Problem.Solver.MomEq.bodyForce[2]
	local z0 = 2*0.146
	
	if(pos[2] <= z0 and pos[1] <= z0/2 + 1.1*Problem.Mesh.hchar) then
	    local rho = rhoStar*((K0p - 1)/K0*rhoStar*g*(z0 - pos[2]) + 1)^(1/(K0p - 1))
		local p = K0/K0p*((rho/rhoStar)^K0p - 1)
		return {0, 0, p, rho, 0, 0}
	else --]] 
	return {1, 0, 0, rhoStar, 0, 0}
    --end
end

function Problem.IC:initBoundaryStates(pos)
	local rhoStar = Problem.Material.rhoStar
	return {0, 0, 0, rhoStar, 0, 0}
end

function Problem.Solver.MomEq.BC:BoundaryV(pos, initPos, states, t)
	return {0, 0}
end

function Problem.Solver.MomEq.BC:FluidInputV(pos, initPos, states, t)
	return {0, 0}
end
Problem = {
    id = "BoussinesqWC",
	simulationTime = 4,
	verboseOutput = true,
	
	Mesh = {
		hchar = 0.1,
		alpha = 1.2,
		omega = 0.7,
		gamma = 0.7,
		boundingBox = {-0.05, -0.05, -0.05, 4.05, 4.05, 100},
		mshFile = "examples/3D/rayleigh/geometry.msh"
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
		gamma = 0,
		K0 = 2200000,
		K0p = 7.6,
		rhoStar = 1000,
		k = 0.6,
		cv = 1,
		alpha = 69e-6,
		Tr = 650
	},
	
	IC = {
		LateralFixed = true,
		TopFixed = true,
		BottomFixed = true
	},
	
	Solver = {
	    id = "CDS",
		adaptDT = true,
		maxDT = 0.01,
		initialDT = 1e-8,
		securityCoeff = 0.5,
		
		MomEq = {
			bodyForce = {0, 0, -9.81},
			BC = {
			
			}
		},
		
		ContEq = {
			strongContinuity = false,
			BC = {

			}
		},
		
		HeatEq = {
			BC = {
			
			}
		}
	}
}

function Problem.IC:initStates(pos)
    local rhoStar = Problem.Material.rhoStar
    local K0 = Problem.Material.K0
	local K0p = Problem.Material.K0p
	local g = -Problem.Solver.MomEq.bodyForce[3]
	local z0 = 1
	
	local rho = rhoStar*((K0p - 1)/K0*rhoStar*g*(z0 - pos[3]) + 1)^(1/(K0p - 1))
	local p = K0/K0p*((rho/rhoStar)^K0p - 1)
	return {0, 0, 0, p, rho, 0, 0, 0, 650}
end

function Problem.IC:initTopStates(pos)
	local rhoStar = Problem.Material.rhoStar
    local K0 = Problem.Material.K0
	local K0p = Problem.Material.K0p
	local g = -Problem.Solver.MomEq.bodyForce[3]
	local z0 = 1
	
	local rho = rhoStar*((K0p - 1)/K0*rhoStar*g*(z0 - pos[3]) + 1)^(1/(K0p - 1))
	local p = K0/K0p*((rho/rhoStar)^K0p - 1)
	return {0, 0, 0, p, rho, 0, 0, 0, 300}
end

function Problem.IC:initBottomStates(pos)
	local rhoStar = Problem.Material.rhoStar
    local K0 = Problem.Material.K0
	local K0p = Problem.Material.K0p
	local g = -Problem.Solver.MomEq.bodyForce[3]
	local z0 = 1
	
	local rho = rhoStar*((K0p - 1)/K0*rhoStar*g*(z0 - pos[3]) + 1)^(1/(K0p - 1))
	local p = K0/K0p*((rho/rhoStar)^K0p - 1)
	return {0, 0, 0, p, rho, 0, 0, 0, 1000}
end

function Problem.Solver.MomEq.BC:LateralV(pos, initPos, states, t)
	return {0, 0, 0}
end

function Problem.Solver.MomEq.BC:TopV(pos, initPos, states, t)
	return {0, 0, 0}
end

function Problem.Solver.MomEq.BC:BottomV(pos, initPos, states, t)
	return {0, 0, 0}
end

function Problem.Solver.MomEq.BC:LateralQ(pos, initPos, states, t)
	return {0}
end

function Problem.Solver.MomEq.BC:TopT(pos, initPos, states, t)
	return {300}
end

function Problem.Solver.MomEq.BC:BottomT(pos, initPos, states, t)
	return {1000}
end

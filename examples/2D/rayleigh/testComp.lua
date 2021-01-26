Problem = {
    id = "BoussinesqWC",
	simulationTime = 50,
	verboseOutput = false,
	
	Mesh = {
		hchar = 0.02,
		alpha = 1.2,
		omega = 0.7,
		gamma = 0.7,
		boundingBox = {-0.05, -0.05, 4.05, 100},
		mshFile = "examples/2D/rayleigh/geometry.msh"
	},
	
	Extractors = {
		{
			kind = "GMSH",
			outputFile = "results.msh",
			timeBetweenWriting = 0.4,
			whatToWrite = {"T", "ke", "p"},
			writeAs = "NodesElements" 
		}
	},
	
	Material = {
		mu = 1e-3,
		k = 0.6,
		alpha = 69e-6,
		Tr = 650,
		cv = 1,
		gamma = 0,
		K0 = 2200000,
		K0p = 7.6,
		rhoStar = 1000
	},
	
	IC = {
		TopFixed = true,
		BottomFixed = true,
		LeftFixed = true,
		RightFixed = true
	},
	
	Solver = {
	    id = "CDS",
		adaptDT = true,
		maxDT = 0.005,
		initialDT = 1e-8,
		securityCoeff = 10,
		
		MomEq = {
			bodyForce = {0, -9.81},
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
	local K0 = Problem.Material.K0
	local K0p = Problem.Material.K0p
	local rhoStar = Problem.Material.rhoStar
	local g = -Problem.Solver.MomEq.bodyForce[2]

	local rho, p
	rho = rhoStar*((K0p - 1)/K0*rhoStar*g*(1 - pos[2]) + 1)^(1/(K0p - 1))
	p = K0/K0p*((rho/rhoStar)^K0p - 1)

	return {0, 0, p, rho, 0, 0, 1150}
end

function Problem.IC:initTopStates(pos)
    local K0 = Problem.Material.K0
	local K0p = Problem.Material.K0p
	local rhoStar = Problem.Material.rhoStar
	local g = -Problem.Solver.MomEq.bodyForce[2]

	local rho, p
	local rho, p
	rho = rhoStar*((K0p - 1)/K0*rhoStar*g*(1 - pos[2]) + 1)^(1/(K0p - 1))
	p = K0/K0p*((rho/rhoStar)^K0p - 1)
	
	return {0, 0, p, rho, 0, 0, 300}
end

function Problem.IC:initBottomStates(pos)
	local K0 = Problem.Material.K0
	local K0p = Problem.Material.K0p
	local rhoStar = Problem.Material.rhoStar
	local g = -Problem.Solver.MomEq.bodyForce[2]

	local rho, p
	rho = rhoStar*((K0p - 1)/K0*rhoStar*g*(1 - pos[2]) + 1)^(1/(K0p - 1))
	p = K0/K0p*((rho/rhoStar)^K0p - 1)
	
	return {0, 0, p, rho, 0, 0, 2000}
end

function Problem.Solver.HeatEq.BC:TopT(pos, initPos, states, t) 
	return {300}
end

function Problem.Solver.MomEq.BC:TopV(pos, initPos, states, t) 
	return {0, 0}
end

function Problem.Solver.HeatEq.BC:BottomT(pos, initPos, states, t) 
	return {2000}
end

function Problem.Solver.MomEq.BC:BottomV(pos, initPos, states, t) 
	return {0, 0}
end

function Problem.Solver.HeatEq.BC:LeftQ(pos, initPos, states, t) 
	return {0}
end

function Problem.Solver.MomEq.BC:LeftV(pos, initPos, states, t) 
	return {0, 0}
end

function Problem.Solver.HeatEq.BC:RightQ(pos, initPos, states, t) 
	return {0}
end

function Problem.Solver.MomEq.BC:RightV(pos, initPos, states, t) 
	return {0, 0}
end


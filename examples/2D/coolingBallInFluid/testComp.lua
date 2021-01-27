Problem = {
    id = "BoussinesqWC",
	simulationTime = 6,
	verboseOutput = false,
	
	Mesh = {
		hchar = 0.02,
		alpha = 1.2,
		omega = 0.7,
		gamma = 0.7,
		boundingBox = {-0.05, -0.05, 4.05, 100},
		mshFile = "examples/2D/coolingBallInFluid/geometry.msh"
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
		k = 0.6,
		alpha = 69e-6,
		Tr = 650,
		cv = 4186,
		gamma = 0,
		K0 = 2200000,
		K0p = 7.6,
		rhoStar = 1000
	},
	
	IC = {
		BoundaryFixed = true
	},
	
	Solver = {
	    id = "CDS",
		adaptDT = true,
		maxDT = 0.0025,
		initialDT = 1e-8,
		securityCoeff = 0.1,
		
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

function Problem.IC:initStates(pos)
	local T = 300
	if(pos[2] > 1 and pos[1] > 0 and pos[1] < 4) then
        T = 1000
    end
	
	if(pos[2] <= 1) then
		local K0 = Problem.Material.K0
		local K0p = Problem.Material.K0p
		local rhoStar = Problem.Material.rhoStar
		local g = -Problem.Solver.MomEq.bodyForce[2]

		local rho, p
		local rho, p
		rho = rhoStar*((K0p - 1)/K0*rhoStar*g*(1 - pos[2]) + 1)^(1/(K0p - 1))
		p = K0/K0p*((rho/rhoStar)^K0p - 1)
		
		return {0, 0, p, rho, 0, 0, T}
	else
		local rhoStar = Problem.Material.rhoStar
		return {0, 0, 0, rhoStar, 0, 0, T}
	end
end

function Problem.Solver.HeatEq.BC:BoundaryQ(pos, initPos, states, t) 
	return {0}
end

function Problem.Solver.MomEq.BC:BoundaryV(pos, initPos, states, t) 
	return {0, 0}
end
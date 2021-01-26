Problem = {
    id = "WCompNewtonNoT",
	simulationTime = 1.2,
	verboseOutput = false,
	
	Mesh = {
		hchar = 0.25,
		alpha = 1.3,
		omega = 0.7,
		gamma = 0.7,
		boundingBox = {-2, -1, 12, 100},
		mshFile = "examples/2D/hydrostatic/geometry.msh"
	},
	
	Extractors = {
		{
			kind = "GMSH",
			outputFile = "results.msh",
			timeBetweenWriting = 0.2,
			whatToWrite = {"p", "ke"},
			writeAs = "NodesElements" 
		},
		{
			kind = "Point",
			outputFile = "pLine_vert.txt",
			timeBetweenWriting = 0.1,
			whatToWrite = "p",
			points = {{5, 0.1}, {5, 0.5}, {5, 1}, {5, 1.5}, {5, 2}, {5, 2.5}, {5, 3}, {5, 3.5}, {5, 4}, {5, 4.5}, {5, 4.9}} 
		},
		{
			kind = "Point",
			outputFile = "pLine_horiz.txt",
			timeBetweenWriting = 0.1,
			whatToWrite = "p",
			points = {{0.1, 2.5}, {1, 2.5}, {2, 2.5}, {3, 2.5}, {4, 2.5}, {5, 2.5}, {6, 2.5}, {7, 2.5}, {8, 2.5}, {9, 2.5},{9.9, 2.5}} 
		},
		{
			kind = "Mass",
			outputFile = "mass.txt",
			timeBetweenWriting = 0.1
		}
	},
	
	Material = {
		mu = 1e-3,
		gamma = 0,
		K0 = 2200,
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
		securityCoeff = 0.5,
		
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

	local rho, p
	if(pos[2] <= 5) then
		rho = rhoStar*((K0p - 1)/K0*rhoStar*g*(5 - pos[2]) + 1)^(1/(K0p - 1))
	    p = K0/K0p*((rho/rhoStar)^K0p - 1)
	else
		rho = rhoStar
		p = 0
	end
	return {0, 0, p, rho, 0, 0}
end


function Problem.Solver.MomEq.BC:BoundaryV(pos, initPos, states, t)
	return {0, 0}
end

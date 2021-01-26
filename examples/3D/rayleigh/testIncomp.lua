Problem = {
    id = "IncompNewtonNoT",
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
		rho = 1000,
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
	    id = "PSPG",
		adaptDT = true,
		coeffDTincrease = 1.5,
		coeffDTDecrease = 2,
		maxDT = 0.001,
		initialDT = 0.001,
		
		MomContEq = {
			minRes = 1e-6,
			maxIter = 10,
			bodyForce = {0, 0, -9.81},
			BC = {

			}
		}
		
		HeatEq = {
			minRes = 1e-6,
			maxIter = 10,
			BC = {
			
			}
		}
	}
}

function Problem.IC:initStates(pos)
    return {0, 0, 0, 0, 650}
end

function Problem.IC:initTopStates(pos)
	return {0, 0, 0, 0, 300}
end

function Problem.IC:initBottomStates(pos)
	return {0, 0, 0, 0, 1000}
end

function Problem.Solver.MomContEq.BC:LateralV(pos, initPos, states, t)
	return {0, 0, 0}
end

function Problem.Solver.MomContEq.BC:TopV(pos, initPos, states, t)
	return {0, 0, 0}
end

function Problem.Solver.MomContEq.BC:BottomV(pos, initPos, states, t)
	return {0, 0, 0}
end

function Problem.Solver.MomContEq.BC:LateralQ(pos, initPos, states, t)
	return {0}
end

function Problem.Solver.MomContEq.BC:TopT(pos, initPos, states, t)
	return {300}
end

function Problem.Solver.MomContEq.BC:BottomT(pos, initPos, states, t)
	return {1000}
end

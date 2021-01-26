Problem = {
    id = "IncompNewtonNoT",
	simulationTime = 10.1,
	verboseOutput = true,
	
	Mesh = {
		hchar = 0.25,
		alpha = 1.2,
		omega = 0.7,
		gamma = 0.7,
		boundingBox = {-2, -2, -2, 2, 2, 2},
		mshFile = "examples/3D/cubeToSphere/geometry.msh"
	},
	
	Extractors = {
		{
			kind = "GMSH",
			outputFile = "results.msh",
			timeBetweenWriting = 0.1,
			whatToWrite = {"p", "ke", "normals", "curvatures"},
			writeAs = "NodesElements" 
		} 
	},
	
	Material = {
		mu = 1.0,
		rho = 100,
		gamma = 1.0
	},
	
	IC = {

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
			bodyForce = {0, 0, 0},
			BC = {

			}
		}
	}
}

function Problem.IC:initStates(pos)
	return {0, 0, 0, 0}
end

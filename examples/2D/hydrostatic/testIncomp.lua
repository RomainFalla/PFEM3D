Problem = {
    id = "IncompNewtonNoT",
	simulationTime = 4,
	verboseOutput = false,
	
	Mesh = {
		hchar = 0.25,
		alpha = 1.2,
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
		rho = 1000,
		gamma = 0
	},
	
	IC = {
		BoundaryFixed = true
	},
	
	Solver = {
	    id = "PSPG",
		adaptDT = true,
		coeffDTincrease = 1.5,
		coeffDTDecrease = 2,
		maxDT = 0.1,
		initialDT = 0.1,
		
		MomContEq = {
			minRes = 1e-6,
			maxIter = 10,
			bodyForce = {0, -9.81},
			BC = {

			}
		}
	}
}

function Problem.IC:initStates(pos)
	return {0, 0, 0}
end


function Problem.Solver.MomContEq.BC:BoundaryV(pos, initPos, states, t)
	return {0, 0}
end

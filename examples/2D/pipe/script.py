from pfemSolverw import MeshCreateInfo
from pfemSolverw import SolverCreateInfo
from pfemSolverw import SolverIncompCreateInfo
from pfemSolverw import SolverIncompressible
from pfemSolverw import VectorDouble
from pfemSolverw import VectorString
from pfemSolverw import VectorVectorDouble


meshInfos = MeshCreateInfo()
meshInfos.hchar = 0.1
meshInfos.alpha = 1.2
meshInfos.gamma = 0.7
meshInfos.omega = 0.7
boundingBox = VectorDouble(4)
boundingBox[0] = -1
boundingBox[1] = -1
boundingBox[2] = 5
boundingBox[3] = 2
meshInfos.boundingBox = boundingBox
meshInfos.mshFile = "../../examples/2D/pipe/geometry.msh"

solverInfos = SolverCreateInfo()
solverInfos.gravity = 0
solverInfos.strongPatFS = True
solverInfos.adaptDT = True
solverInfos.maxDT = 0.01
solverInfos.initialDT = 0.01
solverInfos.endTime = 6
solverInfos.IBCfile = "../../examples/2D/pipe/IBC_Incomp.lua"
solverInfos.meshInfos = meshInfos

solverIncompInfos = SolverIncompCreateInfo()
solverIncompInfos.rho = 1000
solverIncompInfos.mu = 1e-3
solverIncompInfos.picardRelTol = 1e-6
solverIncompInfos.picardMaxIter = 10
solverIncompInfos.coeffDTincrease = 1.5
solverIncompInfos.coeffDTdecrease = 2
solverIncompInfos.solverInfos = solverInfos

solver = SolverIncompressible(solverIncompInfos)

whatToWrite = VectorString(2)
whatToWrite[0] = "u"
whatToWrite[1] = "p"
solver.addGMSHExtractor("results.msh", 0.1, whatToWrite,"NodesElements")
solver.solveProblem(False)


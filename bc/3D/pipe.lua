function initFluid(pos) 
	return {1, 0, 0, 0}, false 
end

function initBoundary(pos)
	return {0, 0, 0, 0}, false
end

function initFluidInput(pos)
	return {1, 0, 0, 0}, true
end

function Boundary(pos, initPos, t) 
	return {initPos[1], initPos[2], initPos[3]}
end

function FluidInput(pos, initPos, t) 
	return {1, 0, 0}
end
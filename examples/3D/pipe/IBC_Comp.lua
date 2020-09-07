function initFluid(pos) 
	return {1, 0, 0, 0, 1000, 0, 0, 0}, false 
end

function initBoundary(pos)
	return {0, 0, 0, 0, 1000, 0, 0, 0}, false
end

function initFluidInput(pos)
	return {1, 0, 0, 0, 1000, 0, 0, 0}, true
end

function Boundary(pos, initPos, t) 
	return {0, 0, 0}
end

function FluidInput(pos, initPos, t) 
	return {0, 0, 0}
end

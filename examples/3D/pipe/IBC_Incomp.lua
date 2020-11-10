BoundaryFixed = true
FluidInputFixed = true

function initStates(pos) 
	return {1, 0, 0, 0} 
end

function initBoundaryStates(pos)
	return {0, 0, 0, 0}
end

function Boundary(pos, initPos, t) 
	return {0, 0, 0}
end

function FluidInput(pos, initPos, t) 
	return {1, 0, 0}
end

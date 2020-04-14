function initFluid(pos) 
	return {0, 0, 0, 1000, 0, 0}, false 
end

function initBoundary(pos)
	return {0, 0, 0, 1000, 0, 0}, false
end

function Boundary(pos, initPos, t) 
	return {initPos[1], initPos[2]}
end
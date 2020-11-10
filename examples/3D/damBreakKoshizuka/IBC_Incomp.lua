BoundaryFixed = true

function initStates(pos) 
	return {0, 0, 0, 0}, false 
end

function Boundary(pos, initPos, t) 
	return {0, 0, 0}
end

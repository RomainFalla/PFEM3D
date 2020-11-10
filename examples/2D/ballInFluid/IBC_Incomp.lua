BoundaryFixed = true
DiskBoundaryFixed = false

function initStates(pos) 
	return {0, 0, 0} 
end

function initDiskBoundaryStates(pos) 
	return {1, 0, 0}	
end

function Boundary(pos, initPos, t) 
	return {0, 0}
end

function DiskBoundary(pos, initPos, t)
	return {1, 0}
end

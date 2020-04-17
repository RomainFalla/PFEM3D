function initFluid(pos) 
	return {0, 0, 0}, false 
end

function initBoundary(pos)
	return {0, 0, 0}, false
end

function initDiskBoundary(pos) 
	return {0, 0, 0}, false	
end

function Boundary(pos, initPos, t) 
	return {0, 0}
end

function DiskBoundary(pos, initPos, t)
	return {0, 0}
end
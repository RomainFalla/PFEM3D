function initFluid(pos) 
	return {0, 0, 0}, false 
end

function initMovingBoundary(pos)
	return {0, 0, 0}, false
end

function initStaticBoundary(pos) 
	return {0, 0, 0}, false	
end

function MovingBoundary(pos, initPos, t)
	phi = math.acos(initPos[1])
	if(initPos[2] < 0) then
		phi = - phi
	end
		
	return {math.cos(0.5*t + phi), math.sin(0.5*t + phi)}
end

function StaticBoundary(pos, initPos, t) 
	return {pos[1], pos[2]}
end
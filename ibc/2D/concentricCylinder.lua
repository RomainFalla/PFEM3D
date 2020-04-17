function initFluid(pos) 
	return {0, 0, 0}, false 
end

function initMovingBoundary(pos)
	phi = math.acos(initPos[1])
	if(initPos[2] < 0) then
		phi = - phi
	end
		
	return {-0.5*math.sin(0.5*t + phi), 0.5*math.cos(0.5*t + phi), 0}
end

function initStaticBoundary(pos) 
	return {0, 0, 0}, false	
end

function MovingBoundary(pos, initPos, t)
	phi = math.acos(initPos[1])
	if(initPos[2] < 0) then
		phi = - phi
	end
		
	return {-0.5*math.sin(0.5*t + phi), 0.5*math.cos(0.5*t + phi)}
end

function StaticBoundary(pos, initPos, t) 
	return {0, 0}
end
DownFixed = true
UpFixed = true
LeftFixed = true
RightFixed = true

function initStates(pos) 
	return {10*math.exp(-((pos[1]-0.5)^2 + (pos[2]-0.5)^2)/0.05)} 
    --return {10*math.exp(-((pos[2]-0.5)^2)/0.05)}, false
end

function DownQ(pos, initPos, t) 
	return {0}
end

function RightT(pos, initPos, t) 
	return {100}
end

function UpQ(pos, initPos, t) 
	return {0}
end

function LeftQ(pos, initPos, t) 
	return {23700.0}
end

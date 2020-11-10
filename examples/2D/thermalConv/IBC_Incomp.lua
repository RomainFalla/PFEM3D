TopFixed = true
BottomFixed = true
LeftFixed = true
RightFixed = true

function initStates(pos)
    return {0, 0, 0, 300};
end

function TopQ(pos, initPos, t) 
	return {0}
end

function TopV(pos, initPos, t) 
	return {0, 0}
end

function BottomQ(pos, initPos, t) 
	return {0}
end

function BottomV(pos, initPos, t) 
	return {0, 0}
end

function LeftT(pos, initPos, t) 
	return {310}
end

function LeftV(pos, initPos, t) 
	return {0, 0}
end

function RightT(pos, initPos, t) 
	return {290}
end

function RightV(pos, initPos, t) 
	return {0, 0}
end

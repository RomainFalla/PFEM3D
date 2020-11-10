TopFixed = true
BottomFixed = true
LeftFixed = true
RightFixed = true

function initStates(pos)
    return {0, 0, 0, 650};
end

function initTopStates(pos)
    return {0, 0, 0, 300};
end

function initBottomStates(pos)
    return {0, 0, 0, 1000};
end

function TopT(pos, initPos, t) 
	return {300}
end

function TopV(pos, initPos, t) 
	return {0, 0}
end

function BottomT(pos, initPos, t) 
	return {1000}
end

function BottomV(pos, initPos, t) 
	return {0, 0}
end

function LeftQ(pos, initPos, t) 
	return {0}
end

function LeftV(pos, initPos, t) 
	return {0, 0}
end

function RightQ(pos, initPos, t) 
	return {0}
end

function RightV(pos, initPos, t) 
	return {0, 0}
end

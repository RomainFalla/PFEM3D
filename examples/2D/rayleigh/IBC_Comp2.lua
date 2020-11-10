TopFixed = true
BottomFixed = true
LeftFixed = true
RightFixed = true

function initStates(pos) 
    return {0, 0, 0, rho0, 0, 0};
end

function initTopStates(pos) 
    return {0, 0, 0, rho0, 0, 0};
end

function initBottomStates(pos) 
    return {0, 0, 0, rho0, 0, 0};
end

function Top(pos, initPos, t) 
	return {0, 0}
end

function Bottom(pos, initPos, t) 
	return {0, 0}
end

function Left(pos, initPos, t) 
	return {0, 0}
end

function Right(pos, initPos, t) 
	return {0, 0}
end

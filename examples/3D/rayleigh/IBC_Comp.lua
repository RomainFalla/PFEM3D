TopFixed = true
BottomFixed = true
LateralFixed = true
z0 = 1

function initStates(pos) 
    local rho, p
    rho = rho0*((K0prime - 1)/K0*rho0*g*(z0 - pos[3]) + 1)^(1/(K0prime - 1))
    p = K0/K0prime*((rho/rho0)^K0prime - 1)
    return {0, 0, 0, p, rho, 0, 0, 0, 650};
end

function initTopStates(pos) 
    local rho, p
    rho = rho0*((K0prime - 1)/K0*rho0*g*(z0 - pos[3]) + 1)^(1/(K0prime - 1))
    p = K0/K0prime*((rho/rho0)^K0prime - 1)
    return {0, 0, 0, p, rho, 0, 0, 0, 300};
end

function initBottomStates(pos) 
    local rho, p
    rho = rho0*((K0prime - 1)/K0*rho0*g*(z0 - pos[3]) + 1)^(1/(K0prime - 1))
    p = K0/K0prime*((rho/rho0)^K0prime - 1)
    return {0, 0, 0, p, rho, 0, 0, 0, 1000};
end

function TopT(pos, initPos, t) 
	return {300}
end

function TopV(pos, initPos, t) 
	return {0, 0, 0}
end

function BottomT(pos, initPos, t) 
	return {1000}
end

function BottomV(pos, initPos, t) 
	return {0, 0, 0}
end

function LateralQ(pos, initPos, t) 
	return {0}
end

function LateralV(pos, initPos, t) 
	return {0, 0, 0}
end

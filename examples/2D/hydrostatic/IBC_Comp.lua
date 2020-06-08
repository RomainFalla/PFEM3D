z0 = 5
function initFluid(pos) 
	local rho, p
	rho = rho0*((K0prime - 1)/K0*rho0*g*(z0 - pos[2]) + 1)^(1/(K0prime - 1))
	p = K0/K0prime*((rho/rho0)^K0prime - 1)
	return {0, 0, p, rho, 0, 0}, false 
end

function initBoundary(pos)
    local rho, p
	if(pos[2] <= 5) then
		rho = rho0*((K0prime - 1)/K0*rho0*g*(z0 - pos[2]) + 1)^(1/(K0prime - 1))
	    p = K0/K0prime*((rho/rho0)^K0prime - 1)
	else
		rho = rho0
		p = 0
	end
	return {0, 0, p, rho, 0, 0}, false
end

function Boundary(pos, initPos, t) 
	return {0, 0}
end
z0 = 1

function initFluid(pos) 
	rho = rho0*(((K0prime-1)/K0)*(rho0*g*(z0-pos[2])) + 1)^(1/(K0prime -1))
	p = (K0/K0prime)*((rho/rho0)^K0prime - 1)

	return {0, 0, p, rho, 0, 0}, false 
end

function initBoundary(pos)
	rho = rho0*(((K0prime-1)/K0)*(rho0*g*(z0-pos[2])) + 1)^(1/(K0prime -1))
	p = (K0/K0prime)*((rho/rho0)^K0prime - 1)
	
	if(pos[2] > z0) then
		p = 0
		rho = 1000
	end

	return {0, 0, p, rho, 0, 0}, false
end

function initDiskBoundary(pos) 
	rho = rho0*(((K0prime-1)/K0)*(rho0*g*(z0-pos[2])) + 1)^(1/(K0prime -1))
	p = (K0/K0prime)*((rho/rho0)^K0prime - 1)

	return {0.25, 0, p, rho, 0, 0}, false	
end

function Boundary(pos, initPos, t) 
	return {0, 0}
end

function DiskBoundary(pos, initPos, t)
	return {0, 0}
end
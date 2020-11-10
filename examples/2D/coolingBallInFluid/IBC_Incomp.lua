BoundaryFixed = true

function initStates(pos)
    local T = 300
	if(pos[2] > 1 and pos[1] > 0 and pos[1] < 4) then
        T = 1000
    end
        
    return {0, 0, 0, T};
end

function BoundaryQ(pos, initPos, t) 
	return {0}
end

function BoundaryV(pos, initPos, t) 
	return {0, 0}
end

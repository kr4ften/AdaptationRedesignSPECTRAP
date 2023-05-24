-- SIMION workbench program used for tuning vertical beamline

simion.workbench_program()

-- Don't u dare forget to adjust this or u will be sorry
local beamID = 'beam1' -- Ensures csv outputs and inputs align, check ur CSV files

-- Designate splat volume
local xmin = -68
local xmax = 68

local ymin = -1450
local ymax = -1330

local zmin = 48
local zmax = 200

local SplatCount = 0
local SplatGood = 0

-- START OF READ CSV
-- Reads CSV files containing initial conditions for SIMION
-- NOTE: CSV-files typically assume flight in z-direction rearange if you want
-- to start them in another direction

local offset = {5.5, 58.5, -277} -- Position to start ions around

local function parseBeam(filename)
	print("Parsing beam!")
    local function parseLine(line, sep)
        local res = {}
        local pos = 1
        sep = sep
        while true do -- loop through characters in current line
            local c = string.sub(line, pos, pos)
            if (c == "") then break  -- If at end of line then break loop
            else
                local startp, endp = string.find(line, sep, pos)
                if (startp) then
                    table.insert(res, string.sub(line,pos,startp-1))
                    pos = endp + 1
                else 
                    table.insert(res, string.sub(line, pos))
                    break
                end
            end
        end
        return res
    end

    local x = {}
    local y = {}
    local l = {}
    local vx = {}
    local vy = {}
    local vz = {}
    --local E = {}

    local i = 0

    for cline in io.lines(filename) do 
        if i > 0 then --skipping first line
            local parsed = parseLine(cline, ';')
            x[i] = parsed[1]
            y[i] = parsed[2]
            l[i] = parsed[3]
            vx[i] = parsed[4]
            vy[i] = parsed[5]
            vz[i] = parsed[6]
        end
        i = i + 1
    end

    return x, y, l, vx, vy, vz
end

local x0, y0, l0, vx0, vy0, vz0 = parseBeam('input/' .. beamID .. '.csv')
local initiate_ion = 0 -- must be reset if rerunning, ions will run out otherwise
print("Beam parsed!") -- If you did not see this printed something is wrong

function segment.initialize()
    initiate_ion = initiate_ion + 1
    local c = initiate_ion
	--print("Running Initialize!", c)

    -- set read positions in mm
    ion_px_mm = x0[c] * 1e3 + offset[1]
    ion_py_mm = y0[c] * 1e3 + offset[2]
    ion_pz_mm = l0[c] * 1e3 + offset[3]
    -- set read velocities in mm/usec (1 mm/usec = 1e3 m/s)
    ion_vx_mm = vx0[c] * 1e-3
    ion_vy_mm = vy0[c] * 1e-3
    ion_vz_mm = vz0[c] * 1e-3
	
end

-- END OF READ CSV

-- Save ions to csv file
local function ion2csv(x, y, z, t, vx, vy, vz, m, run)
    local row = {}

    local file = io.open("output/result-" .. run .. ".csv", "w")
	file:write("Vertical Beamline - " .. beamID)
    file:write("Format: x [m]; y [m]; z [m]; t [s]; vx[m/s]; vy [m/s]; vz[m/s]; m [amu] \n")

    for i, x_el in pairs(x)
    do
        row = {tostring(x[i] * 1e-3), tostring(y[i] * 1e-3), tostring(z[i] * 1e-3), tostring(t[i] * 1e-6), tostring(vx[i] * 1e3), tostring(vy[i] * 1e3), tostring(vz[i] * 1e3), tostring(m)}
        file:write(table.concat(row, ';'), "\n")
    end
	io.close(file)

end

-- Check that ions are inside designated splat volume
local function check_bounds(x, y, z)
	if x > xmin and x < xmax and y > ymin and y < ymax and z > zmin and z < zmax then
		return true
	else
		return false
	end
end

-- Global-ish lists for storing splat locations
local SplatY = {}
local SplatX = {}
local SplatZ = {}
local SplatTime = {}

local SplatVX = {}
local SplatVY = {}
local SplatVZ = {}

local SplatM

local run = 0
local N = 1

local function SaveIons(run)
	-- Splatting in -y-direction
	ion2csv(SplatX, SplatY, SplatZ, SplatTime, SplatVX, SplatVY, SplatVZ, SplatM, run)
	-- Splatting in +z-direction
    --ion2csv(SplatX, SplatY, SplatTime, SplatVX, SplatVY, SplatVZ, SplatM, run)
end

-- Could be used for multiple reruns, useful if implementing optimiser probably
local function restart()
	initiate_ion = 0
	SplatY = {}
	SplatX = {}
	SplatZ = {}
	
	SplatTime = {}

	SplatVX = {}
	SplatVY = {}
	SplatVZ = {}
	
	SplatBad = 0
	SplatCount = 0
end

function segment.terminate()
	local function avg(tbl)
		local sum = 0
		if #tbl < 1 then
			return 0
		else
			for i, el in pairs(tbl) do
				sum = sum + el
			end
		end
		return sum / #tbl
	end
	
	local function std(tbl)
		local a = avg(tbl)
		local sum = 0
		if #tbl < 2 then
			return 0
		else
			for i, el in pairs(tbl) do
				sum = sum + (el - a)^2
			end
		end
		return sqrt(sum / (#tbl - 1))
	end
	
	local function alpha(a, b)
		local sum = 0
		local count = 0
		local offset = avg(a)
		if #a < 1 then
			return 0
		else
			for i = 1,#a do
				if a[i] and b[i] then
					sum = sum + (a[i] - offset) * b[i]/(sqrt(SplatVX[i]^2 + SplatVY[i]^2 + SplatVZ[i]^2))
					count = count + 1
				end
			end
		end
		return sum / count
	end
	
	SplatCount = SplatCount + 1
	
	-- Save Data from current run
	if check_bounds(ion_px_mm, ion_py_mm, ion_pz_mm) then
		-- Only save good ions
		SplatGood = SplatGood + 1
		
		SplatX[ion_number] = ion_px_mm
		SplatY[ion_number] = ion_py_mm
		SplatZ[ion_number] = ion_pz_mm

		SplatTime[ion_number] = ion_time_of_flight

		SplatVX[ion_number] = ion_vx_mm
		SplatVY[ion_number] = ion_vy_mm
		SplatVZ[ion_number] = ion_vz_mm
		
		SplatM = ion_mass
	else
		print("SPLAT OUTSIDE")
	end
	
	if SplatCount == initiate_ion then
		print("Splats recorded" .. tostring(SplatCount))
		print("Inside region:" .. tostring(SplatGood) .. "/" .. tostring(SplatCount))
		
		-- Save result from fly'm
		SaveIons(tostring(run))

		print("- Result from run #".. tostring(run))
		-- The most useful things to note:
		print("avg_x: " .. tostring(avg(SplatX)) .. " # std_x: " .. tostring(std(SplatX)))
		print("avg_z: " .. tostring(avg(SplatZ)) .. " # std_z: " .. tostring(std(SplatZ)))
		print("avg_vx: " .. tostring(avg(SplatVX)) .. " # std_vx: " .. tostring(std(SplatVX)))
		print("avg_vz: " .. tostring(avg(SplatVZ)) .. " # std_vz: " .. tostring(std(SplatVZ)))
		
		-- Was hoping to be able to tune for these to but was too difficult to do by hand
		local xdiag = alpha(SplatX, SplatVX)
		local zdiag = alpha(SplatZ, SplatVZ)
		print("r12 = " .. tostring(xdiag))
		print("r34 = " .. tostring(zdiag))
		print("delta = " .. tostring(xdiag - zdiag))
		
		
		run = run + 1
		-- Start next run
		if run < N then
			sim_rerun_flym = 1
			restart()
		else
			sim_rerun_flym = 0
			run = 0
			restart()
		end
	end	
end


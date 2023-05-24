--FILE RECORDING ENERGY;TIME;EMITTANCE;FINAL BUNCHLENGTH;IONSSPLAT WHEN BOTH VOLTAGES ARE CHOSEN ON PDT1 and PDT2, USED FOR DIFFERENT TIME CONSTANTS
simion.workbench_program()

--The local variables
local switched = 0
local t_switched = 1
local flag = 1
local n_splats = 0
local particle_count = 0
local avg_switchedtime
local count = 1

--The local arrays
local x = {}
local xprime = {}
local vx = {}
local y = {}
local yprime = {}
local vy = {}
local z = {}
local zprime = {}
local vz = {}


local tofs = {}
local velocity_xyz = {}
local splat_time = {}
local energy_xyz = {}
local pos_current = {}
local switch_time = {}

--The adjustable variables
adjustable switch_z = 200 -- adjust z-coordinate to the middle of the crown-part in workbench
adjustable first_tube = 2900 -- SET FINAL VOLTAGE
adjustable second_tube = 3700 -- SET FINAL VOLTAGE
adjustable n_particles = 350 -- set amount of ions
adjustable n_charge = 13 -- set the ions' charge
adjustable tau = 0.1 -- DEFAULT 0.1 usec, change if required!

-- define splat region
adjustable x_low = -2
adjustable x_high = 32
adjustable y_low = -2
adjustable y_high = 32
adjustable z_low = 403
adjustable z_high = 407

-- 3 = Beam trajectory in z-direction in workbench, change if the trajectory is defined differently
local mode = 3

-- INITIALIZE

local offset = {15, 15, -470} -- SET ION's WORKBENCH ORIGIN

-- Read input beam
local function parseBeam(filename)
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

local x0, y0, l0, vx0, vy0, vz0 = parseBeam('GeneratedBeam.csv')  -- file from Python in csv-file!
local initiate_ion = 100 -- use ions after the first 100 in the csv-file!

--Create ions in SIMION
function segment.initialize()
    initiate_ion = initiate_ion + 1
    local c = initiate_ion
    -- set read positions in mm
    ion_px_mm = x0[c] * 1e3 + offset[1]
    ion_py_mm = y0[c] * 1e3 + offset[2]
    ion_pz_mm = l0[c] * 1e3 + offset[3]
    -- set read velocities in mm/usec (1 mm/usec = 1e3 m/s)
    ion_vx_mm = vx0[c] * 1e-3
    ion_vy_mm = vy0[c] * 1e-3
    ion_vz_mm = vz0[c] * 1e-3
end


--FUNCTIONS - NOT BUILT IN SIMION

-- Compute emittance (approximation), the emittance in the thesis is done in a seperate Python file
function compute_emittance(a, aprime, va, v0) -- will not be used, use Python 
    function avg(arr)
        local res = 0
        for i = 1,#arr
            do
                res = res + arr[i]
            end
        if #arr > 0 
        then
            res = res / #arr
        end
        return res
    end

    local a_avg = avg(a)
    local aprime_avg = avg(aprime)
    local t = {}
    for i = 1, #a do
        t[i] = (a[i] - a_avg)^2
    end
    local da2_avg = avg(t)
    local t = {}
    for i = 1, #a do
        t[i] = (aprime[i] - aprime_avg)^2
    end
    local daprime2_avg = avg(t)
    local t = {}
    for i = 1, #a do
        t[i] = (a[i] - a_avg) * (aprime[i] - aprime_avg) 
    end
    local da_daprime_avg = avg(t)

    local m = da2_avg * daprime2_avg - da_daprime_avg^2
    if m < 0 then m = 0 end
    local emit = sqrt(m) * 1000 -- (mm * mrad)

    return emit
end

-- Used for emittance
function get_prime(va, vb)
    local aprime = {}
    for i = 1,#va
    do 
        aprime[i] = va[i]/vb[i]
    end
    return aprime
end

-- Check region for splats
function check_bounds(px, py, pz)
    if (px > x_low) and (px < x_high) and (py > y_low) and (py < y_high) and (pz > z_low) and (pz < z_high)
    then
        return true
    else
        return false
    end
end

-- Average for arrays
local function average(array)
    local sum = 0
    local n = 0
    if #array > 0 then 
        for i = 1, #array
        do  
            if array[i] ~= nil
            then
                sum = sum + array[i]
                n = n + 1
            end
        end

        return sum / n
    else
        return 0
    end
end 


-- FUNCTIONS - BUILT IN SIMION

-- Simulate the pulsing
function segment.fast_adjust()
    if average(pos_current) < switch_z then
        voltage1 = first_tube
        voltage2 = second_tube
    else
        switched = 1
        if t_switched == 1 then 
            t_switched = 0
            switch_time = tofs
            avg_switchedtime = average(switch_time)
        end
        voltage1 = first_tube * exp(-(average(tofs)- avg_switchedtime)/tau)
        voltage2 = second_tube * exp(-(average(tofs)- avg_switchedtime)/tau)
    end
    adj_elect01 = voltage2 -- the electrodes' names in fast adjust in Simion
    adj_elect02 = voltage1
end

-- Update colours and fill lists when ions splat
function segment.other_actions()
    pos_current[ion_number] = ion_pz_mm --z-positions of the ions
    tofs[ion_number] = ion_time_of_flight -- record flight time of the ions


    if switched == 1 then -- update PE and change colour when pulsed tubes
        ion_color = 2
        flag = 1
    elseif flag == 1 then
        flag = 0
        sim_update_pe_surface = 1
    end

    if ion_splat ~= 0 then -- if splat
        if check_bounds(ion_px_mm, ion_py_mm, ion_pz_mm) -- inside of defined volume
        then 
            particle_count = #x + 1

            splat_time[particle_count] = ion_time_of_flight -- usec 
            
            x[particle_count] = ion_px_mm -- mm
            y[particle_count] = ion_py_mm
            z[particle_count] = ion_pz_mm

            vx[particle_count] = ion_vx_mm -- mm/usec
            vy[particle_count] = ion_vy_mm
            vz[particle_count] = ion_vz_mm 

            velocity_xyz[particle_count] = sqrt(abs(vx[particle_count])^2 + abs(vy[particle_count])^2 + abs(vz[particle_count])^2) -- the total velocity

            energy_xyz[particle_count] = speed_to_ke(velocity_xyz[particle_count], ion_mass) -- eV
        else
            return -- splats outside of splat region
        end  
    else
        return -- no splat occured
    end
end


--FUNCTION - CSV-FILES VISUALIZE DATA

--save the energy
function write_csvenergy(energysp)
    print("Current run: ", count)
    local data = {}
    for i = 1, #energysp
    do
        data[i] = {energysp[i]}
    end
    local file = io.open("energylastruntau0.1.csv", "w")
    file:write("Format: energy \n")
    for i = 1, #energysp
    do
        file:write(table.concat(data[i], ';'), "\n")
    end
    io.close(file)
end

--save the splat time
function write_csvtime(splattime)
    local data = {}
    for i = 1, #splattime
    do
        data[i] = {splattime[i]}
    end
    local file = io.open("timelastruntau0.1.csv", "w")
    file:write("Format: time \n")
    for i = 1, #splattime
    do
        file:write(table.concat(data[i], ';'), "\n")
    end
    io.close(file)
end

-- output of all parameters required to calculate emittance
function write_csvemittance(a, va, b, vb, vc, splattime)
    local data = {}
    for i = 1, #a
    do
        data[i] = {a[i], va[i], b[i], vb[i], vc[i], splattime[i]}
    end
    local file = io.open("emittancelastruntau0.1.csv", "w")
    file:write("Format: x, vx, z, vz, vy, splat_time \n")
    for i = 1, #a
    do
        file:write(table.concat(data[i], ';'), "\n")
    end
    io.close(file)
end

-- save the final bunchlength and the amount of splats outside region
function write_onetime(bunchlts, ionsplat) 
    print('Done!')
    local data = {}
    data = {bunchlts, ionsplat}
    local file = io.open("onetimevaluelastruntau0.1.csv", "w")
    file:write("Format: ... \n")
    file:write(table.concat(data, ';'), "\n")
    io.close(file)
end

-- TERMINATE SEGMENT

function segment.terminate()
    if ion_number > 1 -- performed once
    then
        return
    end

    -- Quick check
    print("Amount of splats outside of the splat region: ", n_particles - particle_count)

    if mode == 1 then -- ion trajectory in x
        yprime = get_prime(vy, vx)
        zprime = get_prime(vz, vx)

        local emity = compute_emittance(y, yprime, vy, vx)
        local emitz = compute_emittance(z, zprime, vz, vx)
		
        write_csvenergy(energy_xyz)
        write_csvtime(splat_time)
        write_csvemittance(y, vy, z, vz, vx, splat_time)

        print("Emittance: (eps_y, eps_z) = (", emity, ", ", emitz, ")")

    elseif mode == 2 then -- ion trajectory in y
        xprime = get_prime(vx, vy)
        zprime = get_prime(vz, vy)

        local emitx = compute_emittance(x, xprime, vx, vy)
        local emitz = compute_emittance(z, zprime, vz, vy)

        write_csvenergy(energy_xyz)
        write_csvtime(splat_time)
        write_csvemittance(x, vx, z, vz, vy, splat_time)

        print("Emittance: (eps_x, eps_z) = (", emitx, ", ", emitz, ")")

    elseif mode == 3 then -- ion trajectory in z
        xprime = get_prime(vx, vz)
        yprime = get_prime(vy, vz)

        local emitx = compute_emittance(x, xprime, vx, vz)
        local emity = compute_emittance(y, yprime, vy, vz)

        write_csvenergy(energy_xyz)
        write_csvtime(splat_time)
        write_csvemittance(x, vx, y, vy, vz, splat_time)

        print("Emittance: (eps_x, eps_y) = (", emitx, ", ", emity, ")") -- quick check of emittance
    end

    table.sort(splat_time)

    local bunchlet =  average(velocity_xyz)*(splat_time[#splat_time] - splat_time[1]) --bunch length
    local ions_splat = (n_particles - particle_count)/n_particles -- amount of splats outside region

    write_onetime(bunchlet, ions_splat)

end



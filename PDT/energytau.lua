--Energy used to take out the mean and standard deviation, USED TO FOR DIFFERENT TIME CONSTANTS
simion.workbench_program()

--The local variables
local switched = 0
local t_switched = 1
local flag = 1
local particle_count = 0
local avg_switchedtime
local count = 1 -- NEW! BC OTHERWISE V1 IS BACKWARDS
local name = "energyresulttau0.1-x.csv" -- energy file for every ion hitting electrode (2), every run
local save_name = name


--The local arrays
local vx = {}
local vy = {}
local vz = {}


local tofs = {}
local velocity_xyz = {}
local splat_time = {}
local energy_xyz = {}
local pos_current = {}
local switch_time = {}

--The adjustable variables
adjustable switch_z = 200 -- adjust z-coordinate to the middle of the crown-part in workbench
adjustable first_tube = 2900 -- DECIDED FROM SIMULATION WITH BUNCHLENGTH!
adjustable second_tube = 2900 -- CHANGES!
adjustable n_particles = 100 -- set amount of ions
adjustable n_charge = 13 -- set the ions' charge
adjustable tau = 0.1 -- DEFAULT 0.1 usec, CHANGE
adjustable n_reruns = 32 -- amount of reruns
local dummy = n_reruns

-- define splat region
adjustable x_low = -2
adjustable x_high = 32
adjustable y_low = -2
adjustable y_high = 32
adjustable z_low = 403
adjustable z_high = 407

--The lists saving data between different runs
local voltage_2 = {}
local precentage_of_ions = {}


--INITIALIZE

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

local x0, y0, l0, vx0, vy0, vz0 = parseBeam('GeneratedBeam.csv')  -- file from Python in csv-file! (input beam)
local initiate_ion = 0 -- start in csv-file

-- Create ions in SIMION
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
        voltage1 = first_tube * exp(-(average(tofs)- avg_switchedtime)/tau) --- NEED TO RECORD THE TIME WHEN SWITCHED!!!!! TO GET SWITCH TIME
        voltage2 = second_tube * exp(-(average(tofs)- avg_switchedtime)/tau)
    end
    adj_elect01 = voltage2 -- the electrodes' names in fast adjust in Simion
    adj_elect02 = voltage1
end

-- Update colours and fill lists when ions splat
function segment.other_actions()
    pos_current[ion_number] = ion_pz_mm --z-positions of the ions
    tofs[ion_number] = ion_time_of_flight -- record flight time of the ions

    if dummy == n_runs then
        dummy = dummy - 1 -- will be false
        ion_color = 0 -- resets the color once every new run, black color
    end

    -- TAKE OUT A LIST SWITCH_TIME WITH TIME OF IONS WHEN SWITCHED
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
            particle_count = #vx + 1

            splat_time[particle_count] = ion_time_of_flight -- usec 

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
-- Creates an energy splat list for all ions, a new file for every rerun
function write_csvenergy(energysp)
    local numb = tostring(count)
    name = name:gsub("x", numb)
    print("Current run: ", numb)
    local data = {}
    for i = 1, #energysp
    do
        data[i] = {energysp[i]}
    end
    local file = io.open(name, "w")
    file:write("Format: energy \n")
    for i = 1, #energysp
    do
        file:write(table.concat(data[i], ';'), "\n")
    end
    io.close(file)
    name = save_name
end

-- Saves the voltage on tube (2) for every rerun (and amount of splats) gives out one list after all reruns has finished
function write_csv(voltage, prec_splat)
    local data = {}
    for i = 1, #voltage
    do
        data[i] = {voltage[i], prec_splat[i]}
    end
    local file = io.open("voltage2tau0.1.csv", "w")
    file:write("Format: voltage2ndtube, percentageSplats \n")
    for i = 1, #voltage
    do
        file:write(table.concat(data[i], ';'), "\n")
    end
    io.close(file)
end


-- TERMINATE SEGMENT

function segment.terminate()
    if ion_number > 1 -- performed once
    then
        return
    end

    Quick checks:
    print("Amount of splats outside of the splat region: ", n_particles - particle_count)
    print("Average energy when splatting: ", average(energy_xyz))

    table.sort(velocity_xyz)
    table.sort(splat_time)

    write_csvenergy(energy_xyz)

    voltage_2[count] = second_tube  -- saves voltage on the second tube
    precentage_of_ions[count] = (n_particles - particle_count)/n_particles  -- precentage of splats, 0% if no splats

    -- Next run
    n_reruns = n_reruns - 1
    sim_rerun_flym = (n_reruns > 0) and 1 or 0 -- if n_reruns > 0 gives 1 the program rerun, if n_runs is 0 the program don't rerun
    second_tube = second_tube + 50 -- change with step of two

    --When all runs are over
    if n_reruns == 0 then
        write_csv(voltage_2, precentage_of_ions)
    end

    --UPDATE VARIABLES AND ARRAYS AFTER EACH RUN
    dummy = n_reruns
    flag = 1
    t_switched = 1
    switched = 0
    particle_count = 0
    count = count + 1 -- UPDATE COUNT
    initiate_ion = 0 -- SAME IONS USED IN EVERY RERUN

    vx = {}
    vy = {}
    vz = {}

    tofs = {}
    velocity_xyz = {}
    splat_time = {}
    energy_xyz = {}
    pos_current = {}
    switch_time = {}

end

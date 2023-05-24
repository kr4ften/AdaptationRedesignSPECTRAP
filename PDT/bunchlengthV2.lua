--- Bunchlength as a function of voltage on PDT2, used for voltage 3100V on PDT1, THIS FILE IS ALSO USED FOR THE CHANGE IN POS OF ELECTRODE 2
simion.workbench_program()

--The local variables
local switched = 0
local t_switched = 1
local flag = 1
local particle_count = 0
local avg_switchedtime
local count = 1


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
adjustable first_tube = 3100 -- TEST FOR 3100V on PDT1
adjustable second_tube = 3250 -- PDT2 will change!
adjustable n_particles = 100  -- set amount of the ions
adjustable n_charge = 13 -- set charge on ion
adjustable tau = 0.01 -- time constant, default 0.01 usec
adjustable n_reruns = 15 -- set amount of reruns
local dummy = n_reruns 

--Define splat region, redefine when changed positions on electrode (2)!
adjustable x_low = -2
adjustable x_high = 32
adjustable y_low = -2
adjustable y_high = 32
adjustable z_low = 687
adjustable z_high = 700

--The lists that saves the wanted data between different runs
local bunchlength = {}
local voltages = {}

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

    if dummy == n_runs then
        dummy = dummy - 1 -- will be false
        ion_color = 0 -- resets the color once every new run, black color
    end

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
-- bunchlength and voltage on both tubes, gives one list after all reruns finishes
function write_csv(voltage, bunchlt)
    local data = {}
    for i = 1, #voltage
    do
        data[i] = {voltage[i], bunchlt[i]}
    end
    local file = io.open("resultbunchlt3100pos.csv", "w")
    file:write("Format: voltage2ndTube, bunchlentgh \n")
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

    -- quick checks
    print("Amount of splats outside of the splat region: ", n_particles - particle_count)
    print("Average energy when splatting: ", average(energy_xyz))  -- just a quick check that the energy given is expected
    print("Current run: ", count)

    table.sort(splat_time)

    -- Saved data between runs
    voltages[count] = second_tube  -- voltage1 = voltage2 for these simulation!
    bunchlength[count] = average(velocity_xyz)*(splat_time[#splat_time]-splat_time[1]) --- gives error if no ions splat in splat-region, v_avg*delta_t in splat

    -- Next run
    n_reruns = n_reruns - 1
    sim_rerun_flym = (n_reruns > 0) and 1 or 0 -- if n_reruns > 0 gives 1 the program rerun, if n_runs is 0 the program don't rerun
    second_tube = second_tube + 50 -- PDT2 change with 50 V

    --When all runs are over
    if n_reruns == 0 then 
        print('Done!')
        write_csv(voltages, bunchlength)
    end

    --UPDATE VARIABLES AND ARRAYS AFTER EACH RUN
    dummy = n_reruns
    flag = 1
    t_switched = 1
    switched = 0
    particle_count = 0
    count = count + 1
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

include("compute_orbits.jl")

matlab_exe = "/Applications/MATLAB_R2023a.app/bin/matlab"

k = sqrt(20)

for period = 5:26
    compute_solutions(period,k)
    run(`$matlab_exe -batch "addpath('matlab'); process_orbit_data"`)
end
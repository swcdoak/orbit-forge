


function process_orbit_data()
    xs = readmatrix('data/data.csv');
    period = size(xs, 2);

    % sort into canonical representation (min entry first)
    for i = 1:size(xs, 1)
        x = xs(i, :);
        [~, min_el_idx] = min(x);
        xs(i, :) = circshift(x, -(min_el_idx-1));
    end

    % remove duplicate solutions
    unique_sol = uniquetol(xs, 1e-4, "ByRows", true);

    orbits = [];
    for i = 1:size(unique_sol, 1)
        x = unique_sol(i, :);
        orbit = struct;
        orbit.points = [x; x(2:end) x(1); x(3:end) x(1:2)];

        if size(uniquetol(orbit.points', 1e-4, "ByRows", true), 1) == size(orbit.points, 2)
            orbits = [orbits orbit];
        end
    end

    size(orbits)

    load("data/orbit_data.mat", "Orbits")
    Orbits{period} = orbits;
    save("data/orbit_data.mat", "Orbits", "-append")
    
end






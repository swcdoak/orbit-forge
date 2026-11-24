using NLsolve, Primes, ProgressBars, DelimitedFiles


# ================================================
#             NECKLACE (SEED SOLUTIONS) FUNCTIONS
# ================================================

function allfactors(N)
    pfacs = factor(Vector, N)
    factors = [1, pfacs[1]]
    for fi in pfacs[2:end]
        factors = [1, fi]*factors'
        factors = unique(vec(factors))
    end
    return factors
end

function nNecklaces(n)
    ds = allfactors(n)
    N = 0
    for i in 1:length(ds)
        N += totient(ds[i]) * 2^ds[end-i+1]
    end
    return N÷n
end

function Theta_fast(b)
    # Compute j 
    j = findlast(x->x==1, b)[2]
    # Produce b 
    b[j] = 0
    for t in 1:length(b)-j 
        b[t+j] = b[t]
    end
    return j
end

# ================================================
#             SYSTEM FUNCTIONS
# ================================================

function F!(y, x, k)
    PERIOD = length(x)
    for i in 1:PERIOD-3
        y[i] = x[i+3] + k*(x[i+1]^2 + x[i+2]^2 - 2) - x[i]
    end
    y[PERIOD-2] = x[1] + k*(x[PERIOD-1]^2 + x[PERIOD]^2 - 2) - x[PERIOD-2]
    y[PERIOD-1] = x[2] + k*(x[PERIOD]^2 + x[1]^2 - 2) - x[PERIOD-1]
    y[PERIOD] = x[3] + k*(x[1]^2 + x[2]^2 - 2) - x[PERIOD]
    nothing
end

function Jac!(J, x, k)
    J = fill!(J, 0.0)
    PERIOD = length(x)
    for i in 1:PERIOD-3
        J[i, i] = -1; # Fixed point condition
        J[i, i+1] = 2*k*x[i+1];
        J[i, i+2] = 2*k*x[i+2];
        J[i, i+3] = 1;
    end
    
    # Third to last row
    J[PERIOD-2, 1] = 1;
    J[PERIOD-2, PERIOD-2] = -1;
    J[PERIOD-2, PERIOD-1] = 2*k*x[PERIOD-1];
    J[PERIOD-2, PERIOD] = 2*k*x[PERIOD];
    
    # Second to last row
    J[PERIOD-1, 1] = 2*k*x[1];
    J[PERIOD-1, 2] = 1;
    J[PERIOD-1, PERIOD-1] = -1;
    J[PERIOD-1, PERIOD] = 2*k*x[PERIOD];
    
    # Final to last row
    J[PERIOD, 1] = 2*k*x[1];
    J[PERIOD, 2] = 2*k*x[2];
    J[PERIOD, 3] = 1;
    J[PERIOD, PERIOD] = -1;
    nothing
end




# ================================================
#             SOLUTION CODE
# ================================================

function compute_solutions(period, k)
    N =  nNecklaces(period)
    xs = zeros(N, period)

    F_k! = (y, x) -> F!(y, x, k)
    Jac_k! = (J, x) -> Jac!(J, x, k)

    necklace = ones(1, period)
    x0 = 2*necklace .- 1
    solution = nlsolve(F_k!, Jac_k!, x0, method=:newton, ftol=1e-12, iterations=12)
    xs[1, :] = solution.zero
    
    next = 2
    for s=ProgressBar(1:N-1)
        j = period - 1
        while period % j != 0
            j = Theta_fast(necklace)
        end
        
        x0 = 2*necklace .- 1
        solution = nlsolve(F_k!, Jac_k!, x0, method=:newton, ftol=1e-14, iterations=30)

        # Skip if solution did not converge
        if !solution.f_converged
            continue
        end

        xs[next, :] = solution.zero
        next += 1
    end

    xs = xs[1:next-1, :]
    writedlm("data/data.csv",  xs, ',')

    println("Number of converged solutions: ", next-1)
    println("Proportion of seed solutions converged: ", (next-1)/N)
end

# k = sqrt(20)
# period = 7
# compute_solutions(period,k)

 
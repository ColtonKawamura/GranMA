export
    makeJobList2D
    
function makeJobList2D(filename::String, K_values::Vector{T1}, M_values::Vector{T2}, Bv_values::Vector{T3}, w_D_values::Vector{T4}, N_values::Vector{T5}, P_values::Vector{T6}, W_values::Vector{T7}, seeds::Vector{T8}) where {T1, T2, T3, T4, T5, T6, T7, T8}
    # makeJobList2D("testfile.txt", [100],[1],exp10.(-5:.2:2),exp10.(-5:.2:2),[10000],[.01,.001],[10,20,50],[1,2,3,4,5])
    # makeJobList2D("simulation_job_list.txt", [100],[1],exp10.(-5:.2:2),exp10.(-5:.2:2),[5000],[0.0001,0.00021544, 0.00046416, 0.001, 0.0021544, 0.0046416, 0.01, 0.021544, 0.046416, 0.1, ],[5],[1,2,3,4,5])    # Convert input vectors to the required types
    K_values = Int64.(K_values)
    M_values = Int64.(M_values)
    Bv_values = Float64.(Bv_values)
    w_D_values = Float64.(w_D_values)
    N_values = Float64.(N_values)
    P_values = Float64.(P_values)
    W_values = Int64.(W_values)
    seeds = Int64.(seeds)

    # Function to generate the MATLAB command EXAMPLE simulation_2d(K, M, Bv, w_D, N, P, W, seed)
    function generate_matlab_command(K, M, Bv, w_D, N, P, W, seed)
        return "matlab -nodisplay -nosplash -r \"addpath('./src/matlab_functions/'); simulation_2d($K, $M, $Bv, $w_D, $N, $P, $W, $seed); exit\""
    end

    # Generate all combinations of parameters using IterTools.product
    combinations = IterTools.product(K_values, M_values, Bv_values, w_D_values, N_values, P_values, W_values, seeds)

    # Open a file to write the MATLAB commands
    open(filename, "w") do file
        for combo in combinations
            command = generate_matlab_command(combo...)
            println(file, command)
        end
    end

    println("Commands written to $filename")
end
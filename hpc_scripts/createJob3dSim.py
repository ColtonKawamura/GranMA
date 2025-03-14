import numpy as np

def generate_matlab_command(K, M, Bv, w_D, N, P, W, seed, in_path, out_path):
    return (
        f"matlab -nodisplay -nosplash -r \"addpath('./src/matlab_functions/'); "
        f"sim3d({K}, {M}, {Bv}, {w_D}, {N}, {P}, {W}, {seed}, '{in_path}', '{out_path}'); exit\""
    )

def main():
    # Define different values for each variable
    K_values = [100]
    M_values = [1]
    Bv_values = [.01, .1, 1, 10, 100] 
    w_D_values = [.01, .1, 1, 10, 100] 
    N_values = [75600]
    P_values = [0.1, 0.01, 0.001]
    W_values = [15] # width of the simulation box
    seed_values = [1, 2, 3]
    in_path = 'in/3d/tiled_15by300/'  # Example input path
    out_path = 'out/simulation_3d/three_pressures_V1/'  # Example output path
    output_file = "./matlab_commands.txt"

    with open(output_file, "w") as file:
        # Generate and write MATLAB commands for all combinations
        for K in K_values:
            for M in M_values:
                for Bv in Bv_values:
                    for w_D in w_D_values:
                        for N in N_values:
                            for P in P_values:
                                for W in W_values:
                                    for seed in seed_values:
                                        command = generate_matlab_command(K, M, Bv, w_D, N, P, W, seed, in_path, out_path)
                                        file.write(command + "\n")  # Write to file

if __name__ == "__main__":
    main()

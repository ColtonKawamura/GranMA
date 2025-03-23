import numpy as np

def generate_matlab_command(K, M, Bv, w_D, N, P, W, seed, in_path, out_path):
    return (
        f"matlab -nodisplay -nosplash -r \"addpath('./src/matlab_functions/'); "
        f"sim2d({K}, {M}, {Bv}, {w_D}, {N}, {P}, {W}, {seed}, '{in_path}', '{out_path}'); exit\""
    )

def main():
    # Define different values for each variable
    K_values = [100]
    M_values = [1]
    # Bv_values = [.01, .05, .1, .5, 1, 5, 10,50, 100] 
    Bv_values = np.linspace(0.001, 100, 10)
    # w_D_values = [.01, .05, .1, .5, 1, 5, 10,50, 100]
    w_D_values = np.linspace(0.001, 100, 10)
    N_values = [80000]
    P_values = [0.1, 0.01, 0.001]
    W_values = [40] # width of the simulation box
    seed_values = [1]
    in_path = 'in/2d_poly_20by20/'  
    out_path = 'out/simulation_2d/poly_80kby40/'
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
                                        file.write(command + "\n") 

if __name__ == "__main__":
    main()

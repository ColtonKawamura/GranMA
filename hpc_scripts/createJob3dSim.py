import numpy as np
import subprocess

def generate_matlab_command(K, M, Bv, w_D, N, P, W, seed, in_path, out_path):
    return (
        f"matlab -nodisplay -nosplash -r \"addpath('./src/matlab_functions/'); "
        f"sim3d({K}, {M}, {Bv}, {w_D}, {N}, {P}, {W}, {seed}, '{in_path}', '{out_path}'); exit\""
    )

def main():
    # Define different values for each variable
    K_values = [100]
    M_values = [1, 2]
    Bv_values = np.logspace(-3, 0, num=5)  # Example logarithmic range
    w_D_values = [0.001, 0.005, 0.01]
    N_values = [40000, 80000, 160000]
    P_values = [0.05, 0.1, 0.2]
    W_values = [20, 40, 60]
    seed_values = [1, 2, 3]
    in_path = 'in/3d_tiled_2000by40/'  # Example input path
    out_path = 'out/simulation_3d/shear_80kby40/'  # Example output path

    # Generate and execute MATLAB commands for all combinations
    for K in K_values:
        for M in M_values:
            for Bv in Bv_values:
                for w_D in w_D_values:
                    for N in N_values:
                        for P in P_values:
                            for W in W_values:
                                for seed in seed_values:
                                    command = generate_matlab_command(K, M, Bv, w_D, N, P, W, seed, in_path, out_path)
                                    print(command)  # Print for verification
                                    subprocess.run(command, shell=True)

if __name__ == "__main__":
    main()
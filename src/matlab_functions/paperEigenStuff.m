% Damped
load("out/2d_damped_eigenStuff/2D_damped_eigenstuff_N400_40by13_K100_M1.mat", "data")

plotDampedModeDensityPDF(data, [.2, .01, .001], [.1])


% UnDamped
file_name_list1 = [
    "in/packings_processed_eig_PDF/2D_N1600_P0.1_Width40_Seed1.mat", 
    "in/packings_processed_eig_PDF/2D_N1600_P0.01_Width40_Seed1.mat", 
    "in/packings_processed_eig_PDF/2D_N1600_P0.001_Width40_Seed1.mat"
];
 plotModeDensityPDF(file_name_list, true)

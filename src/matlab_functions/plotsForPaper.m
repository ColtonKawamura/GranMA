clear all
close all
addpath("src/matlab_functions")
% 3D
data = load('out/processed/3d_80Kby15_V4_1.mat')
gamma_value = .2;
[omega_gamma_values, attenuation_values] = plotAttenuation(data, gamma_value);
plotWavespeed(data, gamma_value);


% 2D
data = load("out/processed/2d_K100_80kby40.mat")
gamma_value = .5;
[omega_gamma_values, attenuation_values] = plotAttenuation(data, gamma_value);
plotWavespeed(data, gamma_value);
filtered_data = filterData(data, 'gamma', .1 , 'pressure', 0.1, 'omega', .1, 'seed', 1);
[mean_field_amp, mean_field_phase, prime_field_amp, prime_field_phase] = plotMeanField(filtered_data)

%  2D Stiched Attenution
gamma_values = [ .5, .6, .8, 1]
plotAttenuationMult(data, gamma_values, 1.2)
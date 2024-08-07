using MAT
using Glob
using Debugger
using JLD2



function plot_ellipse_pdf(simulation_data, gamma_value)
    closest_gamma_index = argmin(abs.([idx.gamma for idx in simulation_data] .- gamma_value))
    matching_gamma_data = filter(entry -> entry.gamma == gamma_value, simulation_data)

    pressure_list = unique([entry.pressure for entry in matching_gamma_data]) # goes through each entry of simulation_data and get the P value at that entry
    
    normalized_variable = (log.(pressure_list) .- minimum(log.(pressure_list))) ./ (maximum(log.(pressure_list)) .- minimum(log.(pressure_list)))

    for pressure_value in pressure_list
        matching_pressure_data = filter(entry -> entry.pressure == pressure_value, matching_gamma_data) # for every entry in simluation_data, replace (->) that entry with result of the boolean expression
        marker_color = [normalized_variable[idx], 0, 1-normalized_variable[idx]]
        
        # Each omega gamma value spans all seeds
        omega_gamma_list = unique([entry.omega_gamma for entry in matching_pressure_data])
        for omega_gamma_value in omega_gamma_list
            matching_omega_gamma_data = filter(entry -> entry.omega_gamma == omega_gamma_value, matching_pressure_data)
            # mean(matrix, dims=2) means across rows. Dims = 1 is cross columns, dims =3 is into the screen
            #  If you have an array [[a1, a2, a3], [b1, b2, b3], [c1, c2, c3]], splatting converts this into calling hcat([a1, a2, a3], [b1, b2, b3], [c1, c2, c3]).
            mean_asp_rat_counts = mean(hcat([entry.asp_rat_counts for entry in matching_omega_gamma_data]...), dims=2)
            mean_asp_rat_bins = mean(hcat([entry.asp_rat_bins for entry in matching_omega_gamma_data]...), dims=2)
            mean_rot_ang_counts = mean(hcat([entry.rot_ang_counts for entry in matching_omega_gamma_data]...), dims=2)
            mean_rot_ang_bins = mean(hcat([entry.rot_ang_bins for entry in matching_omega_gamma_data]...), dims=2)

        end
        closest_asp_rat_counts = matching_pressure_data[closest_gamma_index].asp_rat_counts
    end


    closest_asp_rat_counts = simulation_data[closest_gamma_index].asp_rat_counts
end
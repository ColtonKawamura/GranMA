function plotAttenuationWithPressure(data, gammaValues)
    % arguments
    %     data (1,1) struct % Required argument: scalar struct containing simulation data fields
    %     gammaValues (1,:) double 
    %     options.shear (1,1) logical = false % Optional flag for shear vs compression
    %     options.z (1,1) logical = false % plots z direction instead
    % end

    % Define the plot limits to match the 1D theory plot curves
    theory_x = (3E-4:1E-5:3);
    theory_y = theory_x ./ sqrt(2) .* ((1 + theory_x.^2) .* (1 + sqrt(1 + theory_x.^2))).^(-0.5);

    figure_attenuation = figure;
    % loglog(theory_x, theory_y, 'k', 'DisplayName', '1-D Theory');
    hold on;
    xlabel('$\hat{\omega}\hat{\gamma}/\hat{P}^{3/2}$', 'FontSize', 20, 'Interpreter', 'latex');
    ylabel('$ \frac{\hat{\alpha}}{\hat{\omega}\hat{P}^{1/2}} $', 'FontSize', 20, 'Interpreter', 'latex');
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    set(get(gca, 'ylabel'), 'rotation', 0);
    grid on;
    box on;


    for gammaValue = gammaValues
        matchingGammaData = filterData(data, 'gamma', gammaValue);
        gammaValue = matchingGammaData.gamma(1);
    
        % Get list of unique pressures for the matching gamma data
        pressureList = unique(matchingGammaData.pressure);
        pressureList = sort(pressureList)
        % pressureList = [min(pressureList), max(pressureList)]; % for now, only plot min and max
    
        for pressureValue = pressureList'
            [pressureValueNormalized, markerColor] = normVarColor(pressureList, pressureValue, true);
            pressureValueData = filterData(matchingGammaData, 'pressure', pressureValue);
            omegaGammaValues = pressureValueData.omega_gamma;
            uniqueOmegaGammaValues = unique(omegaGammaValues);
            uniqueOmegaGammaValues = sort(uniqueOmegaGammaValues);
            
            meanAttenuationList = [];
            for omegaGammaValue = uniqueOmegaGammaValues'
                omegaGammaValueData = filterData(pressureValueData, 'omega_gamma', omegaGammaValue);
                seedList = omegaGammaValueData.seed;
                attenuationList = [];
                for seedValue = seedList'
                    seedData = filterData(omegaGammaValueData, 'seed', seedValue);
                    if isempty(seedData.unwrapped_phase_vector_y{1})
                        fprintf('Empty y-phase vector for: Pressure %f OmegaGamma %f seed %d\n', pressureValue, omegaGammaValue, seedValue);
                        continue
                    end
                    attenuationList = [attenuationList, seedData.alphaoveromega_x];
                end
                meanAttenuationList= [meanAttenuationList, mean(attenuationList)];

            end
            markerSize = exp(gammaValue/max(gammaValues))*3;
            pressureLabel = sprintf('$ %.4f, %.4f $', pressureValue, gammaValue); 
            plot(uniqueOmegaGammaValues/pressureValue^1.5, meanAttenuationList/pressureValue^.5, "-o", 'MarkerSize', markerSize, 'MarkerFaceColor', markerColor, 'Color', markerColor, 'DisplayName', pressureLabel);
        end
        leg = legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex', 'FontSize', 15);
        title(leg, "$  \hat{P}, \hat{\gamma} $")
    end
end
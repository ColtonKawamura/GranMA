function plotEnergyTheory(data, gammaValues)
    % arguments
    %     data (1,1) struct % Required argument: scalar struct containing simulation data fields
    %     gammaValues (1,:) double 
    %     options.shear (1,1) logical = false % Optional flag for shear vs compression
    %     options.z (1,1) logical = false % plots z direction instead
    % end
    ax_energy = figure;
    xlabel('$  \hat{\omega} \hat{\gamma} $', "FontSize", 20, "Interpreter", "latex");
    ylabel('$  \frac{2\hat{c}^2 \left( 1 - \cos\overline{\sigma}_{\Delta \phi_{ij}} \right)}{\hat{d^2\omega}^2}\left(\overline{\frac{A_{\perp}}{A_{\parallel}}}\right)^2 $', "FontSize", 20, "Interpreter", "latex");    set(gca, 'YScale', 'log');
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    grid on;
    box on; 
    hold on;


    for gammaValue = gammaValues
        matchingGammaData = filterData(data, 'gamma', gammaValue);
        gammaValue = matchingGammaData.gamma(1);
    
        % Get list of unique pressures for the matching gamma data
        pressureList = unique(matchingGammaData.pressure);
        pressureList = sort(pressureList);
        pressureList = [min(pressureList), max(pressureList)]; % for now, only plot min and max
    
        for pressureValue = pressureList
            [pressureValueNormalized, markerColor] = normVarColor(pressureList, pressureValue, true);
            pressureValueData = filterData(matchingGammaData, 'pressure', pressureValue);
            omegaGammaValues = pressureValueData.omega_gamma;
            uniqueOmegaGammaValues = unique(omegaGammaValues);
            uniqueOmegaGammaValues = sort(uniqueOmegaGammaValues)
            
            meanEnergyLossList = [];
            for omegaGammaValue = uniqueOmegaGammaValues'
                omegaGammaValueData = filterData(pressureValueData, 'omega_gamma', omegaGammaValue);
                seedList = omegaGammaValueData.seed;
                energyLossList = [];
                for seedValue = seedList'
                    seedData = filterData(omegaGammaValueData, 'seed', seedValue);
                    if isempty(seedData.unwrapped_phase_vector_y{1})
                        fprintf('Empty y-phase vector for: Pressure %f OmegaGamma %f seed %d\n', pressureValue, omegaGammaValue, seedValue);
                        continue
                    end
                    phaseScatter = 1-cos(plotPhase(seedData, 'plot', false));
                    ampRatio = plotAmp(seedData, 'plot', false);
                    energyLoss = 2*phaseScatter*ampRatio.^2/seedData.omega^2;
                    energyLossList = [energyLossList, energyLoss];
                end
                meanEnergyLossList= [meanEnergyLossList, mean(energyLossList)];

            end
            markerSize = exp(gammaValue/max(gammaValues))*3;
            pressureLabel = sprintf('$ %.4f, %.4f $', pressureValue, gammaValue); 
            plot(uniqueOmegaGammaValues, meanEnergyLossList, "-o", 'MarkerSize', markerSize, 'MarkerFaceColor', markerColor, 'Color', markerColor, 'DisplayName', pressureLabel);
        end
        leg = legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex', 'FontSize', 15);
        title(leg, "$  \hat{P}, \hat{\gamma} $")
    end
end
function plotPhaseScatter(data, gammaValues)
    % arguments
    %     data (1,1) struct % Required argument: scalar struct containing simulation data fields
    %     gammaValues (1,:) double 
    %     options.shear (1,1) logical = false % Optional flag for shear vs compression
    %     options.z (1,1) logical = false % plots z direction instead
    % end
    ax_energy = figure;
    xlabel('$  \hat{\omega} $', "FontSize", 20, "Interpreter", "latex");

    ylabel('$ 1-\cos \overline{\sigma}_{\Delta \phi_{\perp}} $', "FontSize", 20, "Interpreter", "latex");    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
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
            uniqueOmegaGammaValues = sort(uniqueOmegaGammaValues);
    
            meanPhaseScatterList = [];
            for omegaGammaValue = uniqueOmegaGammaValues'
                omegaGammaValueData = filterData(pressureValueData, 'omega_gamma', omegaGammaValue);
                seedList = omegaGammaValueData.seed;
                phaseScatterList = [];
                for seedValue = seedList'
                    seedData = filterData(omegaGammaValueData, 'seed', seedValue);
                    if isempty(seedData.unwrapped_phase_vector_y{1})
                        fprintf('Empty y-phase vector for: Pressure %f OmegaGamma %f seed %d\n', pressureValue, omegaGammaValue, seedValue);
                        continue
                    end
                    phaseScatter = plotPhase(seedData, 'plot', false);
                    phaseScatterList = [phaseScatterList, 1-cos(phaseScatter)];
                end
                meanPhaseScatter = mean(phaseScatterList);
                meanPhaseScatterList= [meanPhaseScatterList, meanPhaseScatter];

            end
            markerSize = exp(gammaValue/max(gammaValues))*3;
            pressureLabel = sprintf('$ %.4f, %.4f $', pressureValue, gammaValue); 
            plot(uniqueOmegaGammaValues/gammaValue, meanPhaseScatterList, "-o", 'MarkerSize', markerSize, 'MarkerFaceColor', markerColor, 'Color', markerColor, 'DisplayName', pressureLabel);
        end
        leg = legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex', 'FontSize', 15);
        title(leg, "$  \hat{P}, \hat{\gamma} $")
    end
end
function plotAmpRatio(data, gammaValues)

    ax_energy = figure;
    xlabel('\$  \\hat{\\omega} \$', "FontSize", 20, "Interpreter", "latex");

    %xlabel('\$\\hat{\\omega}\\hat{\\gamma}\$', "FontSize", 20, "Interpreter", "latex");
    ylabel('\$ \\hat{\\gamma} \\left( \\overline{\\frac{A_{\\perp}}{A_{\\parallel}}} \\right) ^2\$', "FontSize", 20, "Interpreter", "latex");
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    grid on;
    box on; 
    hold on;


    markerShapeVector = ["-*", "-o", "-v", "-+", "-.", "-x", "-d"];

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
    
            meanAmpRatioList = [];
            for omegaGammaValue = uniqueOmegaGammaValues'
                omegaGammaValueData = filterData(pressureValueData, 'omega_gamma', omegaGammaValue);
                seedList = omegaGammaValueData.seed;
                AmpRatioList = [];
                for seedValue = seedList'
                    seedData = filterData(omegaGammaValueData, 'seed', seedValue);
                    ampRatio = plotAmp(seedData, 'plot', false);
                    AmpRatioList = [AmpRatioList, ampRatio];
                end
                meanAmpRatio = mean(AmpRatioList);
                meanAmpRatioList= [meanAmpRatio];
                
            end
            markerShape = markerShapeVector(find(pressureList == pressureValue));
            pressureLabel = sprintf('$ %.4f, %.4f $', pressureValue, gammaValue); 
            plot(uniqueOmegaGammaValues, meanAmpRatioList, markerShape, 'Color', markerColor, 'DisplayName', pressureLabel);
        end
    
    end
    leg = legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex', 'FontSize', 15);
    title(leg, "$  \hat{P}, \hat{\gamma} $")
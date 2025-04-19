function plotAttenuationMult(simulationData, gammaValues, meanDiameter, isShear)

if nargin <4
    isShear = false;
end
if nargin < 3
    meanDiameter = 1.2;
end

theory_x = 3E-4:1E-5:3;
theory_y = theory_x ./ sqrt(2) .* ((1 + theory_x.^2) .* (1 + sqrt(1 + theory_x.^2))).^(-0.5);

figure_attenuation = figure;
loglog(theory_x, theory_y, 'k', 'DisplayName', '1-D Theory');
hold on;
xlabel('$\hat{\omega}\hat{\gamma}$', "FontSize", 20, "Interpreter", "latex");
ylabel('$ \frac{\hat{\alpha}}{\hat{\omega}} $', "FontSize", 20, "Interpreter", "latex");
set(gca, 'XScale', 'log');
set(get(gca, 'ylabel'), 'rotation', 0);
grid on;
box on;

markerShapeVector = ["-*", "-o", "-v", "-+", "-.", "-x", "-d"];

for gammaValue = gammaValues
    matchingGammaData = filterData(simulationData, 'gamma', gammaValue);
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

        meanAttenuationList = [];
        for omegaGammaValue = uniqueOmegaGammaValues'
            omegaGammaValueData = filterData(pressureValueData, 'omega_gamma', omegaGammaValue);
            attenuationValues = omegaGammaValueData.alphaoveromega_x;
            meanAttenuation = mean(attenuationValues);
            meanAttenuationList = [meanAttenuationList, meanAttenuation];
        end
        markerShape = markerShapeVector(find(pressureList == pressureValue));
        pressureLabel = sprintf('$ %.4f, %.4f $', pressureValue, gammaValue); 
        plot(uniqueOmegaGammaValues, meanAttenuationList, markerShape, 'Color', markerColor, 'DisplayName', pressureLabel);
    end

end
leg = legend('show', 'Location', 'northeastoutside', 'Interpreter', 'latex', 'FontSize', 15);
title(leg, "$  \hat{P}, \hat{\gamma} $")

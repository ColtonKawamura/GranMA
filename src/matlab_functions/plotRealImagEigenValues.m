function plotRealImagEigenValues(plotData)

    realEigenValues = real(plotData.eigenValues{1});
    imagEigenValues = imag(plotData.eigenValues{1});
    scatter(imagEigenValues, -realEigenValues, 20, realEigenValues, 'filled')
    % scatter(-realEigenValues, imagEigenValues, 20, realEigenValues, 'filled')
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    grid on
    box on
    ax = gca;
    ax.FontSize = 15;
    xlabel('Imaginary ($\gamma$)', 'Interpreter', 'latex', 'FontSize', 20)
    ylabel('Real ($\omega$)', 'Interpreter', 'latex', 'FontSize', 20)
    % title(sprintf('$L_x$ by $L_y$: %.2f by %.2f, $\gamma$: .2f', plotData.Lx(1), plotData.Ly(1)), plotData.damping, 'Interpreter', 'latex', 'FontSize', 16);
    title(sprintf('$L_x$ by $L_y$: %.2f by %.2f, $\\hat{P}$: %.3f, $\\hat{\\gamma}$: %.3f', plotData.Lx(1), plotData.Ly(1), plotData.pressure(1), plotData.damping(1)), 'Interpreter', 'latex', 'FontSize', 16);
end
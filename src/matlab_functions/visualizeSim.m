function visualizeSim(nt, x, x0, y, y0, idx, A)
    if mod(nt,1) == 0
        
        % Plot Δx
        subplot(1,2,1)
        plot(x0(idx), x(idx) - x0(idx), 'o')
        xlabel('$x_0$', 'Interpreter', 'latex', 'FontSize', 20)
        ylabel('$ \Delta x$', 'Interpreter', 'latex', 'FontSize', 20)
        ylim(1.2*[-A,A])
        grid on

        % Plot Δy
        subplot(1,2,2)
        plot(x0(idx), y(idx) - y0(idx), 'o')
        xlabel('$x_0$', 'Interpreter', 'latex', 'FontSize', 20)
        ylabel('$ \Delta y$', 'Interpreter', 'latex', 'FontSize', 20)
        ylim(1.2*[-A,A])
        grid on
        
        drawnow
    end
end

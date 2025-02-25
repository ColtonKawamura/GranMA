function visualizeSim(nt, x, x0, idx, A )
    if mod(nt,1)==0

            plot(x0(idx), x(idx) - x0(idx),'o')
            ylim(1.2*[-A,A])
            xlabel('$x_0$', 'Interpreter', 'latex', 'FontSize', 20)
            ylabel('$ \Delta x$', 'Interpreter', 'latex', 'FontSize', 20)
            grid on
        drawnow
    end
end
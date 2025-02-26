function visualizeSim(nt, x, x0, idx, A, flag_shear)
    if mod(nt,1)==0
            plot(x0(idx), x(idx) - x0(idx),'o')
            if flag_shear == true
                ylim(1.2*[-A*.1,A*.1])
            else
                ylim(1.2*[-A,A])
            end
            xlabel('$x_0$', 'Interpreter', 'latex', 'FontSize', 20)
            ylabel('$ \Delta x$', 'Interpreter', 'latex', 'FontSize', 20)
            grid on
        drawnow
    end
end
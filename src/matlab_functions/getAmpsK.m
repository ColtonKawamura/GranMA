function [breakOut, outAmp, outXinit, nt_out] = getAmpsK(nt, x, x0, idx, w_D, A, outAmp, outXinit, Lx, nt_out)
    breakOut = false;
    if mod(nt,10)==0


            plot(x0(idx), x(idx) - x0(idx), '.')
            ylim(1.2*[-A,A])
            outAmp = NaN;
            outXinit = NaN;
            nt_out = NaN;        % Stop once peak has attenuated past threshold, prevents from jumping, need to see if I need this anymore due to previous ones.
            % if second_peak_Amp < outAmp(1) * .7 % stop sim once 2nd peak amp is some percent of inital amp, .5 works well for low freq, but used .7 for others
            %     display("peak dropped")
            %     breakOut = true;
            %     return
            % end
        drawnow
    end
end
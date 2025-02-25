function [breakOut, outAmp, outXinit, nt_out] = getAmps(nt, x, x0, idx, w_D, A, outAmp, outXinit, Lx, nt_out)
    breakOut = false;
    if mod(nt,1)==0

        % smoothAmp = smooth(x(idx) - x0(idx), 20, 'sgolay');
        % smoothAmp = smoothdata(x(idx)- x0(idx), "SamplePoints", x0(idx));
        smoothAmp = smoothdata(x(idx)- x0(idx), "sgolay");
        % minPeakWidth = 1/w_D * .6; % Criteria is at a little more that a half wavelength, some buffer for weird stuff happening
        % [pks,locs]=findpeaks(smoothAmp, "MinPeakWidth", minPeakWidth);
        [pks,locs]=findpeaks(smoothAmp);
        valid_peaks_idx = find(pks > A * 0.08);
        if ~isempty(valid_peaks_idx) && length(valid_peaks_idx) > 2 % this should match the peak
            second_peak_idx = valid_peaks_idx(end-2); % this is the third peak
            second_peak_Amp = pks(second_peak_idx);
            second_peak_xInit = x0(idx(locs(second_peak_idx)));

            % Prevent from jumping to peak bheind
            % if length(outXinit) > 1 && second_peak_xInit < outXinit(end-1) 
                % return
            % end

            % Prevent from jumping to peak that is greater than mean of the last 3 amplitudes
            if length(outXinit) > 1 && second_peak_Amp > mean(outAmp) && second_peak_xInit < outXinit(end) % length(outXinit) should be same as average window
                return
            end

            % Prevent from jumping back to peak that is too much smaller (jumping ahead to lower peak or back to false peak)
            if length(outXinit) > 5 && second_peak_Amp < mean(outAmp(end-5:end)) && second_peak_xInit < outXinit(end)% length(outXinit) should be same as average window
                return
            end


            % Prevent from hitting back wall
            if second_peak_xInit > .8 * Lx 
                display("peak hit .8 wall")
                breakOut = true;
            end

            plot(x0(idx), x(idx) - x0(idx), '.', x0(idx), smoothAmp, 'r-', second_peak_xInit, second_peak_Amp, 'o', 'MarkerFaceColor', 'r')
            ylim(1.2*[-A,A])
            outAmp = [outAmp, second_peak_Amp];
            outXinit = [outXinit, second_peak_xInit];
            nt_out = [nt_out, nt];
            xlabel('$x_0$', 'Interpreter', 'latex', 'FontSize', 20)
            ylabel('$ \Delta x$', 'Interpreter', 'latex', 'FontSize', 20)
            grid on
            % Stop once peak has attenuated past threshold, prevents from jumping, need to see if I need this anymore due to previous ones.
            % if second_peak_Amp < outAmp(1) * .60 % stop sim once 2nd peak amp is some percent of inital amp, .5 works well for low freq, but used .7 for others
            if second_peak_Amp < A * .09 % stop sim once 2nd peak amp is some percent of inital amp, .5 works well for low freq, but used .7 for others
                display("peak dropped")
                breakOut = true;
                return
            end
        end
        drawnow
    end
end
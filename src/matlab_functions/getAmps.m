function [breakOut, outAmp, outXinit] = getAmps(nt, x, x0, idx, w_D, A, outAmp, outXinit)
    breakOut = false;
    if mod(nt,10)==0

        smoothAmp = smooth(x(idx) - x0(idx), 150, 'sgolay');
        minPeakWidth = 1/w_D * .6;
        [pks,locs]=findpeaks(smoothAmp, "MinPeakWidth", minPeakWidth);
        valid_peaks_idx = find(pks > A * 0.1);
        if ~isempty(valid_peaks_idx) && length(valid_peaks_idx) > 2
            second_peak_idx = valid_peaks_idx(end-2); % this is the third peak
            second_peak_Amp = pks(second_peak_idx);
            second_peak_xInit = x0(idx(locs(second_peak_idx)));
            % plot(x0(idx), x(idx) - x0(idx), '.', x0(idx), smoothAmp, 'r-', x0(idx(peak_index)), x(idx(peak_index)) - x0(idx(peak_index)), 'o', 'MarkerFaceColor', 'r')
            plot(x0(idx), x(idx) - x0(idx), '.', x0(idx), smoothAmp, 'r-', second_peak_xInit, second_peak_Amp, 'o', 'MarkerFaceColor', 'r')
            ylim(1.2*[-A,A])
            outAmp = [outAmp, second_peak_Amp];
            outXinit = [outXinit, second_peak_xInit];
            if second_peak_Amp < outAmp(1) * .5 % stop sim once 2nd peak amp is some percent of inital amp, .5 works well for low freq, but used .7 for others
                display("peak dropped")
                breakOut = true;
                return
            end
        end
        drawnow
    end
end
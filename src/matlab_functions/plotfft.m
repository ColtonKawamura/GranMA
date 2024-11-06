function [f,fft_out]=plotfft(data,fps,deltaf)

if ~exist('deltaf')
    deltaf=1;
end

T=1/fps;
L=length(data);
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(data-mean(data),NFFT)/L;
f = fps/2*linspace(0,1,NFFT/2+1);
fft_out=2*smooth(abs(Y(1:NFFT/2+1)),deltaf);
% Plot single-sided amplitude spectrum.
figure(1), plot(data)
figure(2), plot(f,fft_out,'r.')
end
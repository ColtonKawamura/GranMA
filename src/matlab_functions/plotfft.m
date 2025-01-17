function [f,fft_out]=plotfft(data,fps, marker_color)

deltaf = 1;

T=1/fps;
L=length(data);
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(data-mean(data),NFFT)/L;
f = fps/2*linspace(0,1,NFFT/2+1);
fft_out=2*smooth(abs(Y(1:NFFT/2+1)),deltaf);
% Plot single-sided amplitude spectrum.
% figure, plot(data)
% figure(100), plot(f, fft_out, strcat('-o', marker_color)), hold on
figure(100), stem(f, fft_out), hold on
xlabel("$ \omega ", "Interpreter", "latex");
xlabel("$ A $", "Interpreter", "latex");

grid on
end
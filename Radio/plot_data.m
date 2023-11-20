close all

%% Load data
load("RX06.mat");
% contains rx_capture (data), center_frequency [Hz], sample_rate [Hz], capture_time [s]


singal_lenght = capture_time * sample_rate;


x = sample_rate / singal_lenght * (-singal_lenght/2 : singal_lenght/2 - 1);
y = abs(fftshift(fft(rx_capture)));

figure();
plot(abs(fftshift(fft(rx_capture))))
figure();
plot(x, y)
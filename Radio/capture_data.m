clear all
close all
clc

%% Filename
prefix = "R_noise_";
num = "06";

%% Device setup
env_device = 'Pluto';
env_id = 'usb:0';

%% Capture setup
center_frequency = 1626.2708e6; % [Hz]
sample_rate = 10e6; % [Hz]
capture_time = 0.3; % [s]

%% Receiver declaration
rx = sdrrx(env_device,'RadioID', env_id,'CenterFrequency', center_frequency, 'BasebandSampleRate', sample_rate, 'OutputDataType','double', ...
    'ShowAdvancedProperties',true, 'FrequencyCorrection',0);

%% Capture data
rx_capture = capture(rx, capture_time, 'Seconds');
rx.info
release(rx)

save("Data\" + prefix + num + ".mat");
write_complex_binary(rx_capture, "Data\" + prefix + num + "_binary")

plot(abs(fftshift(fft(rx_capture))))
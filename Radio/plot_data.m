close all

%% Settings
file = "RX02.mat"; % contains rx_capture (data), center_frequency [Hz], sample_rate [Hz], capture_time [s]
max_doppler_shift = 36e3; % [Hz]

%% Prepare
load(file);
max_doppler_shift = max_doppler_shift * 1.1; % tolerance
singal_lenght = capture_time * sample_rate;

%% Create FFT and adjust axis
x = (sample_rate / singal_lenght * (-singal_lenght/2 : singal_lenght/2 - 1)) + center_frequency;
y = abs(fftshift(fft(rx_capture)));

%% Plot whole FFT
figure();
plot(x, y)
title("FFT")
xlabel('f [Hz]')

%% Plot FFT around center frequency
figure();
plot(x, y)
axis([center_frequency-max_doppler_shift  center_frequency+max_doppler_shift 0 inf])
title("FFT around center frequency (+- max Doppler shift)")
xlabel('f [Hz]') 

%% Plot sliced segments
segments = [6 15];
y_axis_max = 500;

num_segments = segments(1) * segments(2);
segment_size = floor(singal_lenght / num_segments);
segment_time = round(capture_time * 1000 / num_segments, 1);

figure();
for n = 1:num_segments
    xs = (sample_rate / segment_size * (-segment_size/2 : segment_size/2 - 1)) + center_frequency;
    ys = abs(fftshift(fft(rx_capture((n-1) * segment_size + 1 : n * segment_size))));

    subplot(segments(1), segments(2), n)
    plot(xs, ys)
    axis([center_frequency-max_doppler_shift  center_frequency+max_doppler_shift 0 y_axis_max])
    [t, s] = title(['Segment ',num2str(n),'/', num2str(num_segments)], [num2str((n-1) * segment_time), ' - ', num2str(n * segment_time), ' ms']);
    xlabel('f [Hz]') 
end

close all
format long

%% General
% Settings
file = "RX01"; % contains rx_capture (data), center_frequency [Hz], sample_rate [Hz], capture_time [s]
max_doppler_shift = 36e3; % [Hz]
segments = [6 20];
y_axis_max = 1000; % FFT only
num_segments_dominant = 120;

% Basic computation
load("Data\" + file + ".mat");
rx_capture = rx_capture.^4; % square twice
max_doppler_shift = max_doppler_shift * 1.1; % tolerance
singal_lenght = capture_time * sample_rate;
x_axis_min = center_frequency - max_doppler_shift;
x_axis_max = center_frequency + max_doppler_shift;

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
axis([x_axis_min x_axis_max 0 inf])
title ("FFT around center frequency (+- max Doppler shift)")
xlabel('f [Hz]')

%% Plot PSD
[p, f] = pspectrum(rx_capture, sample_rate);
figure();
plot(f, pow2db(p))
axis([min(f) max(f) min(pow2db(p)) max(pow2db(p))])
title("PSD")
xlabel('relative f [Hz]')
ylabel("PSD [dB]")

%% Plot FFT of sliced segments
num_segments = segments(1) * segments(2);
segment_size = floor(singal_lenght / num_segments);
segment_time = round(capture_time * 1000 / num_segments, 1);

figure();
for n = 1:num_segments
    xs = (sample_rate / segment_size * (-segment_size/2 : segment_size/2 - 1)) + center_frequency;
    ys = abs(fftshift(fft(rx_capture((n-1) * segment_size + 1: n * segment_size))));
    %ys = ys / max(ys(round(0.45 * segment_size): round(0.55 * segment_size))); % normalize

    subplot(segments(1), segments(2), n)
    plot(xs, ys)
    axis([x_axis_min x_axis_max 0 y_axis_max])
    title(['FFT ', num2str(n), '/', num2str(num_segments), ' (', num2str((n-1) * segment_time), ' - ', num2str(n * segment_time), ' ms)']);
end

%% Plot PSD of sliced segments
num_segments = segments(1) * segments(2);
segment_size = floor(singal_lenght / num_segments);
segment_time = round(capture_time * 1000 / num_segments, 1);

figure();
for n = 1:num_segments
    [p, f] = pspectrum(rx_capture((n-1) * segment_size + 1: n * segment_size), sample_rate);
    p = pow2db(p);

    subplot(segments(1), segments(2), n)
    plot(f, p)
    axis([-max_doppler_shift max_doppler_shift -60 -35])
    title(['PSD ', num2str(n), '/', num2str(num_segments) , ' (', num2str((n-1) * segment_time), ' - ', num2str(n * segment_time), ' ms)']);
end

%% Plot dominant frequency in time from FFT
num_segments = num_segments_dominant;

segment_size = floor(singal_lenght / num_segments);
segment_time = round(capture_time * 1000 / num_segments, 1);
dominant_frequencies = zeros(n, 1);
dominant_frequencies_timestamps = zeros(n, 1);
dominant_frequencies_strengths = zeros(n, 1);
for n = 1:num_segments
    xs = (sample_rate / segment_size * (-segment_size/2 : segment_size/2 - 1)) + center_frequency;
    ys = abs(fftshift(fft(rx_capture((n-1) * segment_size + 1: n * segment_size))));
    [argval, argidx] = max(ys); 
    dominant_frequencies_strengths(n) = argval; 
    dominant_frequencies(n) = xs(argidx);
    dominant_frequencies_timestamps(n) = (n-1) * segment_time;
end

figure(); hold on; grid on;
yyaxis left
plot(dominant_frequencies_timestamps, dominant_frequencies)
plot(dominant_frequencies_timestamps, ones(num_segments, 1) * center_frequency)
ylabel('f [Hz]')
yyaxis right
plot(dominant_frequencies_timestamps, dominant_frequencies_strengths)
ylabel('FFT amplitude')
xlabel('t [ms]') 
title(['Dominant frequency (FFT), ', num2str(num_segments), ' segments (', num2str(segment_time), ' ms each)']);

%% Plot dominant frequency in time from PSD
segment_size = floor(singal_lenght / num_segments);
segment_time = round(capture_time * 1000 / num_segments, 1);
dominant_frequencies = zeros(n, 1);
dominant_frequencies_timestamps = zeros(n, 1);
dominant_frequencies_strengths = zeros(n, 1);
for n = 1:num_segments
    [p, f] = pspectrum(rx_capture((n-1) * segment_size + 1: n * segment_size), sample_rate);
    p = pow2db(p);
    [argval, argidx] = max(p);
    dominant_frequencies_strengths(n) = argval;
    dominant_frequencies(n) = f(argidx);
    dominant_frequencies_timestamps(n) = (n-1) * segment_time;
end

figure(); hold on; grid on;
yyaxis left
plot(dominant_frequencies_timestamps, dominant_frequencies)
plot(dominant_frequencies_timestamps, zeros(num_segments, 1))
ylabel('relative f [Hz]')
yyaxis right
plot(dominant_frequencies_timestamps, dominant_frequencies_strengths)
ylabel('PSD [dB]')
xlabel('t [ms]') 
title(['Dominant frequency (PSD), ', num2str(num_segments), ' segments (', num2str(segment_time), ' ms each)']);

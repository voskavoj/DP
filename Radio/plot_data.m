close all
format long

%% Plot settings
plot_signal =       true    ;
plot_fft =          true    ;
plot_fft_cropped =  true    ;
plot_psd =          true    ;
plot_fft_seg =      true    ;
plot_psd_seg =      true    ;
plot_fft_dom =      true    ;
plot_psd_dom =      true    ;

%% General
% Settings
% file contains rx_capture (data), center_frequency [Hz], sample_rate [Hz], capture_time [s]
file_prefix =   "R_157" ; 
file_num =      "02"    ;
do_power =      true    ;
do_filter =     true    ;
segments =      [6 5]   ;
y_axis_max =    400     ; % FFT only
max_doppler_shift = 36e3; % [Hz]
num_segments_dominant = 120;

% Basic computation
load("Data\" + file_prefix + "_" + file_num + ".mat");
if do_power
    rx_capture = rx_capture.^4; % square twice
end
if do_filter
    lowpass(rx_capture, max_doppler_shift*2, sample_rate) % plot data
    rx_capture = lowpass(rx_capture, max_doppler_shift*2, sample_rate);
end
singal_lenght = capture_time * sample_rate;
x_axis_min = center_frequency - max_doppler_shift;
x_axis_max = center_frequency + max_doppler_shift;

%% Plot signal
if plot_signal
    t = (0 : singal_lenght - 1) / (singal_lenght - 1) * capture_time * 1000;
    figure();
    plot(t, rx_capture)
end

%% Create FFT and adjust axis
x = (sample_rate / singal_lenght * (-singal_lenght/2 : singal_lenght/2 - 1)) + center_frequency;
y = abs(fftshift(fft(rx_capture)));

%% Plot whole FFT
if plot_fft
    figure();
    plot(x, y)
    title("FFT")
    xlabel('f [Hz]')
end

%% Plot FFT around center frequency
if plot_fft_cropped
    figure();
    plot(x, y)
    axis([x_axis_min x_axis_max 0 inf])
    title ("FFT around center frequency (+- max Doppler shift)")
    xlabel('f [Hz]')
end

%% Plot PSD
if plot_psd
    [p, f] = pspectrum(rx_capture, sample_rate);
    figure();
    plot(f, pow2db(p))
    axis([min(f) max(f) min(pow2db(p)) max(pow2db(p))])
    title("PSD")
    xlabel('relative f [Hz]')
    ylabel("PSD [dB]")
end

%% Plot FFT of sliced segments
if plot_fft_seg
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
end

%% Plot PSD of sliced segments
if plot_psd_seg
    num_segments = segments(1) * segments(2);
    segment_size = floor(singal_lenght / num_segments);
    segment_time = round(capture_time * 1000 / num_segments, 1);
    
    figure();
    for n = 1:num_segments
        [p, f] = pspectrum(rx_capture((n-1) * segment_size + 1: n * segment_size), sample_rate);
        p = pow2db(p);
    
        subplot(segments(1), segments(2), n)
        plot(f, p)
        axis([-max_doppler_shift max_doppler_shift -60 0])
        title(['PSD ', num2str(n), '/', num2str(num_segments) , ' (', num2str((n-1) * segment_time), ' - ', num2str(n * segment_time), ' ms)']);
    end
end

%% Plot dominant frequency in time from FFT
if plot_fft_dom
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
end

%% Plot dominant frequency in time from PSD
if plot_psd_dom
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
end

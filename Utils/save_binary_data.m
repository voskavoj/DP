filename = "RX06"; % contains rx_capture (data)

load("Data\" + filename + ".mat");
% rx_capture = rx_capture.^4; % square twice
v = write_complex_binary(rx_capture, filename + "_binary");
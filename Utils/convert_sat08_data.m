format long;
t = readtable("Data\t00.txt", 'ReadVariableNames', false);
t.Properties.VariableNames = ["Epoch","Localtime","SatID","BeamID", "Lat", "Lon", "Alt", "Conf","Freq", "RawIQ"];
t(1,:)

steps = zeros(max(size(t)), 1);
for i = 2:max(size(t))
    steps(i) = t.Localtime(i) - t.Localtime(i-1);
end
% plot(t.Localtime, steps)
t.Localtime(end) / 1000 / 60
yyaxis left
plot(t.Freq)
yyaxis right
plot(t.SatID)

max(t.Freq) - min (t.Freq)


A = cell2mat(t.RawIQ(1));
B = str2num(A);
figure();
plot(abs(fftshift(fft(B))))

[p, f] = pspectrum(B);
figure();
plot(f, pow2db(p))

figure()
plot(abs(B))

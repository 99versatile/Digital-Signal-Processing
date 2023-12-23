yfile = readmatrix('./problem3/ycomplex_signal.csv');

time = yfile(:, 1);
y = yfile(:, 2);
Ts = time(2) - time(1);
Fs = 1/Ts;
Ns = 1024;
n = 0: 1: Ns-1;


M = 0;
T0 = 0;
D = 0;

figure(1)
% recieved signal y(t)
subplot(2,1,1)
stem(time, y)
xlim([0 time(Ns-1)])
% discrete sequence y(nTs)
subplot(2,1,2)
stem(n, y)
xlim([0 Ns])
%%
[ryy, lag] = xcorr(y);
[pksyy, locsyy] = findpeaks(abs(ryy), lag, 'MinPeakDistance', 90, 'MinPeakHeight', 2e-6);
figure(2)
hold on
stem(lag, abs(ryy))
stem(locsyy, pksyy)
xlabel('l')
ylabel('|ryy|')
hold off
%%
xfile = readmatrix("problem3/problem3_d/xcomplex_signal.csv");

xtime = xfile(:, 1);
x = xfile(:, 2);

[ryx, lagyx] = xcorr(y, x);
[pksyx, locsyx] = findpeaks(abs(ryx), lag, 'MinPeakDistance', 90, 'MinPeakHeight', 1e-6);
figure(3)
hold on
stem(lagyx, abs(ryx))
stem(locsyx, pksyx)
xlabel('l')
ylabel('|ryx|')
hold off
%%
mryy = mean(diff(locsyy));
vryy = var(diff(locsyy));
sryy = std(diff(locsyy));

mryx = mean(diff(locsyx));
vryx = var(diff(locsyx));
sryx = std(diff(locsyx));
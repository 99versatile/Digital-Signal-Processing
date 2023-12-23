N = 100;
n = 0: 1: N-1;
Ts = 10;
Fs = 1/Ts;

x = sin(2*pi*Fs*n);
figure(1)
stem(n, x)
xlabel('n')
ylabel('x[n]')

%% 
w1 = sqrt(60) .* transpose(rand(N, 1)-0.5);
w2 = sqrt(6/(10.^(1/10))) .* transpose(rand(N, 1)-0.5);
w3 = sqrt(3/5) .* transpose(rand(N, 1)-0.5);

figure(2)
subplot(3,1,1)
stem(n, w1)
xlabel('n')
ylabel('w1[n]')

subplot(3,1,2)
stem(n, w2)
xlabel('n')
ylabel('w2[n]')

subplot(3,1,3)
stem(n, w3)
xlabel('n')
ylabel('w3[n]')
%%
y1 = x + w1;
y2 = x + w2;
y3 = x + w3;

[ryy1, lag1] = xcorr(y1);
[ryy2, lag2] = xcorr(y2);
[ryy3, lag3] = xcorr(y3);
%%
[pks1, locs1] = findpeaks(ryy1, lag1, 'MinPeakDistance', 6, 'MinPeakHeight',0);
m1 = mean(diff(locs1));
var1 = var(diff(locs1));
std1 = std(diff(locs1));

figure(3)
subplot(2,1,1)
stem(n, y1)
xlabel('n')
ylabel('y1[n]')

subplot(2,1,2)
hold on
stem(lag1, ryy1)
stem(locs1, pks1)
xlabel('lag1 (l)')
ylabel('ryy1[l]')
hold off
%%
[pks2, locs2] = findpeaks(ryy2, lag2, 'MinPeakDistance', 6, 'MinPeakHeight',0);
m2 = mean(diff(locs2));
var2 = var(diff(locs2));
std2 = std(diff(locs2));

figure(4)
subplot(3,1,1)
stem(n, y2)
xlabel('n')
ylabel('y2[n]')

subplot(2,1,2)
hold on
stem(lag2, ryy2)
stem(locs2, pks2)
xlabel('lag2 (l)')
ylabel('ryy2[l]')
hold off
%%
[pks3, locs3] = findpeaks(ryy3, lag3, 'MinPeakDistance', 6, 'MinPeakHeight',0);
m3 = mean(diff(locs3));
var3 = var(diff(locs3));
std3 = var(diff(locs3));

figure(5)
subplot(2,1,1)
stem(n, y3)
xlabel('n')
ylabel('y3[n]')

subplot(2,1,2)
hold on
stem(lag3, ryy3)
stem(locs3, pks3)
xlabel('lag3 (l)')
ylabel('ryy3[l]')
hold off
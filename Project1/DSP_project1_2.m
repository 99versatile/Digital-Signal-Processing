file = readmatrix('./problem2/sunspot_number.csv');

time = file(:, 1);
ssn = file(:, 2);

stem(time, ssn)
xlabel('year')
ylabel('sun spots number')
xlim([1699 2023])

%%
[rssn, timelag] = xcorr(ssn);
[pksssn, locsssn] = findpeaks(rssn, timelag, 'MinPeakDistance', 6, 'MinPeakHeight',0);

figure(1)
hold on
stem(timelag, rssn)
stem(locsssn, pksssn)
xlabel('l')
ylabel('rssn')
hold off
ml = mean(diff(locsssn));
vl = var(diff(locsssn));
sl = std(diff(locsssn));

%%
m = mean(ssn);
nssn = ssn - m;

stem(time, nssn)
xlabel('year')
ylabel('sun spots number')
xlim([1699 2023])

%%
[rnssn, ntimelag] = xcorr(nssn);
[pksnssn, locsnssn] = findpeaks(rnssn, ntimelag, 'MinPeakDistance', 6, 'MinPeakHeight',0);

figure(2)
hold on 
stem(ntimelag, rnssn)
stem(locsnssn, pksnssn)
xlabel('l')
ylabel('rnssn')
hold off

mnl = mean(diff(locsnssn));
vnl = var(diff(locsnssn));
snl = std(diff(locsnssn));
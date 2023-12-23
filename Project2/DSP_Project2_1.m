%% Load original voice signal and plot its time series %%
load 10000.mat; Fs = 11025;% Sampling Rate of Signal
voice=voice-mean(voice);% Remove DC offset 
maxv=max(voice); time=1:length(voice); 
figure;plot(time/Fs,voice); grid on;
axis tight; xlabel('time'); title('Time series of voice');
%% Plot the frequency series of original voice signal %%
voice_len = length(voice);
f=1:voice_len/2;%Plotpositive frequency only 
fn=f/voice_len;%Normalization of the frequency 
figure;
vdft =fft(voice,voice_len); % Spectrum Computation
% vspec =vdft.*conj(vdft)/voice_len; %Power Spectrum 
semilogy(fn,abs(vspec(1:(length(voice)/2)))) % You may plot in log scale
% plot(fn,abs(vspec(1:(length(voice)/2)))) 
title('Frequency Spectrum of Voice');grid on;
%% Load Vocal Signal with one frequency component added %%
load voice_noise_1.mat; Fs = 11025;% Sampling Rate of Signal
signal_1=signal;
maxs_1=max(signal_1); time1=1:length(signal_1); 
figure;plot(time1/Fs,signal_1); grid on;
axis tight; xlabel('time'); title('Time series of Voice noise 1 signal');
%% Plot power spectrum of the vocal signal (log scale) %%
signal_len = length(signal_1);
f1=1:signal_len/2;%Plotpositive frequency only 
fn1=f1/signal_len;%Normalization of the frequency 
figure;
sdft_1 =fft(signal_1,signal_len); % Spectrum Computation
sspec_1 =sdft_1.*conj(sdft_1)/signal_len; %Power Spectrum 
semilogy(fn1,abs(sspec_1(1:(length(fn1))))) % You may plot in log scale
% plot(fn,abs(sspec_1(1:(length(fn))))) 
title('Frequency Spectrum of Voice noise 1 signal');grid on;
%% Load Vocal Signal with Two Frequency Components (Example 2) %%
load voice_noise_2.mat; Fs = 11025;
signal_2=signal;
maxs_2=max(signal_2); time2=1:length(signal_2);
figure;plot(time2/Fs, signal_2); grid on;
axis tight; xlabel('time'); title('Time series of Voice noise 2 signal');
%% Plot power spectrum of the vocal signal (log scale) %%
signal2_len = length(signal_2);
f2=1:signal2_len/2;%Plotpositive frequency only 
fn2=f2/signal2_len;%Normalization of the frequency 
figure;
sdft_2 =fft(signal_2,signal2_len); % Spectrum Computation
sspec_2 =sdft_2.*conj(sdft_2)/signal2_len; %Power Spectrum 
semilogy(fn2,abs(sspec_2(1:(length(fn2))))) % You may plot in log scale
% plot(fn,abs(sspec_1(1:(length(fn))))) 
title('Frequency Spectrum of Voice noise 2 signal');grid on;
%% Apply a Notch Filter to filter out the one frequency noise component %%
w0 = 0.2*pi;
FIR_Filter1 =[1 -2*cos(w0) 1]; 
Output_signal_1=conv(signal_1, FIR_Filter1, 'same'); 
odft_1 = fft(Output_signal_1, length(Output_signal_1));
ospec_1 =odft_1.*conj(odft_1)/length(Output_signal_1); %Power Spectrum 
figure; semilogy(fn1,abs(ospec_1(1:(length(fn1))))) % You may plot in log scale
%% plot the frequency response of the FIR notch filter %%
figure;
freqz(FIR_Filter1, 1)
%% plot the time series of FIR filter output %%
figure;
plot(time/Fs, Output_signal_1); grid on;
axis tight; xlabel('time'); title('Time series of voice');
%% Applying Notch Filter to filter out two frequency components %%
w1 = 0.3*pi;
FIR_Filter_temp = [1 -2*cos(w1) 1];
FIR_Filter2 = conv(FIR_Filter_temp, FIR_Filter1);
Output_signal_2 = conv(signal_2, FIR_Filter2, 'same');
odft_2 = fft(Output_signal_2, length(Output_signal_2));
ospec_2 =odft_2.*conj(odft_2)/length(Output_signal_2); %Power Spectrum 
figure; semilogy(fn2,abs(ospec_2(1:(length(fn2))))) % You may plot in log scale
%% plot the frequency response of the FIR notch filter %%
figure;
freqz(FIR_Filter2, 1)
%% plot the time series of FIR filter output %%
figure;
plot((1:length(signal_2))/Fs, Output_signal_2); grid on;
axis tight; xlabel('time'); title('Time series of voice');
%% Listen to the Output Signal %%
% a = audioplayer(Output_signal_1, Fs);
% play(a);
b = audioplayer(Output_signal_2, Fs);
play(b);
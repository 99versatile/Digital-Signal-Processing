%% Load the original Lenna image %%
Data = load("LENNA.MAT");
img = cell2mat(struct2cell(Data)); imagesc(img);
colormap gray
%% Plot the horizontal scan of image in row 200 %%
figure;
Row200oftheImg = img(200, :); %200th row of the image
plot(Row200oftheImg); 
grid on;
xlabel('Column Index');
title('Row 200 of the image');
%% Lowpass and Highpass Filtering %%
h1 = [0.5, 0.5]; % impulse response of lowpass filter
h2 = [0.5, -0.5]; % impulse response of highpass filter 

y1 = conv2(img, h1, 'same'); % apply LPF to the image by convolution
y2 = conv2(img, h2, 'same'); % apply HPF to the image by convolution
y3 = y1 + y2; % combine the outputs of LPF and HPF
%% Plot the Output of LPF %%
figure;
imagesc(y1); colormap gray; title('after h1 filter') 
%% Plot the Output of HPF %%
figure;
imagesc(y2); colormap gray; title('after h2 filter') 
%% Plot the combined output of LPF and HPF %%
figure;
imagesc(y3); colormap gray; title('after h1+h2 filter')
%% Load the horizontal noise-added image and compare with original image %%
Data1 = load("lenna_noise1.mat");
img1 = cell2mat(struct2cell(Data1)); imagesc(img1);
colormap gray

figure;
Row200oftheImg1 = img1(200, :); plot(Row200oftheImg1); hold on; plot(Row200oftheImg, 'r'); grid on;
xlabel('Column Index');
title('Row 200 of the image');
%% Compare the frequency spectrum plot between original and distorted image: horizontal scan of 200th row %%
N = 256;
Spectrum_temp0=fft(Row200oftheImg, N); % apply DTFT to 200th row of original image
Spectrum_temp1=fft(Row200oftheImg1, N); % apply DTFT to 200th row of noise added image
F = (0:N-1)/N; % Frequency scaling

figure;
semilogy(F, abs(Spectrum_temp0)); % logarithmic plot of the spectrum from original image
hold on; 
semilogy(F, abs(Spectrum_temp1), 'r'); % logarithmic plot of the spectrum from noisy image 
grid on;
xlim([0 .5]);
xlabel('Normalized Frequency');
%% Compare the frequency spectrum plot between original and distorted image: mean of image rows %%
RowsImg = mean(img, 1);
RowsImg1 = mean(img1, 1);
Spectrum_all0=fft(RowsImg, N);
Spectrum_all1=fft(RowsImg1, N);
F = (0:N-1)/N; % Frequency scaling

figure;
semilogy(F, abs(Spectrum_all0)); % logarithmic plot of the spectrum from original image
hold on; 
semilogy(F, abs(Spectrum_all1), 'r'); % logarithmic plot of the spectrum from noisy image 
grid on;
xlim([0 .5]);
xlabel('Normalized Frequency');
%% Generate and Plot the Running Average Filter with length 5 %%
L = 5; % Length of the running average window

RAF = ones(1, L)/L; % moving average filter coefficients

% Plotting the frequency response of length 5 RAF
freqz(RAF, 1);
%% Applying the length 5 RAF to the image and show the result: horizontal scan of 200th row %%
img_filt1 = conv2(img, RAF, 'same');
Row200oftheImgf1 = img_filt1(200, :);

figure;
plot(Row200oftheImg); hold on; plot(Row200oftheImgf1, 'r'); grid on;
legend('original image', 'filtered image', 'Location', 'northeast');
xlabel('Column Index');
title('Row 200 of the image');
%% 
Spectrum_f1 = fft(Row200oftheImgf1, N);

figure;
semilogy(F, abs(Spectrum_temp1)); % logarithmic plot of the spectrum from noisy image
hold on; 
semilogy(F, abs(Spectrum_f1), 'r'); % logarithmic plot of the spectrum from filtered image 
semilogy(F, abs(Spectrum_temp0), 'm'); % logarithmic plot of the spectrum from original image 
grid on;
xlim([0 .5]);
legend('noisy image', 'filtered image', 'original image', 'Location', 'northeast');
xlabel('Normalized Frequency');
%% Show the output image with lengh 5 RAF applied %%
figure;
imagesc(img_filt1); colormap gray; title('after RAF filter') 
%% Design and plot the notch filter without poles %%
w0 = 0.4*pi; % Notch frequency (normalized frequency in the range [0, pi])
H1_n = [1 -2*cos(w0) 1];
H1_d = 2-2*cos(w0);
freqz(H1_n, H1_d);
%% Apply the notch filter without poles %% 
img_filt2 = conv2(img, H1_n ./ H1_d, 'same');
Row200oftheImgf2 = img_filt2(200, :);
Spectrum_f2 = fft(Row200oftheImgf2, N);
%% plot the output: horizontal scan of 200th row %%
figure;
plot(Row200oftheImg); hold on;
plot(Row200oftheImgf2, 'r'); 
legend('original image', 'filtered image', 'Location', 'northeast');
xlabel('Column Index');
title('Row 200 of the image');
%% plot the output: frequency spectrum of row 200 %%
figure;
semilogy(F, abs(Spectrum_temp1)); hold on; 
semilogy(F, abs(Spectrum_f2), 'r'); 
semilogy(F, abs(Spectrum_temp0), 'm'); 
semilogy(F, abs(Spectrum_f1)); 
grid on; xlim([0 .5]);
legend('noisy image', 'filtered image', 'original image', 'RAF filtered', 'Location', 'northeast');
xlabel('Normalized Frequency');
title('Row 200 of the image');
%% plot the output image through notch filter without poles %%
figure;
imagesc(img_filt2); colormap gray; title('after FIR notch filter w/o poles') 
%% Design and plot the notch filter with poles %%
w0 = 0.4*pi; % Notch frequency (normalized frequency in the range [0, pi])
H2_n = [1 -2*cos(w0) 1];
H2_d = [1 -2*0.9*cos(w0) (0.9)^2];
freqz((1+(0.9)^2+2*(0.9)*cos(w0)).*H2_n, (2+2*cos(w0)).*H2_d);
%% Apply the notch filter with poles %%
filt3z = freqz((1+(0.9)^2+2*(0.9)*cos(w0)).*H2_n, (2+2*cos(w0)).*H2_d, size(img));
img_filt3z = conv2(img, real(ifft(filt3z)), 'same');
Row200oftheImgf3 = img_filt3z(200, :);
Spectrum_f3 = fft(Row200oftheImgf3, N);
%% plot the output: horizontal scan of 200th row %%
figure;
plot(Row200oftheImg); hold on;
plot(Row200oftheImgf3, 'r'); 
legend('original image', 'filtered image', 'Location', 'northeast');
xlabel('Column Index');
title('Row 200 of the image');
%% plot the output: frequency spectrum of 200th row %%
figure;
semilogy(F, abs(Spectrum_temp1)); hold on; 
semilogy(F, abs(Spectrum_f3), 'r'); 
semilogy(F, abs(Spectrum_temp0), 'm');
semilogy(F, abs(Spectrum_f1)); grid on; xlim([0 .5]);
legend('noisy image', 'filtered image', 'original image', 'RAF filtered', 'Location', 'northeast');
xlabel('Normalized Frequency');
title('Row 200 of the image');
%% plot the output image through notch filter without poles %%
figure;
imagesc(img_filt3z); colormap gray; title('after FIR notch filter with poles') 
%% Load and plot the unknown image %%
Data_u = load("unknown.MAT");
unknown_img = cell2mat(struct2cell(Data_u)); imagesc(unknown_img);
colormap gray
%% Plot the horizontal scan of row 200 of the unknown image %%
figure;
Row200oftheImg_u = unknown_img(200, :); plot(Row200oftheImg_u); grid on;
xlabel('Column Index');
title('Row 200 of the image');
%% Plot the power spectrum of the unknown image %%
% frequency spectrum of row 200 of distorted unknown image
Spectrum_u=fft(Row200oftheImg_u, N);
F = (0:N-1)/N; % Frequency scaling

figure;
semilogy(F, abs(Spectrum_u)); % logarithmic plot of the spectrum from original image
grid on;
xlim([0 .5]);
xlabel('Normalized Frequency');
%% Plot the power spectrum of the unknown image %%
% frequency spectrum of row average of distorted unknown image
RowsImg_u = mean(unknown_img, 1);
Spectrum_all_u=fft(RowsImg_u, N);
F = (0:N-1)/N; % Frequency scaling

figure;
semilogy(F, abs(Spectrum_all_u)); % logarithmic plot of the spectrum from original image
grid on;
xlim([0 .5]);
xlabel('Normalized Frequency');
%% Generate and Plot the Running Average Filter with length 10 %%
% RAF filter for restoring unknown image
L_u = 10; % Length of the running average window

RAF_u = ones(1, L_u)/L_u; % moving average filter coefficients
% Plotting the frequency response
freqz(RAF_u, 1);
%% Applying the RAF with length 10 and plot the output comparison: horizontal scan of row 200 %%
unknown_img_f = conv2(unknown_img, RAF_u, 'same');
Row200oftheImguf = unknown_img_f(200, :);

figure;
plot(Row200oftheImg_u); hold on; plot(Row200oftheImguf, 'r'); grid on;
xlabel('Column Index');
title('Row 200 of the image');
%% plot the output comparison: frequency spectrum of row 200 %%
Spectrum_u = fft(Row200oftheImg_u, N);
Spectrum_uf = fft(Row200oftheImguf, N);

figure;
semilogy(F, abs(Spectrum_u)); hold on; 
semilogy(F, abs(Spectrum_uf), 'r'); grid on; xlim([0 .5]);
legend('noisy image', 'filtered image', 'Location', 'northeast');
xlabel('Normalized Frequency');
title('Row 200 of the image');
%% Plot the image after RAF filtering %%
figure;
imagesc(unknown_img_f); colormap gray; title('after RAF filter'); 
%% Plot the unknown noise image and the horizontal scan of original and the noise image %%
Data_s = load("lenna_sand.MAT");
sand_img = cell2mat(struct2cell(Data_s)); imagesc(sand_img);
colormap gray

figure;
imagesc(sand_img); colormap gray; title('before RAF filter')

figure;
Row200oftheImg_s = sand_img(200, :); plot(Row200oftheImg_s);
hold on;
plot(Row200oftheImg, 'r');
grid on;
xlabel('Column Index');
title('Row 200 of the image');
%% The row-averaged power spectrum of original and noisy image %%
% frequency spectrum of row average of distorted unknown image
RowsImg_s = mean(sand_img, 1);
Spectrum_all_s=fft(RowsImg_s, N);
Spectrum_s=fft(Row200oftheImg_s, N);
F = (0:N-1)/N; % Frequency scaling

figure;
semilogy(F, abs(Spectrum_temp0)); % logarithmic plot of the spectrum from original image
hold on;
semilogy(F, abs(Spectrum_s), 'r'); % logarithmic plot of the spectrum from noise image
grid on;
legend('original image', 'unknown noise image', 'Location', 'northeast')
xlim([0 .5]);
xlabel('Normalized Frequency');
%% Generating the FIR Running Average Filter to block out random noise %%
% RAF filter for sand noise image
L_s = 20; % Length of the running average window

RAF_s = ones(1, L_s)/L_s; % moving average filter coefficients
% Plotting the frequency response
freqz(RAF_s, 1);
%% Apply the RAF to the noisy image and plot the horizontal scan of output image %%
sand_img_f = conv2(sand_img, RAF_s, 'same');
Row200oftheImgsf = sand_img_f(200, :);

figure;
plot(Row200oftheImg_s); hold on; plot(Row200oftheImgsf, 'r'); plot(Row200oftheImg, 'm'); grid on;
xlabel('Column Index');
title('Row 200 of the image');
%% frequency spectrum of row 200 of distorted unknown image %%
RowsImg_s = mean(sand_img, 1);
Spectrum_s=fft(Row200oftheImg_s, N);
Spectrum_sf=fft(Row200oftheImgsf, N);
F = (0:N-1)/N; % Frequency scaling

figure;
semilogy(F, abs(Spectrum_temp0)); % logarithmic plot of the spectrum from original image
hold on;
semilogy(F, abs(Spectrum_s), 'r'); % logarithmic plot of the spectrum from noise image
semilogy(F, abs(Spectrum_sf), 'm'); % logarithmic plot of the spectrum from filtered image
grid on;
legend('original image', 'unknown noise image', 'filtered image', 'Location', 'northeast')
xlim([0 .5]);
xlabel('Normalized Frequency');
%% Plot the output image from the FIR Running Average Filter %%
figure;
imagesc(sand_img_f); colormap gray; title('after RAF filter'); 
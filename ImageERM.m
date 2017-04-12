clc;

% Original Images are 256 pixels X 256 pixels.

fprintf('Reading Image \n');

% COMPRESSION_PERCENT = 0.1; % Compressed Signal will be approximately
% % n = 256^2 * COMPRESSION_PERCENT dimensional.
snr = .01;

Original_Image_Double = double(imread('Lena.bmp'));

fprintf('Performing Image Compression \n')

Compressed_Image_Double = fft(reshape(Original_Image_Double,[256*256,1]));
LSC = Compressed_Image_Double(abs(Compressed_Image_Double)<6*10^3);
Compressed_Image_Double(abs(Compressed_Image_Double)<6*10^3)=0;
Least_Significant_Indices = (Compressed_Image_Double(:) == 0);
Compressed_Image_Double(Least_Significant_Indices) = [];

n = 256*256-sum(Least_Significant_Indices);
N = 2*n+1000;
m = 1000;
L = [1:815];

f = Compressed_Image_Double;

fprintf('Creating Frames \n');

A = randn(N,2*n+m);
[A,~] = qr(A,0);

DF = sqrt(N/n)*A(:,1:n)';
EF = sqrt(n/N)*A(:,n+1:2*n)' + (n/N)*DF;
M = sqrt(N/m)*A(:,2*n+1:2*n+m)';

fprintf('Reconstructing Erasures \n');

FC = EF' * f;
FC(L) = zeros(size(L'));
f_R = DF*FC;

noise = randn(N,1);
noise = snr * noise ./ norm(noise) * norm(FC);
FC = FC + noise;

LC = setdiff(1:N,L);
FC(L) = -(M(:,L)' * M(:,L))\(M(:,L)' * (M(:,LC) * FC(LC)));

g = f_R + DF(:,L) * FC(L);

fprintf('Plotting Images \n');

C_f = zeros(256*256,1); % Compressed Image.
C_g = zeros(256*256,1); % Reconstructed Image.

C_f(Least_Significant_Indices(:) == 0) = f;
C_g(Least_Significant_Indices(:) == 0) = g;
C_g1(Least_Significant_Indices(:) == 0) = g;
C_g1(Least_Significant_Indices(:) == 1) = LSC;

Uncompressed_g = ifft(C_g);
Uncompressed_g = reshape(Uncompressed_g,[256,256]);

Uncompressed_f = ifft(C_f);
Uncompressed_f = reshape(Uncompressed_f,[256,256]);

Uncompressed_g1 = ifft(C_g1);
Uncompressed_g1 = reshape(Uncompressed_g1,[256,256]);

J_g = uint8(Uncompressed_g);
J_g1 = uint8(Uncompressed_g1);
J_f = uint8(Uncompressed_f);

figure(2);

subplot(2,2,1);
imshow(uint8(Original_Image_Double));
title('Original Image');

subplot(2,2,2);
imshow(J_f);
title('Compressed Image');

subplot(2,2,3);
imshow(J_g1);
title('Reconstructed & Uncompressed Image');

subplot(2,2,4);
imshow(J_g);
title('Compressed Image Reconstructed');
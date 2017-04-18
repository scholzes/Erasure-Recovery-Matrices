clc;

% Original Images are 256 pixels X 256 pixels.

fprintf('Reading Image \n');

COMPRESSION_PERCENT = 0.03; % Compressed Signal will be approximately
% n = 256^2 * COMPRESSION_PERCENT dimensional.
snr = .00;
percenterasures = .01;

Original_Image_Double = double(imread('Lena.bmp'));

fprintf('Performing Image Compression \n')

% Compressed_Image_Double = fft(reshape(Original_Image_Double,[256*256,1]));
% LSC = Compressed_Image_Double(abs(Compressed_Image_Double)<6*10^3);
% Compressed_Image_Double(abs(Compressed_Image_Double)<6*10^3)=0;
% Least_Significant_Indices = (Compressed_Image_Double(:) == 0);
% Compressed_Image_Double(Least_Significant_Indices) = [];

Compressed_Image_Double = fft(reshape(Original_Image_Double,[256*256,1]));
[S,I] = sort(abs(Compressed_Image_Double),'descend');
n = round(COMPRESSION_PERCENT*256*256)
LSC = Compressed_Image_Double(I(n+1:256*256));
Compressed_Image_Double(I(n+1:256*256)) = [];

N = 2*n+1000;
m = 1000;
L = [1:round(percenterasures*N)];

f = Compressed_Image_Double;

fprintf('Creating Frames \n');

A = randn(N,2*n+m);
[A,~] = qr(A,0);

DF = sqrt(N/n)*A(:,1:n)';
EF = sqrt(n/N)*A(:,n+1:2*n)' + (n/N)*DF;
M = sqrt(N/m)*A(:,2*n+1:2*n+m)';

fprintf('Creating More Frames \n');

A1 = randn(N,2*n+m);
A1(:,1:n) = DF';
[A1,~] = qr(A1,0);

DF1 = sqrt(N/n)*A1(:,1:n)';
EF1 = sqrt(n/N)*A1(:,n+1:2*n)' + (n/N)*DF1;
M1 = sqrt(N/m)*A1(:,2*n+1:2*n+m)';

fprintf('Reconstructing Erasures \n');

FC = EF' * f;
FC(L) = zeros(size(L'));
FC1 = FC;
f_R = DF*FC;

% noise = randn(N,1);
% noise = snr * noise ./ norm(noise) * norm(FC);
% FC = FC + noise;

LC = setdiff(1:N,L);
FC(L) = -(M(:,L)' * M(:,L))\(M(:,L)' * (M(:,LC) * FC(LC)));

g = f_R + DF(:,L) * FC(L);

fprintf('Reconstructing More Erasures \n');

FC1(L) = -(M1(:,L)' * M1(:,L))\(M1(:,L)' * (M1(:,LC) * FC1(LC)));

g2 = f_R + DF(:,L) * FC1(L);

fprintf('Plotting Images \n');

C_f = zeros(256*256,1); % Compressed Image.
C_g = zeros(256*256,1); % Reconstructed Image.
C_g2 = zeros(256*256,1);
C_g1 = zeros(256*256,1);

I1 = sort(I(1:n),'ascend');

C_f(I1) = f;
C_g(I1(1:n)) = g;
C_g2(I1(1:n)) = g2;
% C_g(I1(1:n)) = g1;
C_g1(I1(1:n)) = g;
C_g1(I(n+1:256*256)) = LSC;

Uncompressed_g = ifft(C_g);
Uncompressed_g = reshape(Uncompressed_g,[256,256]);

Uncompressed_g2 = ifft(C_g2);
Uncompressed_g2 = reshape(Uncompressed_g2,[256,256]);

Uncompressed_f = ifft(C_f);
Uncompressed_f = reshape(Uncompressed_f,[256,256]);

Uncompressed_g1 = ifft(C_g1);
Uncompressed_g1 = reshape(Uncompressed_g1,[256,256]);

J_g = uint8(Uncompressed_g);
J_g1 = uint8(Uncompressed_g1);
J_f = uint8(Uncompressed_f);
J_g2 = uint8(Uncompressed_g2);

figure;

subplot(2,2,1);
imshow(uint8(Original_Image_Double));
title('Original Image');

subplot(2,2,2);
imshow(J_f);
title('Compressed Image');

subplot(2,2,3);
imshow(J_g2);
title('Reconstructed with Wrong ERM');

subplot(2,2,4);
imshow(J_g);
title('Compressed Image Reconstructed');
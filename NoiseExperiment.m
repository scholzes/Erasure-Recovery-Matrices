clc;

% Original Images are 256 pixels X 256 pixels.

fprintf('Reading Image \n');

COMPRESSION_PERCENT = 0.15; % Compressed Signal will be approximately
% n = 256^2 * COMPRESSION_PERCENT dimensional.
percenterasures = .05;

Original_Image_Double = double(imread('Pepper.bmp'));

fprintf('Performing Image Compression \n')

Compressed_Image_Double = fft(reshape(Original_Image_Double,[256*256,1]));
[S,I] = sort(abs(Compressed_Image_Double),'descend');
n = round(COMPRESSION_PERCENT*256*256)
Compressed_Image_Double(I(n+1:256*256)) = [];

m = 1100;
N = 2*n+m;
L = [1:round(percenterasures*N)];
LC = setdiff(1:N,L);
snr = [.01, .03, .05, .07, .09];

f = Compressed_Image_Double;

fprintf('Creating Frames \n');

A = randn(N,2*n+m);
[A,~] = qr(A,0);

DF = sqrt(N/n)*A(:,1:n)';
EF = sqrt(n/N)*A(:,n+1:2*n)' + (n/N)*DF;
M = sqrt(N/m)*A(:,2*n+1:2*n+m)';

fprintf('Reconstructing Erasures \n');

figure;

C_f = zeros(256*256,1); % Compressed Image.
I1 = sort(I(1:n),'ascend');
C_f(I1) = f;
Uncompressed_f = ifft(C_f);
Uncompressed_f = reshape(Uncompressed_f,[256,256]);
J_f = uint8(Uncompressed_f);

subplot(3,5,3);
imshow(J_f);
title('Compressed Image');

for(j = 1:1:length(snr))

    FC = EF' * f;
    FC(L) = zeros(size(L'));
    noise = randn(size(LC'));
    noise = noise / norm(noise) * snr(j) * norm(FC(LC));
    FC(LC) = FC(LC) + noise;
    f_R = DF*FC;

    FC(L) = -(M(:,L)' * M(:,L)) \ (M(:,L)' * (M(:,LC) * FC(LC)));
    g = f_R + DF(:,L) * FC(L);

    C_f_R = zeros(256*256,1); % Erased Image with noise.
    C_f_R(I1(1:n)) = f_R;
    Uncompressed_f_R = ifft(C_f_R);
    Uncompressed_f_R = reshape(Uncompressed_f_R,[256,256]);
    J_f_R = uint8(Uncompressed_f_R);

    C_g = zeros(256*256,1); % Reconstructed Image.
    C_g(I1(1:n)) = g;
    Uncompressed_g = ifft(C_g);
    Uncompressed_g = reshape(Uncompressed_g,[256,256]);
    J_g = uint8(Uncompressed_g);

    subplot(3,5,5+j);
    imshow(J_f_R);
    title('Erased Image with Noise');

    subplot(3,5,10+j);
    imshow(J_g);
    title('Reconstructed with Noise');
    
    norm(noise)
    norm(f-f_R)

end

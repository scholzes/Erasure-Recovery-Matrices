clc;

% Original Images are 256 pixels X 256 pixels.

fprintf('Reading Image \n');

COMPRESSION_PERCENT = 0.15; % Compressed Signal will be approximately
% n = 256^2 * COMPRESSION_PERCENT dimensional.
% percenterasures = .01;

Original_Image_Double = double(imread('Lena.bmp'));

fprintf('Performing Image Compression \n')

Compressed_Image_Double = fft(reshape(Original_Image_Double,[256*256,1]));
[S,I] = sort(abs(Compressed_Image_Double),'descend');
n = round(COMPRESSION_PERCENT*256*256)
Compressed_Image_Double(I(n+1:256*256)) = [];

m = 1100;
N = 2*n+m;

f = Compressed_Image_Double;

fprintf('Creating Frames \n');

M = (1/sqrt(m)) * randn(m,N);
A = [M',randn(N,2*n)];
[A,~] = qr(A,0);
DF = A(:,m+1:m+n)';
EF = A(:,m+n+1:m+2*n)' + DF;

% A = randn(N,2*n+m);
% [A,~] = qr(A,0);
% 
% DF = sqrt(N/n)*A(:,1:n)';
% EF = sqrt(n/N)*A(:,n+1:2*n)' + (n/N)*DF;
% M = sqrt(N/m)*A(:,2*n+1:2*n+m)';

fprintf('Creating More Frames \n');

A = [DF',randn(N,n+m)];
[A,~] = qr(A,0);

DF1 = DF;
EF1 = A(:,n+m+1:2*n+m)' + DF1;
M1 = sqrt(N/m) * A(:,n+1:n+m)';

fprintf('Reconstructing Erasures \n');

percenterasures = [.01, .02, .03, .04, .05];

figure;

C_f = zeros(256*256,1); % Compressed Image.
I1 = sort(I(1:n),'ascend');
C_f(I1) = f;
Uncompressed_f = ifft(C_f);
Uncompressed_f = reshape(Uncompressed_f,[256,256]);
J_f = uint8(Uncompressed_f);

% subplot(3,5,2);
% imshow(J_f);
% title('Compressed Image');

Data = zeros(2,length(percenterasures));

for(j = 1:1:length(percenterasures))

    L = [1:round(percenterasures(j)*N)];

    FC = EF' * f;
    Data(1,j) = norm(FC(L));
    FC(L) = zeros(size(L'));
    FC1 = FC;
    f_R = DF*FC;

    LC = setdiff(1:N,L);
    FC(L) = -(M(:,L)' * M(:,L)) \ (M(:,L)' * (M(:,LC) * FC(LC)));
    g = f_R + DF(:,L) * FC(L);

    FC1(L) = -(M1(:,L)' * M1(:,L))\(M1(:,L)' * (M1(:,LC) * FC1(LC)));
    Data(2,j) = norm(FC1(L));

    g1 = f_R + DF(:,L) * FC1(L);
    
    C_f_R = zeros(256*256,1); % Reconstructed Image.
    C_f_R(I1(1:n)) = f_R;
    Uncompressed_f_R = ifft(C_f_R);
    Uncompressed_f_R = reshape(Uncompressed_f_R,[256,256]);
    J_f_R = uint8(Uncompressed_f_R);

    C_g = zeros(256*256,1); % Reconstructed Image.
    C_g(I1(1:n)) = g;
    Uncompressed_g = ifft(C_g);
    Uncompressed_g = reshape(Uncompressed_g,[256,256]);
    J_g = uint8(Uncompressed_g);

    C_g1 = zeros(256*256,1); % Reconstructed Image.
    C_g1(I1(1:n)) = g1;
    Uncompressed_g1 = ifft(C_g1);
    Uncompressed_g1 = reshape(Uncompressed_g1,[256,256]);
    J_g1 = uint8(Uncompressed_g1);

    subplot(3,5,j);
    imshow(J_f_R);
    title('No Reconstruction');
    
    subplot(3,5,5+j);
    imshow(J_g);
    title('Correct Erasure Recovery Matrix');

    subplot(3,5,10+j);
    imshow(J_g1);
    title('Incorrect Erasure Recovery Matrix');

end

Data

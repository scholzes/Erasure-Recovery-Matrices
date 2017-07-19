m = 250;
n = 250;
N = 1000;
snr = .05;
Trials = 50;
EC = [10:10:250];

P_Data = zeros(Trials,length(EC));
P_Data2 = zeros(Trials,length(EC));
DFP_Data = zeros(Trials,length(EC));
DFP_Data2 = zeros(Trials,length(EC));
Orig_Data = zeros(Trials,length(EC));
Orig_Data2 = zeros(Trials,length(EC));

for(k=1:1:length(EC))
    
    for(t = 1:1:Trials)
        
        L = [1:1:EC(k)];
        LC = setdiff(1:N,L);

        M = (1/sqrt(m)) * randn(m,N);
        A = [M',randn(N,2*n)];
        [A,~] = qr(A,0);
        F = A(:,m+1:m+n)';
	G = A(:,m+n+1:m+2*n)';
	G = F + G;

	P = M' * ((M * M') \ M);
	Gt = (eye(N) - P);
	T = sqrt(1/n) * randn(n,N);
	G1 = T * Gt;
	F1 = (G * G') \ G;

        f = randn(n,1);
        f = f./norm(f);

        noise_term = randn(length(LC),1);

	% Parseval Reconstruction:
        
        FC = F' * f;
        noise = snr * norm(FC(LC))/norm(noise_term) * noise_term;
        FC(LC) = FC(LC) + noise;
        FC(L) = zeros(size(L'));
        f_R = F*FC;
        
        LC = setdiff(1:N,L);
        FC(L) = -(M(:,L)' * M(:,L))\(M(:,L)' * (M(:,LC) * FC(LC)));
        
        g = f_R + F(:,L) * FC(L);
        
        P_Data(t,k) = norm(f-g);
        P_Data2(t,k) = norm(f-f_R);
	clearvars FC noise f_R g;

	% DFP Reconstruction

        FC = G' * f;
        noise = snr * norm(FC(LC))/norm(noise_term) * noise_term;
        FC(LC) = FC(LC) + noise;
        FC(L) = zeros(size(L'));
        f_R = F*FC;
        
        LC = setdiff(1:N,L);
        FC(L) = -(M(:,L)' * M(:,L))\(M(:,L)' * (M(:,LC) * FC(LC)));
        
        g = f_R + F(:,L) * FC(L);
        
        DFP_Data(t,k) = norm(f-g);
        DFP_Data2(t,k) = norm(f-f_R);
	clearvars FC noise f_R g;

	% Original Construction

        FC = G1' * f;
        noise = snr * norm(FC(LC))/norm(noise_term) * noise_term;
        FC(LC) = FC(LC) + noise;
        FC(L) = zeros(size(L'));
        f_R = F1 * FC;
        
        LC = setdiff(1:N,L);
        FC(L) = -(M(:,L)' * M(:,L))\(M(:,L)' * (M(:,LC) * FC(LC)));
        
        g = f_R + F1(:,L) * FC(L);
        
        Orig_Data(t,k) = norm(f-g);
        Orig_Data2(t,k) = norm(f-f_R);
	clearvars FC noise f_R g;
        
    end
    
    EC(k)

end

X = repmat(EC,Trials,1);
X = reshape(X,[length(EC)*Trials,1]);
P_Y = reshape(P_Data,[length(EC)*Trials,1]);
P_Z = reshape(P_Data2,[length(EC)*Trials,1]);
DFP_Y = reshape(DFP_Data,[length(EC)*Trials,1]);
DFP_Z = reshape(DFP_Data2,[length(EC)*Trials,1]);
Orig_Y = reshape(Orig_Data,[length(EC)*Trials,1]);
Orig_Z = reshape(Orig_Data2,[length(EC)*Trials,1]);

xlswrite('Temporary.xls',[X,P_Y,P_Z,DFP_Y,DFP_Z,Orig_Y,Orig_Z]);


% plot(X,Y,'x')
% hold on;
% plot(EC,median(Data));
% title('Erasure Set Size vs Reconstruction Error');
% xlabel('Erasure Set Size');
% ylabel('Reconstruction Error');
% hold off;

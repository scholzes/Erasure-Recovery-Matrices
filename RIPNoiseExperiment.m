m = 200;
n = 200;
N = 1000;
snr = .05;
Trials = 50;
EC = [5:5:20];

Data = zeros(Trials,length(EC));

for(k=1:1:length(EC))
    
    for(t = 1:1:Trials)
        
        L = [1:1:k];
        
        A = randn(N,2*n+m);
        
        [A,~] = qr(A,0);
        
        DF = sqrt(N/n)*A(:,1:n)';
        EF = sqrt(n/N)*A(:,n+1:2*n)' + (n/N)*DF;
        M = sqrt(N/m)*A(:,2*n+1:2*n+m)';
        
        f = randn(n,1);
        f = f./norm(f,2);
        
        FC = EF' * f;
        noise = randn(N,1);
        noise = snr/norm(noise) * noise;
        FC = FC + noise;
        FC(L) = zeros(size(L'));
        f_R = DF*FC;
        
        LC = setdiff(1:N,L);
        FC(L) = -(M(:,L)' * M(:,L))\(M(:,L)' * (M(:,LC) * FC(LC)));
        
        g = f_R + DF(:,L) * FC(L);
        
        Data(t,k) = norm(f-g);
        
    end
    
    k

end

X = repmat(EC,Trials);
X = reshape(X,[length(EC)*Trials,1]);
Y = reshape(Data,[length(EC)*Trials,1]);
plot(X,Y,'x')
hold on;
plot(EC,median(Data));
hold off;

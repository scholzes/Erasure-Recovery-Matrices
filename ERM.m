m = 250;
n = 250;
N = 1000;
snr = .05;
Trials = 50;

Data = zeros(Trials,length(5:5:75));

figure;
plot([0],[0],'o');
hold on;

for(k=[5:5:75])
    
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
        
        % Data(t,k/5) = norm(f-g);
        plot([k],[norm(f-g)],'x')
        
    end
    
    k

end

hold off;

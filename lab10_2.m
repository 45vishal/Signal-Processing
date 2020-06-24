clc;
clear all;
N=1024;
K=1024;
Q=4;
x2=zeros(K,N);

h_k=zeros(K,N);
X_s=zeros(K*N,1);
% Number of users
n=0:N-1;
for k=0:1:K-1
x1 = randi([0 Q-1],1,N);
X(k+1,:) = qammod(x1,Q,'UnitAveragePower',true);

end
s=0;
su=0;
a=1;
%%no of active users
for k=1:K
    x2(k,:)=ifft(X(k,:));

%     s=s+Xh_k;
   P=N;
X_s=par2ser(x2);
X_sf=fft(X_s,P);
 pow=(abs(X_sf).^2);
 pow1=sum(pow)/32;
peak=max(pow);

avg=mean(pow);
 Papr(k)=10*log10((k*peak/avg));
 a=a+1;
 thpapr(k)=10*log10(k);  
end
%%
plot(Papr);
xlabel('Number of OFDM signals');
ylabel('Papr_{dB}');
title('Papr as a function of K');
hold on 
plot(thpapr,'o');
legend('PracticalPapr','TheorticalPapr');
    function [output]=par2ser(input)
[m,n]=size(input);
for ii=1:m
   output((ii-1)*n+1:ii*n,1)=input(ii,:);
end
end
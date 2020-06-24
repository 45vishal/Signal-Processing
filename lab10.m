clc;
clear all
N=64;
K=128;
Q=64;
s=0;
P=16*N;
su=0;
x2=zeros(K,N);
P=16*N;
p=0:P-1;
wT=2*pi*p/P;
n=0:N-1;
wT=2*pi*p/P;
%H_k=zeros(K,N);
H_K=zeros(K,16*N);
X_K=zeros(K,16*N);
%K=S;
for k=1:K
x1 = randi([0 Q-1],1,N);
X(k,:) = qammod(x1,Q);

end
for k=1:K
for S= [6:26,38:57]

    
    if k==S
    x2(k,:)=ifft(X(k,:));
 %  X_K(k,:)=fft(X(k,:),P);
    end
end
end
X_s=par2ser(x2);
 X_sf=fft(X_s,P);
 pow=((abs(X_sf).^2));
A=20*log10((pow)/max(abs(pow)));



plot(wT/pi,A);
xlabel('Normalized frequency/pi');
ylabel('Power in dB');
title('Power spectrum of x(n)');
 hold on



P=16*N;
p=0:P-1;
wT=2*pi*p/P;
n=0:N-1;
X_s=par2ser(x2);
for k=1:K
    for S= [6:26,38:57]
        h_k(k,:)=(1/N)*exp(i*2*pi*k*n/N);
if k==S

H_k(k,:)=(h_k(k,:));
H_K(k,:)=fft(H_k(k,:),P);
end
su=su+(abs(H_K(k,:)).^2);
end
end

 X_sf=fft(X_s,P);
 pow=((abs(X_sf).^2))
A=20*log10((pow)/max(abs(pow)));
% plot(wT/pi,A);
% xlabel('Normalized frequency/pi');
% ylabel('Power in dB');
% title('Power spectrum of x(n)');
%  hold on
  ylim([-100 0]);

su=0;
p=0:P-1;
wT=2*pi*p/P;
n=0:N-1;

for k=1:K
    for S= [6:26,38:57]
h_k(k,:)=(1/N)*exp(i*2*pi*k*n/N);
if k==S

H_k(k,:)=(h_k(k,:));
H_K(k,:)=fft(H_k(k,:),P);
end
su=su+(abs(H_K(k,:)).^2);
end
end
 wT=linspace(0,2*pi-2*pi/P,P);
S=10*log10(abs(su)/max(abs(su)));
 plot(wT/pi,S);
legend('Power spectrum','Transmultiplexer');
function [output]=par2ser(input)
[m,n]=size(input);
for ii=1:m
   output((ii-1)*n+1:ii*n,1)=input(ii,:);
end
end

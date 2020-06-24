clc;
clear all;
N=8192;
%n=0:N-1;
sigma=0.01;
e1=sqrt(sigma/2)*(randn(1,N)+1i*randn(1,N));
e2=sqrt(sigma/2)*(randn(1,N)+1i*randn(1,N));
e1=e2;
%x=cos(wcT*n+sigma1*n)
wT=linspace(0,2*pi-2*pi/N,N);
wcT=0.125*pi;
delta=0.001;
for n=1:N-1
sigma1(1)=delta*e1(1);
sigma1(n+1)=sigma1(n)+delta*e1(n);
sigma2(1)=delta*e2(1);
sigma2(n+1)=sigma2(n)+delta*e2(n);

x(n+1)=cos(wcT*n+sigma1(n+1))-1i*sin(wcT*n+sigma2(n+1));
end
 w=window('blackmanharris',N)';
xn=x.*w;
X=fft(fftshift(xn,N));
n1=0:N-1;
A=20*log10(abs(X)/max(abs(X)));
figure();
f=1875e5;
f1=10:f;
semilogx(A);
% plot(wT/pi,A);
xlabel('Normalized frequency/pi');
ylabel('Magnitude response');
title('Amplitude spectrum of x(n) for e1(n)~=e2(n)');



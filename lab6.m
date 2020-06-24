clc;
clear all;
close all;
a=1;
for l=9:1:14
L=2.^l;

 N=300;
Q=16;
x = randi([0 Q-1],1,L+N);
x(1:N)=x(N+1:2*N);
l2=length(x);
n=(0:l2-1);
xn1 = qammod(x,Q,'UnitAveragePower',true);
xn=awgn(xn1,30);
%  scatterplot(xn);
%  title('QPSK signals with noise');
w0t=5*10e-5*pi;

alpha=0.1*pi;
n1=0:N-1;

yn=xn.*exp(-1i*(w0t*n+alpha));
est=angle(sum(conj(yn(n1+1)).*yn(n1+1+N)))/N;

Nu=sum(abs(xn1.^2));
De=sum(abs(yn-xn1).^2);
SDR(a)=10*log10(Nu/De);

%subplot(212);
%scatterplot(yn);

a=a+1;
end

%%
 l=12;
L=2.^l;
a=1;
m=1;

for N=10:20:400
    for g=1:100
Q=64;
x = randi([0 Q-1],1,L+N);
x(1:N)=x(N+1:2*N);
l2=length(x);
n=(0:l2-1);
xn1 = qammod(x,Q,'UnitAveragePower',true);
xn=awgn(xn1,30);
% scatterplot(xn);
% title('QPSK signals with noise');
w0t=5*10e-5*pi;

alpha=0.1*pi;
%pfo = comm.PhaseFrequencyOffset('PhaseOffset',30);
%yn1 = pfo(xn);
n1=0:N-1;

yn=xn.*exp(-1i*(w0t*n+alpha));
est=angle(sum(conj(yn(n1+1)).*yn(n1+1+N)))/N;
% scatterplot(yn);
% title('64Qpsk Signal with both FO and PO');
%n=length(yn);
ycomp=yn.*exp(-1i*est*n);
[heq, dmin,Err]=equalize1(xn1,ycomp,1);

ycomp1=ycomp(N+1:end);
ycomp2=conv(ycomp1,heq);
ycomp2=ycomp2(dmin+1:end);
Nu=sum(abs(xn1.^2));
De=sum(abs(yn-xn1).^2);
SDR(l)=10*log10(Nu/De);
 scatterplot(ycomp2);
% title('64_{qpsk} signal after compensation');
Nu2=sum(abs(xn1.^2));
De2=sum(abs(ycomp-xn1).^2);
SNDR1=10*log10(Nu2/De2);

Nu2=sum(abs(xn1.^2));
De2=sum(abs(ycomp2(1:L)-xn1(N+1:end)).^2);
SNDR(g)=10*log10(Nu2/De2);
a=a+1;
    end
    SNDRmean(m)=mean(SNDR);
    m=m+1;
end
n1=10:20:400;
plot(n1,SNDRmean);
xlabel('N: Number of redundant bits');
ylabel('SNDR_{dB}');
title('SNDR for different values of N');
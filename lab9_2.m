%% filter design
clc;
clear all;
S=20;
S1=30;
beta=1/3;
M1=4;
M2=8;
M=100;
Num2=0;
Den2=0;
Num1=0;
Den1=0;
fc=350e6;
N1=8192;

wc=2*pi*fc;
g1=rcosdesign(beta,S1,M1,'sqrt')/sqrt(M1);
G1=fft(g1,1000);
g2=rcosdesign(beta,S1,M2,'sqrt')/sqrt(M2);
%impz(G1);
N=512;

% signal x1 generation
fs=30e6;
L1=2^13;
L2=2^12;
A1=0.05;
A2=1;
Q = 16;
w1=-pi/3;
h1=M1*g1;
x1 = randi([0 Q-1],1,L1);
y1 = qammod(x1,Q,'UnitAveragePower',true);
y1=upsample(y1,M1);% upsample
y12=A1*conv(y1,h1,'same');
l1=length(y12);
n1=(0:l1-1);
a1=exp(1i*n1*w1);
y12=y12.*a1;


% signal x2 generation
%signal x2
w2=pi/6;
x2 = randi([0 Q-1],1,L2);
h2=M2*g2;
y2 = qammod(x2,Q,'UnitAveragePower',true);
y3=upsample(y2,M2);
y4=A2*conv(y3,h2,'same');
l2=length(y4);
n2=(0:l2-1);
a2=exp(1i*n2*w2);
y5=y4.*a2;


ytx=y12+y5; 
ytx1=[randi(1,1,N) ytx ];
ytx1(1:N)=ytx1(N+1:2*N);


%% passband generation
ytxn=upsample(ytx1,M);
len=length(ytxn);
t=(0:len-1)/fs;
glp=rcosdesign(beta,S,M,'sqrt')/sqrt(M);
hlp=M*glp;
ytxnhlp=conv(ytxn,hlp);
m=1:length(ytxnhlp);
b1=sqrt(2)*exp(1i*wc/(M*fs)*m);
ytxn2=ytxnhlp.*b1;
 yc=real(ytxn2);
%% channel
 cn=[1 zeros(1,15) 2.4 zeros(1,15) 1]/4; 
 yrc=conv(yc,cn,'same');
% freqz(cn);
%%
SNR=50;
yrc=awgn(yrc,SNR);
 a00=0.01; a11=1; a22=-0.2/max(abs(yrc));
 a33=0.15/max(abs(yrc).^2);
yrc=a00+a11*yrc+a22*yrc.^2+a33*yrc.^3;

%%
len1=1:length(yrc);
w0th=5*10e-7*pi;
alpha=0.1*pi;
yrc=yrc.*exp(-1i*(w0th*len1+alpha));
%% Band pass filter
[Nbp,fo,mo,w]=firpmord(fc*[0.5 0.9 1.1 1.5], [0 1 0], 0.001*[1 1 1], M*fs);
    Dbp=4*round(Nbp/8); 
    Nbp=2*Dbp;
    hbp=firpm(Nbp,fo,mo,w);
yrc=conv(yrc,hbp,'same');



%% Demodulation
g11=0.99;
 g22=1.02;
 %sigma1=0.5;
 %sigma2=0.55;
 n=0:length(m);
sigma=1;
e1=sqrt(sigma/2)*(randn(1,length(m))+1i*randn(1,length(m)));
e2=sqrt(sigma/2)*(randn(1,length(m))+1i*randn(1,length(m)));
e1=e2;
delta=0.0001;
for n=1:length(m)-1
sigma1(1)=delta*e1(1);
sigma1(n+1)=sigma1(n)+delta*e1(n);
sigma2(1)=delta*e2(1);
sigma2(n+1)=sigma2(n)+delta*e2(n);
b2(n+1)=sqrt(2)*(g11*cos(wc/(M*fs)*n+sigma1(n+1))-1i*g22*sin(wc/(M*fs)*n+sigma2(n+1)));
end
 %b2=sqrt(2)*exp(-1i*wc/(M*fs)*m);
 yrxn=yrc.*b2;
 yrxlp=conv(yrxn,glp);
yrx=downsample(yrxlp,M);
yrxd=yrx(S+1:end);

%% linearization
yrxd2 =sqrt(2)*(yrxd-3/4*a33*(abs(yrxd).^2).*(yrxd));
 %% Estimation and compensation of IQM 
% freqz(yrxd);
 yr=real(yrxd2);
 yi=imag(yrxd2);
 Pr=mean(yr.^2);
 Pi=mean(yi.^2);
 gest=sqrt(Pi/Pr);
 yp=-yr.*yi;
 Pyp=mean(yp);
 phiest=Pyp/Pr;
 % Compensation
 yr=yr*gest;
 yi=yi+yr*phiest;
 yrxd1=yr+j*yi;
%% CFO estimation and compensationn
n5=0:N-1;
est=angle(sum(conj(yrxd1(n5+1)).*yrxd1(n5+1+N)))/N;
j=0:length(yrxd1)-1;
ycomp=yrxd1.*exp(-1i*est*j);
ycomp=ycomp(N+1:end);
%% Equalization  
Neq=6;
 [heq, dmin,Err]=equalize1(ytx,ycomp,Neq);
  yrxeq=conv(ycomp,heq);
 yrxeq=yrxeq(dmin+1:end);
%  figure(2);
% freqz(heq);
%% detection of x2(n)
 %demod 

 j2=0:length(yrxeq)-1;

  a9=exp(-1i*j2*w2);
 yrx2=yrxeq.*a9;
 yrec12=conv(yrx2,g2,'same');
 yrec22=downsample(yrec12,M2);
 yrec32=yrec22/A2;
 xest2=qamdemod(yrec32,Q,'UnitAveragePower',true);
Num2=sum(abs(x2).^2);
Den2=sum(abs(xest2(1:L2)-x2).^2);
SIDR2=10*log10(Num2/Den2);
scatterplot(yrec32);
title('xest2');
%% Detection of signal x1(n)
 j1=0:length(yrxeq)-1;
 a10=exp(-1i*j1*w1);
 yrx1=yrxeq.*a10;
 yrec11=conv(yrx1,g1,'same');
 yrec21=downsample(yrec11,M1);
 yrec31=yrec21/A1;
 xest1=qamdemod(yrec31,Q,'UnitAveragePower',true);
Num1=sum(abs(x1).^2);
Den1=sum(abs(xest1(1:L1)-x1).^2);
SIDR1=10*log10(Num1/Den1);
%SIDR11(SNR,:)=SIDR1;
scatterplot(yrec31);
title('xest1');
%axis([-1.5 1.5 -1.5 1.5])
err1=symerr(xest1(1:L1),x1)
err2=symerr(xest2(1:L2),x2)
%plot(SIDR1);
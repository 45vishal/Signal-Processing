%% filter design for A1=1,A2=0.1
clc;
%clear all;
%S=5;
beta=1/3;
M1=4;
M2=8;
Num2=0;
Den2=0;
Num1=0;
Den1=0;
%SIDR1=size(40);
%SIDR2=size(40);
S=20;
g1=rcosdesign(beta,S,M1,'sqrt')/sqrt(M1);
%freqz(g1);
G1=fft(g1,1000);
g2=rcosdesign(beta,S,M2,'sqrt')/sqrt(M2);
G2=fft(g2,1000);
%impz(g2.*g2);
G1M=M1*G1.*G1;
%G1=downsample(M1*(g1.*g1),M1);
%impz(G1);


% signal x1 generation
fs=30e6;
L1=2^13;
L2=2^12;
A1=1;
A2=0.05;
Q = 64;
w1=-pi/3;
h1=M1*g1;
x1 = randi([0 Q-1],1,L1);
y1 = qammod(x1,Q,'UnitAveragePower',true);
y1=upsample(y1,M1);% upsample
y12=A1*conv(y1,h1);
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
y4=A2*conv(y3,h2);
l2=length(y4);
n2=(0:l2-1);
a2=exp(1i*n2*w2);
y5=y4.*a2;

append=zeros(1,length(y5)-length(y12));
y12=[y12 append];
 ytx=y12+y5;    %transmission of both signals
%yrx=awgn(ytx,100);
% detection of x2(n)
%ytx=ytx(S+1:end);
 %demod 
 
 n6=1:length(ytx);
  a9=exp(-1i*n2*w2);
 yrx2=ytx.*a9;

 yrec12=conv(yrx2,g2);
 yrec22=downsample(yrec12,M2);
 yrec32=yrec22/A2;
 xest2=qamdemod(yrec32,Q,'UnitAveragePower',true);

Num2=Num2+(abs(x2).^2);
%Num2=sum(Num2);
Den2=Den2+(abs(xest2(1:L2)-x2).^2);
SIDR2=10*log10(Num2/Den2);
%end2
%
% Detection of signal x1(n)
n3=1:length(ytx);
 a5=exp(-1i*n2*w1);
 yrx1=ytx.*a5;

 yrec11=conv(yrx1,g1);
 yrec21=downsample(yrec11,M1);
 yrec31=yrec21/A1;
 xest1=qamdemod(yrec31,Q,'UnitAveragePower',true);
Num1=Num1+(abs(x1).^2);
Den1=Den1+(abs(xest1(1:L1)-x1).^2);
SIDR1(S)=10*log10(Num1/Den1);
% find(xest1~=x1)
 [sym1,rate1]=symerr(xest1(1:L1),x1)
 [sym2,rate2]=symerr(xest2(1:L2),x2)
% figure;
% plot(SIDR2)
% %%
scatterplot(yrec31);
scatterplot(yrec32);
%scatterplot(yrec31);
%  yout=fftshift(fft(ytx,L2));
%  A=20*log10(abs(yout));
%  wT=linspace(0,2*pi-2*pi/L2,L2);
% 
%  figure(2);
%  subplot(311);
%   plot(wT/pi,A);
%  title('A1=1,A2=0.1 with Q=64');
%  xlabel('Normalized frequency,\omegaT/pi');
% ylabel('Amplitude spectrum');
 %stem(x1);


%% filter design Q=4
clc;
%clear all;
%S=5;
beta=1/3;
M1=4;
M2=8;
Num2=0;
Den2=0;
Num1=0;
Den1=0;
%SIDR1=size(40);
%SIDR2=size(40);
S=20
g1=rcosdesign(beta,S,M1,'sqrt')/sqrt(M1);
%freqz(g1);
G1=fft(g1,1000);
g2=rcosdesign(beta,S,M2,'sqrt')/sqrt(M2);
G2=fft(g2,1000);
%impz(g2.*g2);
G1M=M1*G1.*G1;
%G1=downsample(M1*(g1.*g1),M1);
%impz(G1);


% signal x1 generation
fs=30e6;
L1=2^13;
L2=2^12;
A1=1;
A2=0.01;
Q = 64;
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


 yrx=y5+y12;    %transmission of both signals
%yrx=awgn(ytx,3);
% detection of x2(n)
 %demod 
  a3=exp(-1i*n2*w2);
 yrx2=yrx.*a3;
 

 yrec12=conv(yrx2,g2,'same');
 yrec22=downsample(yrec12,M2);
 yrec32=yrec22/A2;
 xest2=qamdemod(yrec32,Q,'UnitAveragePower',true);

Num2=Num2+(abs(x2).^2);
%Num2=sum(Num2);
Den2=Den2+(abs(xest2-x2).^2);
%SIDR2(S)=10*log10(Num2/Den2);
%end2
%
% Detection of signal x1(n)
 a3=exp(-1i*n1*w1);
 yrx1=yrx.*a3;

 yrec11=conv(yrx1,g1,'same');
 yrec21=downsample(yrec11,M1);
 yrec31=yrec21/A1;
 xest1=qamdemod(yrec31,Q,'UnitAveragePower',true);
Num1=Num1+(abs(x1).^2);
Den1=Den1+(abs(xest1-x1).^2);
%SIDR1(S)=10*log10(Num1/Den1);
[sym1,rate1]=symerr(xest1,x1);
[sym2,rate2]=symerr(xest2,x2);
% figure;
% plot(SIDR2)
% %%
%scatterplot(y1);
%scatterplot(yrec32);
 yout=fftshift(fft(ytx,L2));
 A=20*log10(abs(yout));
 wT=linspace(0,2*pi-2*pi/L2,L2);
e=symerr(xest1,x1)


 %stem(x1);
 
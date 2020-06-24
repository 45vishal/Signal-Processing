clc;
clear all;
close all;
N=8192;
n=0:N-1;
p=5;
w1=(0.25)*pi;
a=[(2*pi*p)/N];
delta=a;
 x=sin(w1*n)+0.01*sin((w1+delta).*n);
w=window('blackmanharris',N)'; %The “”’ makes w a row vector like x

ww=fft(w);
xw=x.*w;
wT=linspace(0,2*pi-2*pi/N,N);
Xw=(fft(xw));
A=20*log10(abs(Xw)/max(abs(Xw)));

figure();

subplot(411);
title('Amplitude 0.01');
plot(wT/pi,A);
title('DFT with coherent sampling and Blackmanharris');
xlabel('Normalized frequency,\omegaT/pi');
ylabel('Amplitude spectrum');
%%
p=5;
w1=(0.25)*pi;
a=[(2*pi*p)/N];
delta=a;
 x=sin(w1*n)+0.001*sin((w1+delta).*n);
w=window('rectwin',N)'; %The “”’ makes w a row vector like x

ww=fft(w);
xw=x%.*w;
wT=linspace(0,2*pi-2*pi/N,N);
Xw=(fft(xw));
A=20*log10(abs(Xw)/max(abs(Xw)));

%figure;
subplot(412);
plot(wT/pi,A);
title('DFT with Coherent sampling and rectangular window');
xlabel('Normalized frequency,\omegaT/pi');
ylabel('Amplitude spectrum');
%%
p=50;
w1=(0.25+1/N)*pi;
a=[(2*pi*p)/N];
delta=a;
 x=sin(w1*n)+0.001*sin((w1+delta).*n);
w=window('rectwin',N)'; %The “”’ makes w a row vector like x

ww=fft(w);
xw=x%.*w;
wT=linspace(0,2*pi-2*pi/N,N);
Xw=(fft(xw));
A=20*log10(abs(Xw)/max(abs(Xw)));

%figure;
subplot(413);
plot(wT/pi,A);
title('DFT with non-coherent sampling and rectangular window');
xlabel('Normalized frequency,\omegaT/pi');
ylabel('Amplitude spectrum');
%%
p=10;
w1=(0.25)*pi;
a=[(2*pi*p)/N];
delta=a;
 x=sin(w1*n)+0.001*sin((w1+delta).*n);
w=window('blackmanharris',N)'; %The “”’ makes w a row vector like x

ww=fft(w);
xw=x.*w;
wT=linspace(0,2*pi-2*pi/N,N);
Xw=(fft(xw));
A=20*log10(abs(Xw)/max(abs(Xw)));

%figure;
subplot(414);
plot(wT/pi,A);
title('DFT with non coherent sampling and blackmanharris window');
xlabel('Normalized frequency,\omegaT/pi');
ylabel('Amplitude spectrum');
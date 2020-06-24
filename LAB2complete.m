

%% for the first expression of x
N = 2.^20;
x1 = randn(1,N);
%  % p = 1000;
% %x3 = cos(2*pi*p/N);
% %x4 = cos(0.25*pi);

B1 = 1;
y1 = AGCunity(x1);
xq1 = quant(y1,B1);
e1 = xq1-x1;
figure(1)
title('Histogram plots x(n) = randn(1, N) ');
subplot(231)
histogram(e1);
title('B=1');
B2 = 2;
y1 = AGCunity(x1);
xq2 = quant(y1,B2);
e2 = xq2-x1;
subplot(232)
histogram(e2);
title('B=2');
B3 = 4;
y1 = AGCunity(x1);
xq3 = quant(y1,B3);
e3 = xq3-x1;
subplot(233)
histogram(e3);
title('B=4');
B4 = 8;
y1 = AGCunity(x1);
xq4 = quant(y1,B4);
e4 = xq4-x1;
subplot(234)
histogram(e4);
title('B=8');
B5 = 16;
y1 = AGCunity(x1);
xq5 = quant(y1,B5);
e5 = xq5-x1;
subplot(235)
histogram(e5);
title('B=16');
sgtitle('Histogram plots for x(n)=randn(1,N)');
%% for the second expression of x 
N=2^20;
sigma=1/12;

x2 = sqrt(sigma)*rand(1,N)-1/2;
y2 = AGCunity(x2);

B1 = 1;
xq1 = quant(y2,B1);
e1 = xq1-x2;

figure(2)

subplot(231)
histogram(e1);
title('B=1');
B2 = 2;
xq2 = quant(y2,B2);
e2 = xq2-x2;
subplot(232)
histogram(e2);
title('B=2');

B3 = 4;
xq3 = quant(y2,B3);
e3 = xq3-x2;
subplot(233)
histogram(e3);
title('B=4');
B4 = 8;
xq4 = quant(y2,B4);
e4 = xq4-x2;
subplot(234)
histogram(e4);
title('B=8');
B5 = 16;
xq5 = quant(y2,B5);
e5 = xq5-x2;
subplot(235)
histogram(e5);
title('B=16');
sgtitle('Histogram plots for x(n) rand(1, N)-1/2');
%% for the third expession of x
N = 2.^20;
n = 0:N-1;
p = 3;
% p = 4;

WT = (2* pi* p)/N;
x3 = cos(WT*n);
y3 = x3;

B1 = 1;
xq1 = quant(y3,B1);
e1 = xq1-x3;
figure(10)
subplot(231)
histogram(e1);
title('B=1');
B2 = 2;
xq2 = quant(y3,B2);
e2 = xq2-x3;
subplot(232)
histogram(e2);
title('B=2');
B3 = 4;
xq3 = quant(y3,B3);
e3 = xq3-x3;
subplot(233)
histogram(e3);
title('B=4');
B4 = 8;
xq4 = quant(y3,B4);
e4 = xq4-x3;
subplot(234)
histogram(e4);
title('B=8');
B5 = 16;
xq5 = quant(y3,B5);
e5 = xq5-x3;
subplot(235)
histogram(e5);
title('B=16');
sgtitle('Histogram plots for x(n)=cos(2*pi*p/N)n')
%% for the fourth expression of x 
N = 2.^20;
n = 0:N-1;
x4 = cos(0.25*pi*n);

y4=x4;

B1 = 1;
xq1 = quant(y4,B1);
e1 = xq1-x4;

figure(4)
subplot(231)
histogram(e1);
title('B=1');
B2 = 2;
xq2 = quant(y4,B2);
e2 = xq2-x4;
subplot(232)
histogram(e2);
title('B=2');
B3 = 4;
xq3 = quant(y4,B3);
e3 = xq3-x4;
subplot(233)
histogram(e3);
title('B=4');
B4 = 8;
xq4 = quant(y4,B4);
e4 = xq4-x4;
subplot(234)
histogram(e4);
title('B=8');
B5 = 16;
xq5 = quant(y4,B5);
e5 = xq5-x4;
subplot(235)
histogram(e5);
title('B=16');
sgtitle('Histogram plots for x(n)=cos(0.25*pi*n)');

%%  PART2
a=1;
for l=5:1:15
    
N = 2.^l;
P1=0;
P2=0;
n = 0:N-1;
B = 16;

    seq2 = cos(0.5*sqrt(3)*pi*(n));
    xq = quant(seq2,B);
    e = xq - seq2;
    P1 = sum((seq2.^2));
    P2 = sum((e.^2));
    SNR(l) = 10*log10(P1/P2);
    a=a+1;
end

 figure(5);
 subplot(311);
 plot(SNR)
 %subplot(212);
 xlabel('Values of 2^N different samples');
 ylabel('SNR_{dB}');
 title('SNR values for x(n)=cos(0.5*sqrt(3)*pi*n)');
 
% ylim([0 150]);

a=1;
for l=5:1:15
N = 2.^l;
P1=0;
P2=0;
n = 0:N-1;
B = 16;
%for i=1:N-1
    x = cos(0.25*pi*(n));
    seq2=x;
    xq = quant(seq2,B);
    e = xq - seq2;
    P1 = P1+(seq2.^2);
    P2 = P2+(e.^2);
    SNR(l) = 10*log10(P1/P2);
a=a+1;
end

subplot(312);
plot(SNR);
 xlabel('Values of 2^N different samples');
 ylabel('SNR_{dB}');
 title('SNR for the signal cos(0.25*pi*n)');


for l=5:1:15
N = 2.^l;
P1=0;
P2=0;
n = 0:N-1;
B = 16;
%for i=1:N-1
    seq3 = cos(0.5*pi*(n));
    xq = quant(seq3,B);
    e = xq - seq3;
    P1 = P1+(seq3.^2);
    P2 = P2+(e.^2);
    SNR(l) = 10*log10(P1/P2);
end

subplot(313);
plot(SNR);

 xlabel('Values of 2^N different samples');
 ylabel('SNR_{dB}');
 title('SNR for the signal cos(0.5pi*n)');
 sgtitle('SNR value as a function of N');
%% Part-3
N = 2.^13;
n = 0:N-1;


WT = 0.25*pi;
x3 = cos(WT*n);
plot(x3);
y3 =AGCunity( x3);
B1 = 16;
xq1 = quant(y3,B1);
plot(xq1);
wT=linspace(0,2*pi-2*pi/N,N);
 w=window('blackmanharris',N)';
    xw=xq1.*w;
    Xw=fft(xw,N);
    A=20*log10(2*abs(Xw)/N);
    
     w1=window('rectwin',N)';
    xw1=xq1.*w1;
    Xw1=fft(xw1,N);
    A1=20*log10(2*abs(Xw1)./N);

figure(1);
subplot(211);
plot(wT/pi,A);
xlabel('Normalized frequency/pi');
ylabel('Magnitude_{dB}');
title('Blackmanharris window');
subplot(212);
plot(wT/pi,A1);
xlabel('Normalized frequency/pi');
ylabel('Magnitude_{dB}');
title('rectangular window');
sgtitle('Amplitude response for coherent sampling');

%%
N=8192;
n=0:N-1;
WT = (0.25+1/N)*pi;
x3 = cos(WT*n);

y3 = AGCunity(x3);
B1 = 16;
xq1 = quant(y3,B1);
 w=window('blackmanharris',N)';
    xw=xq1.*w;
    Xw=fft(xw,N);
    A2=20*log10(2*abs(Xw)/N);
 w=window('rectwin',N)';
    xw=xq1.*w;
    Xw=fft(xw,N);
    A3=20*log10(2*abs(Xw)/N);
figure(2);
subplot(211);
plot(A2);
xlabel('Normalized frequency/pi');
ylabel('Magnitude_{dB}');
title('Blackmanharris window');
subplot(212);
plot(A3);
xlabel('Normalized frequency/pi');
ylabel('Magnitude_{dB}');
title('Rectangular window');
sgtitle('Amplitude response for non-coherent sampling');



% e1 = xq1-x3;
% figure(8)
% subplot(231)
% histogram(e1);



%% 
%PART 4
N=8192;
n=0:N-1;
WT = (0.5*sqrt(3)*pi);
x3 = cos(WT*n);
wT=linspace(0,2*pi-2*pi/N,N);
% plot(x3);
% hold on;
y3 = AGCunity(x3);
plot(y3);
%legend('x3','y3');
B1 = 16;
xq1 = quant(y3,B1);
 w=window('blackmanharris',N)';
    xw=xq1.*w;
    Xw=fft(xw,N);
    A1=20*log10(2*abs(Xw)/N);

figure;
subplot(211);
plot(wT/pi,A1);
xlabel('Normalized frequency/pi');
ylabel('Magnitude_{dB}');
title('N=8192');
N=8192*16;
n=0:N-1;
wT=linspace(0,2*pi-2*pi/N,N);
WT = (0.5*sqrt(3)*pi);
x3 = cos(WT*n);

y3 = AGCunity(x3);
B1 = 16;
xq1 = quant(y3,B1);
 w=window('blackmanharris',N)';
    xw=xq1.*w;
    Xw=fft(xw,N);
    A1=20*log10(2*abs(Xw)/N);


subplot(212);
plot(wT/pi,A1);
xlabel('Normalized frequency/pi');
ylabel('Magnitude_{dB}');
title('N=16*8192');
sgtitle('Quantized spectra for N=8192 and N=16*8192');

% e1 = xq1-x3;
% figure(8)
% subplot(231)
% histogram(e1);

%%
function y=AGCunity(x)
N=length(x);
n=0:N-1;
t=length(n);
gain_level=1;
output_power_Normal = 10.^(gain_level/10);
% Noisy signal
%9Find the frequency components of the signal using serdes.AGC.
%p_db=10*log10(sum(x(1:t,:).^2)/t);
energy=sum(x(:,1:t).^2);
k1=sqrt((output_power_Normal*N)/energy);
y(1,1:N)= x(1,1:N)*k1;
y= normalize(y, 'range', [-1 1]);
return
end


%%
function xq=quant(x,B)
%%
% B=16;
% N=256;
% n=0:N-1;
% WT = (0.25+1/N)*pi;
% x = cos(WT*n);
Nsamples=length(x);
% plot(x);
% hold on;
d=1/(2^(B-1));
xmax=max(x)-d/2;
xmin=min(x)+d/2;
%quant1=zeros(1:Nsamples-1,1);
for i=xmin:d:xmax
    for j=1:Nsamples
        if((i-d/2)<=x(j)&&(x(j)<=i+d/2))
            quant1(j)=i;
        end
    end 
end
xq = quant1;
% hold on;
% stem(quant1);
return
end
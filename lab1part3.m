%%
delphase=0.1;
A1=1;
A2=1;
fc1=10e6;
fc2=30.5e6;
delta=1/8;
fs=40e6;
T=1/fs;
N=8192;
xa1_filt=0;
xa2_filt=0;
x0=0;
xnsum=0;
fs=40e6;
 wc=0.8125*pi*fs;
 n=0:N-1;
 
for P=1:200%Filter Order

t =(0:N-1)*T;
f=(0:N-1)*fs;
fs=40e6;
n=0:N-1;
 wc=0.8125*pi*fs;
% Amaxln=1
  Amax=0.1;
     %Amax=mag2db(Amaxln);
     e=sqrt(10^(0.1*Amax)-1);
     
 for    k=-5:5
         
     f1k=fc1+(fc1*k*delta);
     f2k=fc2+(fc1*k*delta);
      w1=2*pi*f1k;
      w2=2*pi*f2k;
   H1=sqrt(1/(1+e^2*((w1/wc)^(2*P))));
   H2=sqrt(1/(1+e^2*((w2/wc)^(2*P))));
      xa1_filt=xa1_filt+(H1*sin(2*pi.*f1k*t+(delphase*(-1)^k)));
      xa2_filt=xa2_filt+(H2.*sin(2*pi.*f2k*t));
      xa1=xa1+sin(2*pi.*f1k*t);
     xa2=xa2+sin(2*pi.*f2k*t);
 end
 
 

 %summation of signals
      xn= (A1*xa1_filt)+(A2*xa2_filt);
      x0 = sum((A1*xa1).^2);
fx= xn-xa1;
xnsum=sum((fx).^2);
xa1_filt=0;
 xa2_filt=0;
 xa1=0;
 xa2=0;
%SIDR
    SIDR(P)=10*log10(x0./xnsum);
end
% figure;
%  plot(H11);
 w=window('blackmanharris',N)';
    xw=xn.*w;
    Xw=fft(xw,N);
   % freqz(Xw);
    A=20*log10(abs(Xw));
     figure(3);
     subplot(411);
 plot(wT/pi,A);
  title('DFT with blackmanharris N=8192 with A1=1 ,Amax=0.1dB and phase distortion=0.1');
 xlabel('Normalized frequency,\omegaT/pi');
ylabel('Amplitude spectrum');
    figure(4);  
    subplot(411);
 plot(SIDR)
 title('SIDR with A1=1 and Amax=0.1dB and phase distortion=0.1');
 xlabel('Order of Filter:P');
ylabel('SIDR(dB)');
%%
delphase=0.1;
A1=0.01;
A2=1;
fc1=10e6;
fc2=30.5e6;
delta=1/8;
fs=40e6;
T=1/fs;
N=8192;
xa1_filt=0;
xa2_filt=0;
x0=0;
xnsum=0;
fs=40e6;
 wc=0.8125*pi*fs;
 n=0:N-1;
 
for P=1:200%Filter Order

t =(0:N-1)*T;
f=(0:N-1)*fs;
fs=40e6;
n=0:N-1;
 wc=0.8125*pi*fs;
% Amaxln=1
  Amax=0.1;
     %Amax=mag2db(Amaxln);
     e=sqrt(10^(0.1*Amax)-1);
     
 for    k=-5:5
         
     f1k=fc1+(fc1*k*delta);
     f2k=fc2+(fc1*k*delta);
      w1=2*pi*f1k;
      w2=2*pi*f2k;
   H1=sqrt(1/(1+e^2*((w1/wc)^(2*P))));
   H2=sqrt(1/(1+e^2*((w2/wc)^(2*P))));
      xa1_filt=xa1_filt+(H1*sin(2*pi.*f1k*t+(delphase*(-1)^k)));
      xa2_filt=xa2_filt+(H2.*sin(2*pi.*f2k*t));
      xa1=xa1+sin(2*pi.*f1k*t);
     xa2=A2*xa2+sin(2*pi.*f2k*t);
 end
 
 

 %summation of signals
      xn= (A1*xa1_filt)+(A2*xa2_filt);
      x0 = sum((A1*xa1).^2);
fx= xn-xa1;
xnsum=sum((fx).^2);
xa1_filt=0;
 xa2_filt=0;
 xa1=0;
 xa2=0;
%SIDR
    SIDR(P)=10*log10(x0./xnsum);
end
% figure;
%  plot(H11);
 w=window('blackmanharris',N)';
    xw=xn.*w;
    Xw=fft(xw,N);
   % freqz(Xw);
    A=20*log10(abs(Xw));
     figure(3);
     subplot(412);
 plot(wT/pi,A);
  title('DFT with blackmanharris N=8192 with A1=0.01 ,Amax=0.1dB and phase distortion=0.1 ');
 xlabel('Normalized frequency,\omegaT/pi');
ylabel('Amplitude spectrum');
    figure(4);  
    subplot(412);
 plot(SIDR)
 title('SIDR with A1=0.01 and Amax=0.1dB and phase distortion=0.1');
 xlabel('Order of Filter:P');
ylabel('SIDR(dB)');
%%
delphase=0.1;
A1=1;
A2=1;
fc1=10e6;
fc2=30.5e6;
delta=1/8;
fs=40e6;
T=1/fs;
N=8192;
xa1_filt=0;
xa2_filt=0;
x0=0;
xnsum=0;
fs=40e6;
 wc=0.8125*pi*fs;
 n=0:N-1;
 
for P=1:200%Filter Order

t =(0:N-1)*T;
f=(0:N-1)*fs;
fs=40e6;
n=0:N-1;
 wc=0.8125*pi*fs;
% Amaxln=1
  Amax=0.01;
     %Amax=mag2db(Amaxln);
     e=sqrt(10^(0.1*Amax)-1);
     
 for    k=-5:5
         
     f1k=fc1+(fc1*k*delta);
     f2k=fc2+(fc1*k*delta);
      w1=2*pi*f1k;
      w2=2*pi*f2k;
   H1=sqrt(1/(1+e^2*((w1/wc)^(2*P))));
   H2=sqrt(1/(1+e^2*((w2/wc)^(2*P))));
      xa1_filt=xa1_filt+(H1*sin(2*pi.*f1k*t+(delphase*(-1)^k)));
      xa2_filt=xa2_filt+(H2.*sin(2*pi.*f2k*t));
      xa1=A1*xa1+sin(2*pi.*f1k*t);
     xa2=A2*xa2+sin(2*pi.*f2k*t);
 end
 
 

 %summation of signals
      xn= (A1*xa1_filt)+(A2*xa2_filt);
      x0 = sum((A1*xa1).^2);
fx= xn-xa1;
xnsum=sum((fx).^2);
xa1_filt=0;
 xa2_filt=0;
 xa1=0;
 xa2=0;
%SIDR
    SIDR(P)=10*log10(x0./xnsum);
end
% figure;
%  plot(H11);
 w=window('blackmanharris',N)';
    xw=xn.*w;
    Xw=fft(xw,N);
   % freqz(Xw);
    A=20*log10(abs(Xw));
     figure(3);
     subplot(413);
 plot(wT/pi,A);
  title('DFT with blackmanharris N=8192 with A1=1 ,Amax=0.01dB and phase distortion=0.1');
 xlabel('Normalized frequency,\omegaT/pi');
ylabel('Amplitude spectrum');
    figure(4);  
    subplot(413);
 plot(SIDR)
 title('SIDR with A1=1 and Amax=0.01dB and phase distortion=0.1');
 xlabel('Order of Filter:P');
ylabel('SIDR(dB)');
%%
delphase=0.1;
A1=0.01;
A2=1;
fc1=10e6;
fc2=30.5e6;
delta=1/8;
fs=40e6;
T=1/fs;
N=8192;
xa1_filt=0;
xa2_filt=0;
x0=0;
xnsum=0;
fs=40e6;
 wc=0.8125*pi*fs;
 n=0:N-1;
 
for P=1:200%Filter Order

t =(0:N-1)*T;
f=(0:N-1)*fs;
fs=40e6;
n=0:N-1;
 wc=0.8125*pi*fs;
% Amaxln=1
  Amax=0.01;
     %Amax=mag2db(Amaxln);
     e=sqrt(10^(0.1*Amax)-1);
     
 for    k=-5:5
         
     f1k=fc1+(fc1*k*delta);
     f2k=fc2+(fc1*k*delta);
      w1=2*pi*f1k;
      w2=2*pi*f2k;
   H1=sqrt(1/(1+e^2*((w1/wc)^(2*P))));
   H2=sqrt(1/(1+e^2*((w2/wc)^(2*P))));
      xa1_filt=xa1_filt+(H1*sin(2*pi.*f1k*t+(delphase*(-1)^k)));
      xa2_filt=xa2_filt+(H2.*sin(2*pi.*f2k*t));
      xa1=xa1+sin(2*pi.*f1k*t);
     xa2=xa2+sin(2*pi.*f2k*t);
 end
 
 

 %summation of signals
      xn= (A1*xa1_filt)+(A2*xa2_filt);
      x0 = sum((A1*xa1).^2);
fx= xn-xa1;
xnsum=sum((fx).^2);
xa1_filt=0;
 xa2_filt=0;
 xa1=0;
 xa2=0;
%SIDR
    SIDR(P)=10*log10(x0./xnsum);
end
% figure;
%  plot(H11);
 w=window('blackmanharris',N)';
    xw=xn.*w;
    Xw=fft(xw,N);
   % freqz(Xw);
    A=20*log10(abs(Xw)/max(abs(Xw)));
     figure(3);
     subplot(414);
 plot(wT/pi,A);
  title('DFT with blackmanharris N=8192 with A1=0.01 ,Amax=0.01dB and phase distortion=0.1 ');
 xlabel('Normalized frequency,\omegaT/pi');
ylabel('Amplitude spectrum');
    figure(4);  
    subplot(414);
 plot(SIDR)
 title('SIDR with A1=0.01 and Amax=0.01dB and phase distortion=0.1');
 xlabel('Order of Filter:P');
ylabel('SIDR(dB)');
%% dephase = 0.01
delphase=0.01;
A1=1;
A2=1;
fc1=10e6;
fc2=30.5e6;
delta=1/8;
fs=40e6;
T=1/fs;
N=8192;
xa1_filt=0;
xa2_filt=0;
x0=0;
xnsum=0;
fs=40e6;
 wc=0.8125*pi*fs;
 n=0:N-1;
 
for P=1:200%Filter Order

t =(0:N-1)*T;
f=(0:N-1)*fs;
fs=40e6;
n=0:N-1;
 wc=0.8125*pi*fs;
% Amaxln=1
  Amax=0.1;
     %Amax=mag2db(Amaxln);
     e=sqrt(10^(0.1*Amax)-1);
     
 for    k=-5:5
         
     f1k=fc1+(fc1*k*delta);
     f2k=fc2+(fc1*k*delta);
      w1=2*pi*f1k;
      w2=2*pi*f2k;
   H1=sqrt(1/(1+e^2*((w1/wc)^(2*P))));
   H2=sqrt(1/(1+e^2*((w2/wc)^(2*P))));
      xa1_filt=xa1_filt+(H1*sin(2*pi.*f1k*t+(delphase*(-1)^k)));
      xa2_filt=xa2_filt+(H2.*sin(2*pi.*f2k*t));
      xa1=xa1+sin(2*pi.*f1k*t);
     xa2=xa2+sin(2*pi.*f2k*t);
 end
 
 

 %summation of signals
      xn= (A1*xa1_filt)+(A2*xa2_filt);
      x0 = sum((A1*xa1).^2);
fx= xn-xa1;
xnsum=sum((fx).^2);
xa1_filt=0;
 xa2_filt=0;
 xa1=0;
 xa2=0;
%SIDR
    SIDR(P)=10*log10(x0./xnsum);
end
% figure;
%  plot(H11);
 w=window('blackmanharris',N)';
    xw=xn.*w;
    Xw=fft(xw,N);
   % freqz(Xw);
    A=20*log10(abs(Xw));
     figure(5);
     subplot(411);
 plot(wT/pi,A);
  title('DFT with blackmanharris N=8192 with A1=1 ,Amax=0.1dB and phase distortion=0.01 ');
 xlabel('Normalized frequency,\omegaT/pi');
ylabel('Amplitude spectrum');
    figure(6);  
    subplot(411);
 plot(SIDR)
 title('SIDR with A1=1 and Amax=0.1dB and phase distortion=0.01');
 xlabel('Order of Filter:P');
ylabel('SIDR(dB)');
%%
delphase=0.01;
A1=0.01;
A2=1;
fc1=10e6;
fc2=30.5e6;
delta=1/8;
fs=40e6;
T=1/fs;
N=8192;
xa1_filt=0;
xa2_filt=0;
x0=0;
xnsum=0;
fs=40e6;
 wc=0.8125*pi*fs;
 n=0:N-1;
 
for P=1:200%Filter Order

t =(0:N-1)*T;
f=(0:N-1)*fs;
fs=40e6;
n=0:N-1;
 wc=0.8125*pi*fs;
% Amaxln=1
  Amax=0.1;
     %Amax=mag2db(Amaxln);
     e=sqrt(10^(0.1*Amax)-1);
     
 for    k=-5:5
         
     f1k=fc1+(fc1*k*delta);
     f2k=fc2+(fc1*k*delta);
      w1=2*pi*f1k;
      w2=2*pi*f2k;
   H1=sqrt(1/(1+e^2*((w1/wc)^(2*P))));
   H2=sqrt(1/(1+e^2*((w2/wc)^(2*P))));
      xa1_filt=xa1_filt+(H1*sin(2*pi.*f1k*t+(delphase*(-1)^k)));
      xa2_filt=xa2_filt+(H2.*sin(2*pi.*f2k*t));
      xa1=A1*xa1+sin(2*pi.*f1k*t);
     xa2=A2*xa2+sin(2*pi.*f2k*t);
 end
 
 

 %summation of signals
      xn= (A1*xa1_filt)+(A2*xa2_filt);
      x0 = sum((A1*xa1).^2);
fx= xn-xa1;
xnsum=sum((fx).^2);
xa1_filt=0;
 xa2_filt=0;
 xa1=0;
 xa2=0;
%SIDR
    SIDR(P)=10*log10(x0./xnsum);
end
% figure;
%  plot(H11);
 w=window('blackmanharris',N)';
    xw=xn.*w;
    Xw=fft(xw,N);
   % freqz(Xw);
    A=20*log10(abs(Xw));
     figure(5);
     subplot(412);
 plot(wT/pi,A);
  title('DFT with blackmanharris N=8192 with A1=0.01 ,Amax=0.1dB and phase distortion=0.01 ');
 xlabel('Normalized frequency,\omegaT/pi');
ylabel('Amplitude spectrum');
    figure(6);  
    subplot(412);
 plot(SIDR)
 title('SIDR with A1=0.01 and Amax=0.1dB and phase distortion=0.01');
 xlabel('Order of Filter:P');
ylabel('SIDR(dB)');
%%
delphase=0.01;
A1=1;
A2=1;
fc1=10e6;
fc2=30.5e6;
delta=1/8;
fs=40e6;
T=1/fs;
N=8192;
xa1_filt=0;
xa2_filt=0;
x0=0;
xnsum=0;
fs=40e6;
 wc=0.8125*pi*fs;
 n=0:N-1;
 
for P=1:200%Filter Order

t =(0:N-1)*T;
f=(0:N-1)*fs;
fs=40e6;
n=0:N-1;
 wc=0.8125*pi*fs;
% Amaxln=1
  Amax=0.01;
     %Amax=mag2db(Amaxln);
     e=sqrt(10^(0.1*Amax)-1);
     
 for    k=-5:5
         
     f1k=fc1+(fc1*k*delta);
     f2k=fc2+(fc1*k*delta);
      w1=2*pi*f1k;
      w2=2*pi*f2k;
   H1=sqrt(1/(1+e^2*((w1/wc)^(2*P))));
   H2=sqrt(1/(1+e^2*((w2/wc)^(2*P))));
      xa1_filt=xa1_filt+(H1*sin(2*pi.*f1k*t+(delphase*(-1)^k)));
      xa2_filt=xa2_filt+(H2.*sin(2*pi.*f2k*t));
      xa1=A1*xa1+sin(2*pi.*f1k*t);
     xa2=A2*xa2+sin(2*pi.*f2k*t);
 end
 
 

 %summation of signals
      xn= (A1*xa1_filt)+(A2*xa2_filt);
      x0 = sum((A1*xa1).^2);
fx= xn-xa1;
xnsum=sum((fx).^2);
xa1_filt=0;
 xa2_filt=0;
 xa1=0;
 xa2=0;
%SIDR
    SIDR(P)=10*log10(x0./xnsum);
end
% figure;
%  plot(H11);
 w=window('blackmanharris',N)';
    xw=xn.*w;
    Xw=fft(xw,N);
   % freqz(Xw);
    A=20*log10(abs(Xw));
     figure(5);
     subplot(413);
 plot(wT/pi,A);
  title('DFT with blackmanharris N=8192 with A1=1 ,Amax=0.01dB and phase distortion=0.01 ');
 xlabel('Normalized frequency,\omegaT/pi');
ylabel('Amplitude spectrum');
    figure(6);  
    subplot(413);
 plot(SIDR)
 title('SIDR with A1=1 and Amax=0.01dB and phase distortion=0.01');
 xlabel('Order of Filter:P');
ylabel('SIDR(dB)');
%%
delphase=0.01;
A1=0.01;
A2=1;
fc1=10e6;
fc2=30.5e6;
delta=1/8;
fs=40e6;
T=1/fs;
N=8192;
xa1_filt=0;
xa2_filt=0;
x0=0;
xnsum=0;
fs=40e6;
 wc=0.8125*pi*fs;
 n=0:N-1;
 
for P=1:200%Filter Order

t =(0:N-1)*T;
f=(0:N-1)*fs;
fs=40e6;
n=0:N-1;
 wc=0.8125*pi*fs;
% Amaxln=1
  Amax=0.01;
     %Amax=mag2db(Amaxln);
     e=sqrt(10^(0.1*Amax)-1);
     
 for    k=-5:5
         
     f1k=fc1+(fc1*k*delta);
     f2k=fc2+(fc1*k*delta);
      w1=2*pi*f1k;
      w2=2*pi*f2k;
   H1=sqrt(1/(1+e^2*((w1/wc)^(2*P))));
   H2=sqrt(1/(1+e^2*((w2/wc)^(2*P))));
      xa1_filt=xa1_filt+(H1*sin(2*pi.*f1k*t+(delphase*(-1)^k)));
      xa2_filt=xa2_filt+(H2.*sin(2*pi.*f2k*t));
      xa1=xa1+sin(2*pi.*f1k*t);
     xa2=xa2+sin(2*pi.*f2k*t);
 end
 
 

 %summation of signals
      xn= (A1*xa1_filt)+(A2*xa2_filt);
      
      x0 = sum((A1*xa1).^2);
fx= xn-xa1;
xnsum=sum((fx).^2);
xa1_filt=0;
 xa2_filt=0;
 xa1=0;
 xa2=0;
%SIDR
    SIDR(P)=10*log10(x0./xnsum);
end
% figure;
%  plot(H11);
 w=window('blackmanharris',N)';

    xw=xn.*w;
    Xw=fft(xw,N);
   % freqz(Xw);
    A=20*log10(abs(Xw));
     figure(5);
     subplot(414);
 plot(wT/pi,A);
  title('DFT with blackmanharris N=8192 with A1=0.01 ,Amax=0.01dB and phase distortion=0.01 ');
 xlabel('Normalized frequency,\omegaT/pi');
ylabel('Amplitude spectrum');
    figure(6);  
    subplot(414);
 plot(SIDR)
 title('SIDR with A1=0.01 and Amax=0.01dB and phase distortion=0.01');
 xlabel('Order of Filter:P');
ylabel('SIDR(dB)');



clc; clear ;close all;

 L1=2.^12;

x = randn(1,L1);
cn=0.25*exp(0.1i*pi)*[1 2.4 1];
 figure(1);
freqz(cn);
title('Magnitude response of the channel');

ytx=conv(x,cn); 
yrx=awgn(ytx,30);
 Neq=20;
 [heq, dmin,Er]=equalize1(x,yrx,Neq);
 E(Neq)=Er;
 
 figure(2);
 freqz(heq);

 legend('Frequency response of the equalizer');
 title('y(tx) is the  white gaussian noise');

%%
L1=2.^12;
% n=0:L1-1;
x = randn(1,L1);
cn=0.25*exp(0.1i*pi)*[1 2.4 1];

[Nhlp,fo,mo,w]=firpmord([0.45 0.55],[1 0],[0.01 0.001],2);
hlp=firpm(2*round(Nhlp/2),fo,mo,w);

ytx=conv(x,hlp);
ytx=conv(ytx,cn); 
yrx=awgn(ytx,30);
 Neq=20
   [heq, dmin,Er]=equalize1(ytx,yrx,Neq);
 
 figure(3);
 freqz(heq);

title('ytx is a filtered white gaussian noise');


%%
 L1=2.^12;
% n=0:L1-1;
x = randn(1,L1);
cn=0.25*exp(0.1i*pi)*[1 2.4 1];

[Nhlp,fo,mo,w]=firpmord([0.45 0.55],[1 0],[0.01 0.001],2);
hlp=firpm(2*round(Nhlp/2),fo,mo,w);

ytx=conv(x,hlp);
ytx=conv(ytx,cn); 
yrx=awgn(ytx,30);
for Neq=1:50
   [heq, dmin,Er]=equalize1(ytx,yrx,Neq);
 E(Neq)=Er;
end
 figure(4);
 plot(E);
xlabel('Order of filter: Neq');
ylabel('Err(min)');
title('ytx is a filtered white gaussian noise');
legend('Err_{min}');

%%

 L1=2.^12;

x = randn(1,L1);
cn=0.25*exp(0.1i*pi)*[1 2.4 1];
 

ytx=conv(x,cn); 
yrx=awgn(ytx,30);
for  Neq=1:50
 [heq, dmin,Er]=equalize1(x,yrx,Neq);
 E(Neq)=Er;
end
 figure(4);
 plot(E);
xlabel('Order of filter: Neq');
ylabel('Err(min)');
title('ytx is a  white gaussian noise');
legend('Err_{min}');


function [heq, dmin,Err]=equalize1(ytx,yrx,Neq)

%%

a=1;
Err=[];

%First, you need to check the lengths of the input signals
 L=length(ytx);
if length(yrx)<L
    yrx(end+1:L)=0;
elseif length(yrx)>L
    yrx=yrx(1:L);
end


%A is a Toeplitz matrix 
part1 = [yrx zeros(1,Neq)].';  
part2 = [yrx(1) zeros(1,Neq)].';

Atop = toeplitz(part1, part2);


ERR=[];
h=[];
for delay=0:Neq 
    ytx_delayed=[zeros(1,delay) ytx zeros(1,Neq-delay)].';
    heq_save(delay+1,:) = (Atop'*Atop)\Atop'*ytx_delayed; % Eqaution 7 from Lab Memo
    

yrxeq=conv(yrx,heq_save(delay+1,:));
 Er(a)=sum(abs(yrxeq-ytx_delayed.').^2);
  a=a+1;

 
[~, index] = min(Er);
end
%freqz(yrxeq);
heq = heq_save(:,index); % return filter

dmin = index-1; % return dmin
Err=Er(dmin+1);
 return

end


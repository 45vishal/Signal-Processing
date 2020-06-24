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
part1 = [yrx zeros(1,Neq)].'; %you need to put zeros here as well they were missing in your code
part2 = [yrx(1) zeros(1,Neq)].';

Atop = toeplitz(part1, part2);


ERR=[];
h=[];
for delay=0:Neq % it is also OK to go through 0:2*Neq_i as well but you need to modify x_del accordingly 
    ytx_delayed=[zeros(1,delay) ytx zeros(1,Neq-delay)].';
    heq_save(delay+1,:) = (Atop'*Atop)\Atop'*ytx_delayed; % Eqaution 7 from Lab Memo
    %You need to calculate error here as you have done down below
    %.....

yrxeq=conv(yrx,heq_save(delay+1,:));
%yrxeq=[yrxeq append];
 Er(a)=sum(abs(yrxeq-ytx_delayed.').^2);
%  Err(a)=sum(Er);
  a=a+1;

 %dmin=find(Err==dmin1);
[dmin, index] = min(Er);
end
heq = heq_save(:,index); % return filter
dmin = index-1; % return dmin
Err=Er(dmin+1);
 return

end
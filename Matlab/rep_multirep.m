function [vcos,vsin]=rep_multirep(MAS,m,n,i,ksi)
%m- number of harmonics
%n- number of molecules


u0=MAS(:,1); %к-ты при единицах
u=MAS(:,2:(m+1));%к-ты при cos
v=MAS(:, (m+2):(2*m+1)); %к-ты при sin
 

ui=cat(2, u0(i), u(i,:));
vi=cat(2,0,v(i,:));
uk=cat(2,u0(ksi), u(ksi,:));
vk=cat(2,0,v(ksi,:));

B1=ui'*uk;
B2=ui'*vk;
B3=vi'*uk;
B4=vi'*vk;
%______________________________________________for B1&B4
for s=1:(m+1)
    k1(s)=sum(diag(flipdim(B1,2), m+1-s)) - sum(diag(flipdim(B4,2), m+1-s)); %minus because of minus in formulae
end;
k11(1)=trace(B1)+trace(B4);
for s=2:(m+1)
    k11(s)=sum(diag(B1+B1',s-1)) + sum(diag(B1+B1',s-1));
end;
k1=(k1+k11)/2;%vector for cos
%______________________________________________for B2&B3
for s=1:(m+1)
    k2(s)=sum(diag(flipdim(B2,2), m+1-s)) + sum(diag(flipdim(B3,2), m+1-s)) ; %because of symmetry 
end;
k22(1)=trace(B2) + trace(B3);
for s=1:(m+1)
    k22(s)= sum(diag(-B2+B2',s-1)) + sum(diag(-B3'+B3,s-1));
end;
k2=(k2+k22)/2;%vector for sin
%-----------------------------------------------

% vectors 1x3, containing fourier coefficients

vcos=k1;
vsin=k2;













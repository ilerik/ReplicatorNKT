function [RESMAS, f] = rep_ffunc(MAS,m,n,A,d,g)
%m- number of harmonics
%n- number of molecules
%A=n x n - replicator matrix
%q= coefficient matrix is in A

u0=MAS(:,1); %к-ты при единицах
u=MAS(:,1:(m+1));%к-ты при cos
v=MAS(:, (m+2):(2*m+1)); 
z=zeros(n,1);
v=cat(2,z,v);%к-ты при sin

vcos(:,:,:)=zeros(n,n,m+1);
vsin(:,:,:)=zeros(n,n,m+1);

for i=1:n %calculating the 3-dimentional arrays of cross-
    for k=1:n
        if A(i,k)~=0
            [vc vs]=rep_multirep(MAS,m,n,i,k);
            vcos(i,k,:)=vc;
            vsin(i,k,:)=vs;
        else
            vcos(i,k,:)=zeros(1,m+1);
            vsin(i,k,:)=zeros(1,m+1);
        end;
    end;
end;
% f=sum(sum(A.*vcos(:,:,1),1))*2*pi;%for normal replicator equation
f=exp(-sum(u0));%for exponential variation

UU=zeros(n,m+1);
VV=zeros(n,m+1);
for i=1:n
%     i
%     A(i,:);
    c(:,:)=vcos(i,:,:);
    s(:,:)=vsin(i,:,:);
%     A(i,:)*c
%     A(i,:)*s
%     v
    vector_i=(0:m).^2;
%     UU(i,:)=A(i,:)*c - f*u(i,:)- d(i)*vector_i.*u(i,:);%for normal
%     replicator equation
%     VV(i,:)=A(i,:)*s - f*v(i,:)- d(i)*vector_i.*v(i,:);
    UU(i,:)=A(i,:)*c*f - g(i)*u(i,:)- d(i)*vector_i.*u(i,:);%for exponential var.
    VV(i,:)=A(i,:)*s*f - g(i)*v(i,:)- d(i)*vector_i.*v(i,:);
end;

RESMAS(:,1:m+1)=UU;
RESMAS(:,(m+2):(2*m+1))=VV(:,2:m+1);

%RESMAS;
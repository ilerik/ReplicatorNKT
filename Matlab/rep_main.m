function rep_main

m=30;% the number of harmonics
n=3;% the number of molecules
A=[1.2 0.5   0;
   0   0.4 1.4; 
   0.9  0.5   0.5];% the replicator matrix
d=1*[0.1 0.1 0.1];%greater the diffusion, biger must be NT
g=0.1*[1,1,1];

T=150;%must be a whole number
NT=40;

MAS(1:n, 1:(2*m+1))=0; 
Const=1/n/2/pi;
MAS(:,1)=16*Const;%guaranties integral=1

MAS(:,2)=0.3*Const; %initial 
MAS(1,32)=0.3;
MAS(2,3)=0.3;
MAS(3,6)=0.2;

rep_show(MAS,m,n,1);

% [RESMAS f]=rep_4isl(MAS,NT,T,m,n,d,A,g);
% y=1:1:size(f,1);
% figure(3);
% plot(y,f','b-');
% rep_show(RESMAS,m,n,2);

F=[];
h = waitbar(0,'Please wait...');
for cycle_t=1:T
    [MAS f]=rep_4isl(MAS,NT,1,m,n,d,A,g);
    MASFUN=rep_show(MAS,m,n,2);
    DIM_ARRAY(cycle_t,:,:)=MASFUN;
    F=cat(1,F,f);
    waitbar(cycle_t / T);
end;
close(h);
figure(5);
TEMP=permute(DIM_ARRAY(:,1,:),[1 3 2]);
surf(TEMP);
hold on;
TEMP=permute(DIM_ARRAY(:,2,:),[1 3 2]);
surf(TEMP);
TEMP=permute(DIM_ARRAY(:,3,:),[1 3 2]);
surf(TEMP);
colormap hsv;
figure(6);
plot(F);

% figure(4);
% for i=1:n
%     surfc()
    








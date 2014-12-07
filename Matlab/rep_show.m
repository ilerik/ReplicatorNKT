function MASFUN=rep_show(MAS,m,n,fig)

u0=MAS(:,1);% the constants
u=MAS(:,2:(m+1));%coeff for cosinuses
v=MAS(:, (m+2):(2*m+1)); % coeff for sinuses


x=0:0.1:2*pi;
cos_jx=zeros(m,size(x,2));
sin_jx=zeros(m,size(x,2));
const_x=zeros(n,size(x,2));

for j=1:m %forms matrices cos_jx and sin_jx that are needed for the graph. Each column contains the set of cos(jx) for a particular x
    cos_jx(j,:)=cos(j*x);
    sin_jx(j,:)=sin(j*x);
end;
for j=1:n
    const_x(j,:)=u0(j);
end;


MASFUN=u*cos_jx + v*sin_jx +const_x; %forms the matrix of functions. k-th row contains the u(k) for a needed x

figure(fig);
for k=1:n
    plot (x,MASFUN(k,:),'r-');
    hold on;
end;

%sum(u0); were just for proof

function [uu, ff]=rep_4isl(MAS,NT,T,m,n,d,A,g)
%solves the system numericaly, using the Rhunge-Kutta method 
h=T/NT;

[~, fu]=rep_ffunc(MAS,m,n,A,d,g);
f(1)=fu;
for i=2:NT+1;
    
    [ru, fu]=rep_ffunc(MAS,m,n,A,d,g);
    f(i)=fu;

    K0=ru*h;
    
    [ru1, fu1]=rep_ffunc(MAS+0.5*K0,m,n,A,d,g);
    K1=ru1*h;
    
    [ru2, fu2]=rep_ffunc(MAS+0.5*K1,m,n,A,d,g);
    K2=ru2*h;
    
    [ru3, fu3]=rep_ffunc(MAS+K2,m,n,A,d,g);
    K3=ru3*h;
    
    MAS=MAS+(1/6)*(K0+2*K1+2*K2 +K3);

end;

ff=f';
uu=MAS;
%% Particle-based simulation for the heterogeneous receptors at the receiver
fr=0.15;
NN=9;
N=2*NN+1;
theta=pi/2-asin(2*(-NN:NN)/N);
phi=4*pi*(-NN:NN)/(1+sqrt(5));

x=sin(theta).*cos(phi);
y=sin(theta).*sin(phi);
z=cos(theta);

cx=x*10^-5;
cy=y*10^-5;
cz=z*10^-5;
cx=cx';
cy=cy';
cz=cz';

r=10*10^-6;
S=(fr/N)*4*pi*r^2;
h=S/(2*pi*r);
lt=sqrt(2*r^2-2*r*(r-h));
Ntx=1000;
dt=10^-3;
dt2=10^-7;
numT=5*10^7+1;
numT1=5*10^3+1;
MC=10;
kd=0.8;
Pr=1-exp(-kd*dt);
Pr2=1-exp(-kd*dt2);
Xrx=0;
Yrx=0;
Zrx=0;
D=79.4*10^-12;
d_Rx=zeros(MC,numT-1);

dit=20*10^-6;
for p=1:MC
    Xtx=pp(p,1)*23*10^-6;
    Ytx=pp(p,2)*23*10^-6;
    Ztx=pp(p,3)*23*10^-6;
    Xpos=zeros(1,Ntx);
    Xpos(1,:)=Xtx;
    Ypos=zeros(1,Ntx);
    Ypos(1,:)=Ytx;
    Zpos=zeros(1,Ntx);
    Zpos(1,:)=Ztx;
    for i=2:numT1
        Xpos_ini=Xpos;
        Ypos_ini=Ypos;
        Zpos_ini=Zpos;
        d1=find(Xpos.^2+Ypos.^2+Zpos.^2>=dit^2);
        d2=find(Xpos.^2+Ypos.^2+Zpos.^2<dit^2);
        Xpos1=Xpos(d1);Ypos1=Ypos(d1);Zpos1=Zpos(d1);
        Xpos2=Xpos(d2);Ypos2=Ypos(d2);Zpos2=Zpos(d2);
        l1=rand(1,length(Xpos1));
        r1=find(l1<=Pr);
        Xpos1(r1)=[];
        Ypos1(r1)=[];
        Zpos1(r1)=[];
        Xpos1=Xpos1+sqrt(2*D*dt)*randn(1,length(Xpos1));
        Ypos1=Ypos1+sqrt(2*D*dt)*randn(1,length(Xpos1));
        Zpos1=Zpos1+sqrt(2*D*dt)*randn(1,length(Xpos1));
        if isempty(d2)~=1
        for j=1:dt/dt2
        l2=rand(1,length(Xpos2));
        r2=find(l2<=Pr2);
        Xpos2(r2)=[];
        Ypos2(r2)=[];
        Zpos2(r2)=[];
        Xpos2_ini=Xpos2;
        Ypos2_ini=Ypos2;
        Zpos2_ini=Zpos2;
        Xpos2=Xpos2+sqrt(2*D*dt2)*randn(1,length(Xpos2));
        Ypos2=Ypos2+sqrt(2*D*dt2)*randn(1,length(Xpos2));
        Zpos2=Zpos2+sqrt(2*D*dt2)*randn(1,length(Xpos2));
        
        a1=find((Xpos2-Xrx).^2+(Ypos2-Yrx).^2+(Zpos2-Zrx).^2<r^2);
        if isempty(a1)~=1
        xi=Xpos2(a1);
        xa=Xpos2_ini(a1);
        yi=Ypos2(a1);
        ya=Ypos2_ini(a1);
        zi=Zpos2(a1);
        za=Zpos2_ini(a1);
        gamma1=(xa-xi).^2+(ya-yi).^2+(za-zi).^2;
        gamma2=2*(ya-yi).*(xa.*yi-xi.*ya)+2.*(za-zi).*(xa.*zi-xi.*za);
        gamma3=xi.*(ya-yi).*(xi.*ya-2.*yi.*xa+yi.*xi)+xi.*(za-zi).*(xi.*za-2*xa.*zi+xi.*zi)+(xa-xi).^2.*(yi.^2+zi.^2-r^2);
        x1=(-gamma2+sqrt(gamma2.^2-4*gamma1.*gamma3))./(2*gamma1);
        x2=(-gamma2-sqrt(gamma2.^2-4*gamma1.*gamma3))./(2*gamma1);
        index1=(xa-x1).*(xi-x1);
        index2=(xa-x2).*(xi-x2);
        p1=index1<0;
        p2=index2<0;
        x=p1.*x1+p2.*x2;
        y=(x-xi).*(ya-yi)./(xa-xi)+yi;
        z=(x-xi).*(za-zi)./(xa-xi)+zi;
        
        a2=(x-cx).^2+(y-cy).^2+(z-cz).^2<=lt^2;
        af=sum(a2);
       
        c1=find(af==1);
        c2=find(af==0);
        ll=floor((i-2)*dt/dt2+j);
        d_Rx(p,ll)=sum(af);
        
        Xpos2(a1(c2))=Xpos2_ini(a1(c2));
        Ypos2(a1(c2))=Ypos2_ini(a1(c2));
        Zpos2(a1(c2))=Zpos2_ini(a1(c2));
        
        index=a1(c1);
        
        Xpos2(index)=[];
        Ypos2(index)=[];
        Zpos2(index)=[];
        Xpos2_ini(index)=[];
        Ypos2_ini(index)=[];
        Zpos2_ini(index)=[];  
        end
        end
        end
        Xpos=[Xpos1,Xpos2];
        Ypos=[Ypos1,Ypos2];
        Zpos=[Zpos1,Zpos2];
    end
end
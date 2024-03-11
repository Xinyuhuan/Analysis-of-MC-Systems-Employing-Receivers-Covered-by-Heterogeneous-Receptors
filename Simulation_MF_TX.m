%% Particle-based simulation when the transmitter is the MF-based transmitter
fr=0.05;
NN=2;
N=2*NN+1;
theta=pi/2-asin(2*(-NN:NN)/N);
phi=4*pi*(-NN:NN)/(1+sqrt(5));

cx=sin(theta).*cos(phi)*10^-5;
cy=sin(theta).*sin(phi)*10^-5;
cz=cos(theta)*10^-5;

cx=cx';
cy=cy';
cz=cz';
r_rx=10*10^-6;
S=(fr/N)*4*pi*r_rx^2;
hh=S/(2*pi*r_rx);
lt=sqrt(2*r_rx^2-2*r_rx*(r_rx-hh));
Ntx=200;
eta=5;
dt=10^-3;
dt2=10^-5;
numT=2*10^3+1;
numT2=2*10^3+1;
numT3=2*10^5+1;
rtx=5*10^-6;
MC=100;
D=9*10^-12;
Do=79.4*10^-12;
d_Rx=zeros(MC,numT2-1);
kf=30*10^-6;
Pr=kf*sqrt(pi*dt/D);
h_kd=0.8;
Prd2=1-exp(-h_kd.*dt);
Prd3=1-exp(-h_kd.*dt2);
Xrx=0;
Yrx=0;
Zrx=0;
timeRange=(0:(numT2-1))*dt;
timeRange2=(0:(numT3-1))*dt2;

dit=20*10^-6;
for p=1:MC
    Xpos=zeros(1,Ntx);
    Ypos=zeros(1,Ntx);
    Zpos=zeros(1,Ntx);
    Xemit=[];
    Yemit=[];
    Zemit=[];
    temit=[];
    for i=2:numT
        Xpos_ini=Xpos;
        Ypos_ini=Ypos;
        Zpos_ini=Zpos;
        Xpos=Xpos+sqrt(2*D*dt)*randn(1,length(Xpos_ini));
        Ypos=Ypos+sqrt(2*D*dt)*randn(1,length(Ypos_ini));
        Zpos=Zpos+sqrt(2*D*dt)*randn(1,length(Zpos_ini));
        
        a=sqrt(Xpos.^2+Ypos.^2+Zpos.^2);
        e1=find(a>rtx);
        l1=rand(1,length(e1));
        r1=find(l1>Pr);
        Xpos(e1(r1))=Xpos_ini(e1(r1));
        Ypos(e1(r1))=Ypos_ini(e1(r1));
        Zpos(e1(r1))=Zpos_ini(e1(r1));
        r11=find(l1<=Pr);
        xi1=Xpos_ini(e1(r11));
        xa1=Xpos(e1(r11));
        yi1=Ypos_ini(e1(r11));
        ya1=Ypos(e1(r11));
        zi1=Zpos_ini(e1(r11));
        za1=Zpos(e1(r11));
        gamma11=(xa1-xi1).^2+(ya1-yi1).^2+(za1-zi1).^2;
        gamma21=2*(ya1-yi1).*(xa1.*yi1-xi1.*ya1)+2.*(za1-zi1).*(xa1.*zi1-xi1.*za1);
        gamma31=xi1.*(ya1-yi1).*(xi1.*ya1-2.*yi1.*xa1+yi1.*xi1)+xi1.*(za1-zi1).*(xi1.*za1-2*xa1.*zi1+xi1.*zi1)+(xa1-xi1).^2.*(yi1.^2+zi1.^2-rtx^2);
        x11=(-gamma21+sqrt(gamma21.^2-4*gamma11.*gamma31))./(2*gamma11);
        x21=(-gamma21-sqrt(gamma21.^2-4*gamma11.*gamma31))./(2*gamma11);
        index1=(xa1-x11).*(xi1-x11);
        index21=(xa1-x21).*(xi1-x21);
        p11=index1<0;
        p21=index21<0;
        x1=p11.*x11+p21.*x21;
        y1=(x1-xi1).*(ya1-yi1)./(xa1-xi1)+yi1;
        z1=(x1-xi1).*(za1-zi1)./(xa1-xi1)+zi1;
        xv=repmat(x1,1,eta);
        yv=repmat(y1,1,eta);
        zv=repmat(z1,1,eta);
        Xemit=[Xemit,xv];
        Yemit=[Yemit,yv];
        Zemit=[Zemit,zv];
        
        dli=(xa1-xi1).^2+(ya1-yi1).^2+(za1-zi1).^2;
        dla=(x1-xi1).^2+(y1-yi1).^2+(z1-zi1).^2;
        ddt=dt.*dla./dli;
        ddtf=ddt+(i-2)/1000;
        tv=repmat(ddtf,1,eta);
        temit=[temit,tv];
        

        Xpos(e1(r11))=[];
        Ypos(e1(r11))=[];
        Zpos(e1(r11))=[];
        Xpos_ini(e1(r11))=[];
        Ypos_ini(e1(r11))=[];
        Zpos_ini(e1(r11))=[];
    end
    Xs=Xemit+pp(p,1)*20*10^-6;
    Ys=Yemit+pp(p,2)*20*10^-6;
    Zs=Zemit+pp(p,3)*20*10^-6;
    ts=temit;
    for j=2:numT2
        d1=find(Xs.^2+Ys.^2+Zs.^2>=dit^2);
        d2=find(Xs.^2+Ys.^2+Zs.^2<dit^2);
        Xs1=Xs(d1);Ys1=Ys(d1);Zs1=Zs(d1);ts1=ts(d1);
        Xs2=Xs(d2);Ys2=Ys(d2);Zs2=Zs(d2);ts2=ts(d2);
        f1=find(ts1<=timeRange(j-1));
        g1=rand(1,length(f1));
        w1=find(g1>Prd2);
        Xs1(f1(w1))=Xs1(f1(w1))+sqrt(2*Do*dt)*randn(1,length(f1(w1)));
        Ys1(f1(w1))=Ys1(f1(w1))+sqrt(2*Do*dt)*randn(1,length(f1(w1)));
        Zs1(f1(w1))=Zs1(f1(w1))+sqrt(2*Do*dt)*randn(1,length(f1(w1)));
        w11=find(g1<=Prd2);
        Xs1(f1(w11))=[];
        Ys1(f1(w11))=[];
        Zs1(f1(w11))=[];
        ts1(f1(w11))=[];
        
        f2=find(ts1<=timeRange(j)&ts1>timeRange(j-1));
        g2=rand(1,length(f2));
        w2=find(g2>Prd2);
        Xs1(f2(w2))=Xs1(f2(w2))+sqrt(2*Do*(timeRange(j)-ts(f2(w2)))).*randn(1,length(f2(w2)));
        Ys1(f2(w2))=Ys1(f2(w2))+sqrt(2*Do*(timeRange(j)-ts(f2(w2)))).*randn(1,length(f2(w2)));
        Zs1(f2(w2))=Zs1(f2(w2))+sqrt(2*Do*(timeRange(j)-ts(f2(w2)))).*randn(1,length(f2(w2)));
        w22=find(g2<=Prd2);
        Xs1(f2(w22))=[];
        Ys1(f2(w22))=[];
        Zs1(f2(w22))=[];
        ts1(f2(w22))=[];
      
        f3=find(timeRange(j)<ts1);
        Xs1(f3)=Xs1(f3);
        Ys1(f3)=Ys1(f3);
        Zs1(f3)=Zs1(f3);
        if isempty(d2)~=1
            for k=1:dt/dt2
                Xs2_ini=Xs2;
                Ys2_ini=Ys2;
                Zs2_ini=Zs2;
                in=floor((j-2)*dt/dt2+k);
                f4=find(ts2<=timeRange2(in));
                g4=rand(1,length(f4));
                w4=find(g4>Prd3);
                Xs2(f4(w4))=Xs2(f4(w4))+sqrt(2*Do*dt2)*randn(1,length(f4(w4)));
                Ys2(f4(w4))=Ys2(f4(w4))+sqrt(2*Do*dt2)*randn(1,length(f4(w4)));
                Zs2(f4(w4))=Zs2(f4(w4))+sqrt(2*Do*dt2)*randn(1,length(f4(w4)));
                w44=find(g4<=Prd3);
                Xs2(f4(w44))=[];
                Ys2(f4(w44))=[];
                Zs2(f4(w44))=[];
                Xs2_ini(f4(w44))=[];
                Ys2_ini(f4(w44))=[];
                Zs2_ini(f4(w44))=[];
                ts2(f4(w44))=[];
        
                f5=find(ts2<=timeRange2(in+1)&ts2>timeRange2(in));
                g5=rand(1,length(f5));
                w5=find(g5>Prd3);
                Xs2(f5(w5))=Xs2(f5(w5))+sqrt(2*Do*(timeRange2(in+1)-ts(f5(w5)))).*randn(1,length(f5(w5)));
                Ys2(f5(w5))=Ys2(f5(w5))+sqrt(2*Do*(timeRange2(in+1)-ts(f5(w5)))).*randn(1,length(f5(w5)));
                Zs2(f5(w5))=Zs2(f5(w5))+sqrt(2*Do*(timeRange2(in+1)-ts(f5(w5)))).*randn(1,length(f5(w5)));
                w55=find(g5<=Prd3);
                Xs2(f5(w55))=[];
                Ys2(f5(w55))=[];
                Zs2(f5(w55))=[];
                Xs2_ini(f5(w55))=[];
                Ys2_ini(f5(w55))=[];
                Zs2_ini(f5(w55))=[];
                ts2(f5(w55))=[];
      
                f6=find(timeRange2(in+1)<ts2);
                Xs2(f6)=Xs2(f6);
                Ys2(f6)=Ys2(f6);
                Zs2(f6)=Zs2(f6);
                
                a1=find((Xs2-Xrx).^2+(Ys2-Yrx).^2+(Zs2-Zrx).^2<r_rx^2);
                if isempty(a1)~=1
                   xi=Xs2(a1);
                   xa=Xs2_ini(a1);
                   yi=Ys2(a1);
                   ya=Ys2_ini(a1);
                   zi=Zs2(a1);
                   za=Zs2_ini(a1);
                   gamma1=(xa-xi).^2+(ya-yi).^2+(za-zi).^2;
                   gamma2=2*(ya-yi).*(xa.*yi-xi.*ya)+2.*(za-zi).*(xa.*zi-xi.*za);
                   gamma3=xi.*(ya-yi).*(xi.*ya-2.*yi.*xa+yi.*xi)+xi.*(za-zi).*(xi.*za-2*xa.*zi+xi.*zi)+(xa-xi).^2.*(yi.^2+zi.^2-r_rx^2);
                   x1=(-gamma2+sqrt(gamma2.^2-4*gamma1.*gamma3))./(2*gamma1);
                   x2=(-gamma2-sqrt(gamma2.^2-4*gamma1.*gamma3))./(2*gamma1);
                   index=(xa-x1).*(xi-x1);
                   index2=(xa-x2).*(xi-x2);
                   p1=index<0;
                   p2=index2<0;
                   x=p1.*x1+p2.*x2;
                   y=(x-xi).*(ya-yi)./(xa-xi)+yi;
                   z=(x-xi).*(za-zi)./(xa-xi)+zi;
        
                   a2=(x-cx).^2+(y-cy).^2+(z-cz).^2<=lt^2;
                   af=sum(a2);
     
                   c1=find(af==1);
                   c2=find(af==0);
                   d_Rx(p,in)=sum(af);
        
                   Xs2(a1(c2))=Xs2_ini(a1(c2));
                   Ys2(a1(c2))=Ys2_ini(a1(c2));
                   Zs2(a1(c2))=Zs2_ini(a1(c2));
        
                   Xs2(a1(c1))=[];
                   Ys2(a1(c1))=[];
                   Zs2(a1(c1))=[];
                   Xs2_ini(a1(c1))=[];
                   Ys2_ini(a1(c1))=[];
                   Zs2_ini(a1(c1))=[];
                   ts2(a1(c1))=[];
               end
            end
       end
       Xs=[Xs1,Xs2];
       Ys=[Ys1,Ys2];
       Zs=[Zs1,Zs2];
       ts=[ts1,ts2];
    end
end

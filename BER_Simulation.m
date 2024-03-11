%% Simulation for ploting the bit error rate when the transmitter is MF-based
%% TX and receiver is AP-based RX
delta=0.001;
t=0.001:delta:10;
rtx=5*10^-6;
Da=27*10^-12;
kf=30*10^-6;
x=0:1:50*10^6;
y=Da./rtx.*cos(x.*rtx)-Da./(x.*rtx^2).*sin(x.*rtx)+kf.*sin(x.*rtx)./(x.*rtx);
lambdan=zeros(1,length(x));
for i=2:length(x)
    a1=y(i);
    a2=y(i-1);
    if a1*a2<0
        lambdan(i)=x(i-1);
    else
        lambdan(i)=0;
    end
end
index=find(lambdan==0);
lambdan(index)=[];

frt=zeros(1,length(t));
for i=1:length(t)
    a=4*rtx^2*kf*lambdan.^3./(2*lambdan.*rtx-sin(2*lambdan.*rtx)).*sin(lambdan.*rtx)./(lambdan.*rtx).*exp(-Da.*lambdan.^2.*t(i));
    frt(i)=sum(a);
end



NN=100;
N=2*NN+1;
D=79.4*10^-12;
r0=20*10^-6;
rr=10*10^-6;
k=0.8;
fr=0.1;
beta=(r0-rr)*sqrt(k./D);

theta=pi/2-asin(2*(-NN:NN)/N);
phi=4*pi*(-NN:NN)/(1+sqrt(5));
x=sin(theta).*cos(phi);
y=sin(theta).*sin(phi);
z=cos(theta);
u=zeros(1,N*(N-1)/2);
o=1;
for i=1:N
    for j=i+1:N
        a1=[x(i), y(i), z(i)];
        a2=[x(j), y(j), z(j)];
        u(o)=norm(a1-a2);
        o=o+1;
    end
end

H=1./u+1/2*log(u)-1/2*log(2+u);
H=sum(H);

sig=2*sqrt(fr/N);
C0i=pi/(N*sig)*(1+sig/pi*log(sig/2)+sig/pi*(log(4)-3/2+4/N*H));
C0=1/C0i*rr;

we=D*C0/(rr*(rr-C0));
gamma=(we*rr+D)/(D*rr);
zeta=gamma^2*D-k;
omega=gamma*sqrt(D);
z1=r0-rtx-rr;
z2=r0+rtx-rr;
xi21=exp(gamma*z1+zeta*t).*erfc(z1./sqrt(4*D*t)+omega*sqrt(t))-1/(2*sqrt(k))*exp(-z1*sqrt(k/D)).*((omega-sqrt(k)).*erf(z1./sqrt(4*D*t)-sqrt(k*t))-(omega+sqrt(k)).*(1-exp(2*z1*sqrt(k/D)).*erfc(z1./sqrt(4*D*t)+sqrt(k*t))))-exp(-z1*sqrt(k/D));
xi22=exp(gamma*z2+zeta*t).*erfc(z2./sqrt(4*D*t)+omega*sqrt(t))-1/(2*sqrt(k))*exp(-z2*sqrt(k/D)).*((omega-sqrt(k)).*erf(z2./sqrt(4*D*t)-sqrt(k*t))-(omega+sqrt(k)).*(1-exp(2*z2*sqrt(k/D)).*erfc(z2./sqrt(4*D*t)+sqrt(k*t))))-exp(-z2*sqrt(k/D));
Hst=rr*we/(2*rtx*r0*zeta)*(xi21-xi22);

Hmf=conv(frt,Hst)*delta;
tf=0.002:delta:20;
gh=[0.8,1.6,2.4,3.2,4,4.8,5.6,6.4,7.2,8];
inde=gh/0.001-1;
inde=floor(inde);
F=Hmf(inde)*1000;

b=zeros(1,10);
for i=1:9
    b(i)=F(11-i)-F(10-i);
end
b(end)=F(1);

MC=30000;
bf=zeros(length(b),MC);
for i=1:length(b)
    bf(i,:)=random('Poisson',b(i),1,MC);
end
eta=[28,43,54,62,72,80,94,112];
n=[10,9,8,7,6,5,4,3,2,1];
prf=zeros(length(n),length(eta));

for i=1:length(n)
    a=abs(dec2bin(0:(2^n(i)-1), n(i)))-48;
    [col,~]=size(a);
    bbf=zeros(col,length(eta));
    for j=1:col
        chi=a(j,:)*bf(i:10,:);
        for m=1:length(eta)
            if a(j,end)==1
                ber=chi<eta(m);
            else
                ber=chi>=eta(m);
            end
            bbf(j,m)=sum(ber)/MC;
        end
    end
    prf(i,:)=mean(bbf);
end
semilogy(eta,mean(prf),'ko')

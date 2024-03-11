%% Molecule Hitting Rate at AP-based receiver
delta=0.0001;
t=0.001:delta:3;
rtx=5*10^-6;
Da=27*10^-12;
kf=30*10^-6;
x=0:1:32*10^6;
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
% x=pp(1:9,1);
% y=pp(1:9,2);
% z=pp(1:9,3);
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
xi11=exp(gamma*z1+zeta*t).*erfc(omega*sqrt(t)+z1./sqrt(4*D*t));
xi12=exp(gamma*z2+zeta*t).*erfc(omega*sqrt(t)+z2./sqrt(4*D*t));
hst=rr*we./(2*rtx*r0)*(xi11-xi12);

hmf=conv(frt,hst)*delta;
tf=0.002:delta:6;
plot(tf,hmf);
xlim([0 1.5])
hold on

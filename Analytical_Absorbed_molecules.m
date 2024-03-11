%% This code plot the asymptotic number of absorbed molecules
NN=56;
N=2*NN+1;
D=79.4*10^-12;
r0=23*10^-6;
rr=10*10^-6;
k=0.8;
fr=0.15;
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

% C0i=pi/sig*(1+sig/pi*(log(2*sig)-3/2)-sig^2/pi^2*((pi^2+21)/36));
% C0=1/C0i*rr;
w=D*C0/(rr*(rr-C0));

gamma=(w*rr+D)/(D*rr);
zeta=gamma^2*D-k;
F=rr*w*(gamma-sqrt(k/D))/(r0*zeta)*exp(-beta);

t=0:0.01:5;
Fs=zeros(1,length(t));
Fs(1,:)=F*1000;
plot(t,Fs);
hold on

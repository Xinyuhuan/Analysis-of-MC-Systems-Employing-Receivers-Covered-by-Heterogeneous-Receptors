NN=1:5:101;
N=2*NN+1;
w=zeros(1,length(NN));
for k=1:length(NN)
theta=pi/2-asin(2*(-NN(k):NN(k))/N(k));
phi=4*pi*(-NN(k):NN(k))/(1+sqrt(5));
x=sin(theta).*cos(phi);
y=sin(theta).*sin(phi);
z=cos(theta);
u=zeros(1,N(k)*(N(k)-1)/2);
o=1;
for i=1:N(k)
    for j=i+1:N(k)
        a1=[x(i), y(i), z(i)];
        a2=[x(j), y(j), z(j)];
        u(o)=norm(a1-a2);
        o=o+1;
    end
end

H=1./u+1/2*log(u)-1/2*log(2+u);
H=sum(H);
fr=0.05;
D=100*10^-12;
r0=20*10^-6;
rr=10*10^-6;
sig=2*sqrt(fr/N(k));
C0i=pi/(N(k)*sig)*(1+sig/pi*log(sig/2)+sig/pi*(log(4)-3/2+4/N(k)*H));
C0=1/C0i;
w(k)=D*C0/(rr*(1-C0));
end
sig1=2*sqrt(fr);
C0i_2=pi/sig1*(1+sig1/pi*(log(sig1)-3/2+log(2))-sig1^2/pi^2*((pi^2+21)/36));
C02=1/C0i_2;
w1=D*C02/(rr*(1-C02));
kd=0.8;

w=[w1,w];
gamma=(w*rr+D)/(D*rr);
zeta=gamma.^2*D-kd;

y=w*zeta(1).*(gamma-sqrt(kd/D))./(w(1)*zeta.*(gamma(1)-sqrt(kd./D)))-1;
N=[1,N];
plot(N,y);
hold on
% N=[1,N];
% plot(N,w*10^6);
% hold on

r = 550e-9;
R1=linspace(1e-23,r,nth);

if numel(theta)==1
 [R,phi]=meshgrid(R1,phi1);
end
if numel(phi)==1
 [R,theta]=meshgrid(R1,theta1); 
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
% initial values for sums 

Mu= cos(theta);
Ms= sin(theta);
pin=zeros(numel(Mu),N);
taon=zeros(numel(Mu),N);

E0=1;

Elp=0;
Elt=0;
Elr=0;
%first and second terms of recurrence relationship of tao and pi
pin1=0;  pin2=1;
taon1=0; taon2=Mu;
pin=1;
taon=Mu;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

d1=0;
c1=0;
d=zeros(1,N);
c=zeros(1,N);

%notes Y and X
X=R.*sin(theta).*cos(phi);
Y=R.*sin(theta).*sin(phi);
Z=R.*cos(theta);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

 for u=1:N%a and b loop
    n=u;
        n12=u+1/2;

     
    psix=x.*sqrt(pi./(2.*x)).*besselj(n12,x);
%     psix_1=x.*sqrt(pi./(2.*x)).*besselj(n12-1,x);
    psimx=m.*x.*sqrt(pi./(2.*m.*x)).*besselj(n12,m.*x);
%     psimx_1=m.*x.*sqrt(pi./(2.*m.*x)).*besselj(n12-1,m.*x);
    xix=x.*sqrt(pi./(2.*x)).*besselh(n12,x);
%    xix_1=x.*sqrt(pi./(2.*x)).*besselh(n12-1,x);
    ximx=m.*x.*sqrt(pi./(2.*m.*x)).*besselh(n12,m.*x);

    xiprimemx1=1/2*sqrt(pi/2)*sqrt(1./m.*x).*(m.*x.*besselh(n12-1,m.*x)+besselh(n12,m.*x)-m.*x.*besselh(n12+1,m.*x));
    xiprimex1=1/2*sqrt(pi/2)*sqrt(1./x).*(x.*besselh(n12-1,x)+besselh(n12,x)-x.*besselh(n12+1,x));
    psiprimex1=1/2*sqrt(pi/2)*sqrt(1./x).*(x.*besselj(n12-1,x)+besselj(n12,x)-x.*besselj(n12+1,x));
    psiprimemx1=1/2*sqrt(pi/2)*sqrt(1./mx).*(mx.*besselj(n12-1,mx)+besselj(n12,mx)-mx.*besselj(n12+1,mx));

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
    dut=(m.*psix.*xiprimex1-m.*xix.*psiprimex1);
    dub=(m.*psimx.*xiprimex1-psiprimemx1.*xix);
    d1=dut./dub;
%      d coieffent
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
     cut=(psix.*xiprimex1-xix.*psiprimex1);
     cub=(psimx.*xiprimex1-m.*psiprimemx1.*xix);

     c1=cut./cub;
%      c coeifient
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
   
      d(n)=d1;
      c(n)=c1;
 end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  


for nb=1:N
    nc=nb+1;
    n12=nb+1/2;

    k=2*pi*Ns./l;
rho=k.*R;

if nb>2
  %using the 2 previous polynormial terms
    pin=(2*nc-1)/(nc-1)*Mu'.*pin2-(nc)/(nc-1)*pin1; 
    
    taon=(nc)*pin.*Mu'-(nc+1)*pin2;
pin1=pin2;
pin2=pin;
end
   Jrho=sqrt(pi./(2.*rho)).*besselj(n12,rho);                    
   primeJrho=1/2*sqrt(pi/2)*sqrt(1./rho).*(rho.*besselj(n12-1,rho)+besselj(n12,rho)-rho.*besselj(n12+1,rho));
   
   Melnt1=-pin.*Jrho.*sin(phi);
    Melnp1=-taon.*Jrho.*cos(phi);

Nolnt1=taon.*(primeJrho)./(rho).*sin(phi);
Nolnp1=pin.*(primeJrho)./(rho).*cos(phi);
Nolnr1=pin.*(Jrho./(rho)).*Ms.*sin(phi)*(nb*(nb+1));
% end
En=1i^nb*(2*nb+1)/(nb*(nb+1))*E0;

Elr=En*(c(nb)*0+1i*d(nb).*Nolnr1)+Elr;
Elt=En*(c(nb).*Melnt1+1i*d(nb).*Nolnt1)+Elt;
Elp=En*(c(nb).*Melnp1+1i*d(nb).*Nolnp1)+Elp;
end
%Es below equals Et
Elr=Eir+Elr;
Elt=Eit+Elt;
Elp=Eip+Elp;
El_mods=(Elr.*conj(Elr)+Elt.*conj(Elt)+Elp.*conj(Elp)).^(1/2);
Ei_mod=sqrt(abs(Eir).^2+abs(Eit).^2+abs(Eip).^2);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


figure(1)
hold on
surface(Y*1e9,Z*1e9,El_mods)
shading flat
colormap hot
ylabel('Z (nm)')
xlabel('Y (nm)')
colorbar
axis equal


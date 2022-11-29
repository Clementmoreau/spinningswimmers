function out = fig_flow_draw(label)


global Yinit Yinit2
global Ntheta Nphi 
global alpha beta


NQ=20;

filenum=num2str(label,'%02d');
filename=strcat('fig_flow',filenum,'.eps');

delta_theta=pi/Ntheta;
delta_phi=2*pi/Nphi;

[phi,theta] = meshgrid(-2*pi:delta_phi:2*pi,0:delta_theta:pi);

U=beta*cos(phi).*cos(phi) + 0.5*(1-beta)+alpha*cos(theta).*cos(phi).*sin(phi);
V= 0.25*beta*sin(2*theta).*sin(2*phi)-0.5*alpha*sin(theta).*cos(2*phi);

[phiQ,thetaQ] = meshgrid(-2*pi:delta_phi*NQ:2*pi,0:delta_theta*NQ:pi);

UQ=beta*cos(phiQ).*cos(phiQ) + 0.5*(1-beta)+alpha*cos(thetaQ).*cos(phiQ).*sin(phiQ);
VQ= 0.25*beta*sin(2*thetaQ).*sin(2*phiQ)-0.5*alpha*sin(thetaQ).*cos(2*phiQ);

starty = 0.1:0.1:0.9;

startx1 = (Yinit/pi)*ones(size(starty));
starty1= starty;

startx2 = (Yinit/pi)*ones(size(starty));
starty2= starty;

startx3 = -(Yinit2/pi)*ones(size(starty));
starty3= starty;

startx4 = -(Yinit2/pi)*ones(size(starty));
starty4= starty;

% -- start to figure generation --
dpi = '-r400';
figure;

hold on
q=quiver(phiQ/pi,thetaQ/pi,UQ/pi,VQ/pi,'LineWidth',1);

str1=streamline(phi/pi,theta/pi, U/pi, V/pi,startx1,starty1);
str2=streamline(phi/pi,theta/pi, U/pi, V/pi,startx2,starty2);
str3=streamline(phi/pi,theta/pi,-U/pi,-V/pi,startx3,starty3);
str4=streamline(phi/pi,theta/pi,-U/pi,-V/pi,startx4,starty4);

set(str1,'LineWidth',2,'Color','r') 
set(str2,'LineWidth',2,'Color','r') 
set(str3,'LineWidth',2,'Color','r') 
set(str4,'LineWidth',2,'Color','r') 

if label==6
Ntheta=800;
Nphi=1600;
NQ=40;
delta_theta=pi/Ntheta;
delta_phi=2*pi/Nphi;
[phi,theta] = meshgrid(-2*pi:delta_phi:2*pi,0:delta_theta:pi);
U=beta*cos(phi).*cos(phi) + 0.5*(1-beta)+alpha*cos(theta).*cos(phi).*sin(phi);
V= 0.25*beta*sin(2*theta).*sin(2*phi)-0.5*alpha*sin(theta).*cos(2*phi);
[phiQ,thetaQ] = meshgrid(-2*pi:delta_phi*NQ:2*pi,0:delta_theta*NQ:pi);
UQ=beta*cos(phiQ).*cos(phiQ) + 0.5*(1-beta)+alpha*cos(thetaQ).*cos(phiQ).*sin(phiQ);
VQ= 0.25*beta*sin(2*thetaQ).*sin(2*phiQ)-0.5*alpha*sin(thetaQ).*cos(2*phiQ);    
startx5 = [-0.75,-0.25, 0.25, 0.75,-0.75,-0.25, 0.25, 0.75];
starty5=  [ 0.8, 0.2, 0.8, 0.2, 0.9, 0.1, 0.9, 0.1];
str5=streamline(phi/pi,theta/pi,U/pi,V/pi,startx5,starty5);
set(str5,'LineWidth',2,'Color','r') 
end

if label==7
starty = 0.1:0.1:0.9;
startx5 = (-0.5)*ones(size(starty));
starty5= starty;
startx6 = (0.5)*ones(size(starty));
starty6= starty;
startx7 = -(-0.5)*ones(size(starty));
starty7= starty;
startx8 = -(0.5)*ones(size(starty));
starty8= starty;
startx9 = (0.0)*ones(size(starty));
starty9= starty;
startx10 = (-1.0)*ones(size(starty));
starty10= starty;
startx11 = -(0.0)*ones(size(starty));
starty11= starty;
startx12 = -(-1.0)*ones(size(starty));
starty12= starty;
str5=streamline(phi/pi,theta/pi, U/pi, V/pi,startx5,starty5);
str6=streamline(phi/pi,theta/pi, U/pi, V/pi,startx6,starty6);
str7=streamline(phi/pi,theta/pi,-U/pi,-V/pi,startx7,starty7);
str8=streamline(phi/pi,theta/pi,-U/pi,-V/pi,startx8,starty8);
str9=streamline(phi/pi,theta/pi, U/pi, V/pi,startx9,starty9);
str10=streamline(phi/pi,theta/pi, U/pi, V/pi,startx10,starty10);
str11=streamline(phi/pi,theta/pi,-U/pi,-V/pi,startx11,starty11);
str12=streamline(phi/pi,theta/pi,-U/pi,-V/pi,startx12,starty12);
set(str5,'LineWidth',2,'Color','r') 
set(str6,'LineWidth',2,'Color','r') 
set(str7,'LineWidth',2,'Color','r') 
set(str8,'LineWidth',2,'Color','r')   
set(str9,'LineWidth',2,'Color','r') 
set(str10,'LineWidth',2,'Color','r') 
set(str11,'LineWidth',2,'Color','r') 
set(str12,'LineWidth',2,'Color','r')   
end

xlabel('$\phi/\pi$','Interpreter','latex')
ylabel('$\theta/\pi$','Interpreter','latex')
xlim([-1 1])
xticks([-1:0.5:1])
ylim([0 1])
yticks([0:0.25:1])

set(gca,'FontSize',24)
%set(gca,'FontSize',12)
set(gca,'TickLabelInterpreter','latex')


hold off
grid on
box on

print(filename,'-depsc',dpi) 

end
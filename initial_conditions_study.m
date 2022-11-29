wperp=1;
wpar=1;

Nth = 25;
Nph = 25;

psi = -pi/2;

figure(7); clf
set(gcf, 'Position',  [1, 250, 1000, 1000])

% Plot the unit sphere
[sx,sy,sz]=sphere(30);
s = surf(sx,sy,sz,'FaceLighting','gouraud','FaceAlpha',0.35,'FaceColor',[1 1 1],'EdgeColor',[0.8 0.8 0.8],'EdgeAlpha',0.4);
% s.LineStyle ='none';
s.LineWidth = 1;

material dull
hold on

% Plot the axis and (xy) unit circle
plot3([0 1.5],[0 0],[0 0],'k','LineWidth',1)
plot3([0 0],[0 1.5],[0 0],'k','LineWidth',1)
plot3([0 0],[0 0],[0 1.5],'k','LineWidth',1)
tt=0:0.1:2*pi+0.1;
plot3(cos(tt),sin(tt),zeros(1,length(tt)),'k','LineWidth',1)

view(120,25)
camlight

tic

for i = 1:Nth
    for j = 1:Nph

        theta = i*pi/Nth;
        phi = 2*j*pi/Nph;

        [alpha,mu,varphi] = Initial_conditions_2(theta,phi,psi,wperp,wpar);

        plot3(sin(theta).*cos(phi),sin(theta).*sin(phi),cos(theta),'ob','LineWidth',1)
        plot3(sin(alpha).*cos(varphi),sin(alpha).*sin(varphi),cos(alpha),'xr','LineWidth',3)

        axis equal off

    end
end

toc




%%%%%%%%%% Figure for wobbling swimmer with translational dynamics %%%%%%%%%
% (without chirality)

% Figure for the rotational dynamics in the general dynamics.

%% Setup.
F1 = figure(1);
delete(findall(gcf,'type','annotation'))
tl = tiledlayout(5,3);
options = odeset('RelTol',1e-9,'AbsTol',1e-9);

% Global parameters.
% G = gamma, the shear rate of the flow.
G = 1;

% Bretherton constants.
B_slender = 0.99;
B_spheroidal = 0.5;
B_oblate = -0.99;

% Initial conditions for each case.
Init_slender = [pi/18,pi/2,-pi/4];
Init_spheroidal = [pi/18,pi/2,-pi/4];
Init_oblate = [pi/18,0,-pi/4];

% Final times.
T_init = 0; 
T_max_slender = [100,60,25,20,15];
T_max_non_slender = [20,20,15,15,15];

%% 1. No spinning (classical Jeffery orbits)
tic
W_par = 0;
W_perp = 0;

% a. Slender object
B = B_slender;
Init = Init_slender;
T_max = T_max_slender(1);%time interval

[tps,y] = ode15s(@(t,y) ode_num(t,y,G,B,W_perp,W_par),[T_init T_max],Init,options);% Solve full ODE
theta_nospin_slender = y(:,1);phi_nospin_slender = y(:,2);psi_nospin_slender = y(:,3);

nexttile;
plot_sphere_nospin(theta_nospin_slender,phi_nospin_slender)

% b. Spheroidal object
B = B_spheroidal;
Init = Init_spheroidal;
T_max = T_max_non_slender(1);%time interval

[tps,y] = ode15s(@(t,y) ode_num(t,y,G,B,W_perp,W_par),[T_init T_max],Init,options);% Solve full ODE
theta_nospin_spheroidal = y(:,1);phi_nospin_spheroidal = y(:,2);psi_nospin_spheroidal = y(:,3);

nexttile;
plot_sphere_nospin(theta_nospin_spheroidal,phi_nospin_spheroidal)

% c. Oblate object
B = B_oblate;
Init = Init_oblate;
T_max = T_max_slender(1);%time interval

[tps,y] = ode15s(@(t,y) ode_num(t,y,G,B,W_perp,W_par),[T_init T_max],Init,options);% Solve full ODE
theta_nospin_oblate = y(:,1);phi_nospin_oblate = y(:,2);psi_nospin_oblate = y(:,3);

nexttile;
plot_sphere_nospin(theta_nospin_oblate,phi_nospin_oblate)

toc
%% 2. Bacterial sublimit
tic

W_par = 15;
W_perp = 1.5;

% Quantities used in the analysis.
w = W_par/W_perp; %spinning ratio
lambda = sqrt(1 + w.^2); 

% a. Slender object
B = B_slender;
Init = Init_slender;
T_max = T_max_slender(2);%time interval
Init_theta = Init(1);Init_phi = Init(2); Init_psi = Init(3);
[ alp_init, mu_init, phibar_init ] = Initial_conditions( Init_theta, Init_phi, Init_psi, W_perp, W_par );
Init_asy = [alp_init, phibar_init, mu_init];

% Solve full ODE
[tps,y] = ode15s(@(t,y) ode_num(t,y,G,B,W_perp,W_par),[T_init T_max],Init,options);
theta = y(:,1);phi = y(:,2);psi = y(:,3);
% Solve the reduced ODE
[ t_asy,y_asy ] = Solve_reduced_ODE(G,B,w,Init_asy,tps);
alp = y_asy(:,1);phibar = y_asy(:,2);mu = y_asy(:,3);

nexttile;
plot_sphere(theta,phi,alp,phibar);

% b. Spheroidal object
B = B_spheroidal;
Init = Init_spheroidal;
T_max = T_max_non_slender(2);%time interval
Init_theta = Init(1);Init_phi = Init(2); Init_psi = Init(3);
[ alp_init, mu_init, phibar_init ] = Initial_conditions( Init_theta, Init_phi, Init_psi, W_perp, W_par );
Init_asy = [alp_init, phibar_init, mu_init];

% Solve full ODE
[tps,y] = ode15s(@(t,y) ode_num(t,y,G,B,W_perp,W_par),[T_init T_max],Init,options);
theta = y(:,1);phi = y(:,2);psi = y(:,3);
% Solve the reduced ODE
[ t_asy,y_asy ] = Solve_reduced_ODE(G,B,w,Init_asy,tps);
alp = y_asy(:,1);phibar = y_asy(:,2);mu = y_asy(:,3);

nexttile;
plot_sphere(theta,phi,alp,phibar);

% c. Oblate object
B = B_oblate;
Init = Init_oblate;
T_max = T_max_slender(2);%time interval
Init_theta = Init(1);Init_phi = Init(2); Init_psi = Init(3);
[ alp_init, mu_init, phibar_init ] = Initial_conditions( Init_theta, Init_phi, Init_psi, W_perp, W_par );
Init_asy = [alp_init, phibar_init, mu_init];

% Solve full ODE
[tps,y] = ode15s(@(t,y) ode_num(t,y,G,B,W_perp,W_par),[T_init T_max],Init,options);
theta = y(:,1);phi = y(:,2);psi = y(:,3);
% Solve the reduced ODE
[ t_asy,y_asy ] = Solve_reduced_ODE(G,B,w,Init_asy,tps);
alp = y_asy(:,1);phibar = y_asy(:,2);mu = y_asy(:,3);

nexttile;
plot_sphere(theta,phi,alp,phibar);
toc

%% 3. General dynamics, low omega
tic

W_par = 15;
W_perp = 7.5;

% Quantities used in the analysis.
w = W_par/W_perp; %spinning ratio
lambda = sqrt(1 + w.^2); 

% a. Slender object
B = B_slender;
Init = Init_slender;
T_max = T_max_slender(3);%time interval
Init_theta = Init(1);Init_phi = Init(2); Init_psi = Init(3);
[ alp_init, mu_init, phibar_init ] = Initial_conditions( Init_theta, Init_phi, Init_psi, W_perp, W_par );
Init_asy = [alp_init, phibar_init, mu_init];

% Solve full ODE
[tps,y] = ode15s(@(t,y) ode_num(t,y,G,B,W_perp,W_par),[T_init T_max],Init,options);
theta = y(:,1);phi = y(:,2);psi = y(:,3);
% Solve the reduced ODE
[ t_asy,y_asy ] = Solve_reduced_ODE(G,B,w,Init_asy,tps);
alp = y_asy(:,1);phibar = y_asy(:,2);mu = y_asy(:,3);

nexttile;
plot_sphere(theta,phi,alp,phibar);

% b. Spheroidal object
B = B_spheroidal;
Init = Init_spheroidal;
T_max = T_max_non_slender(3);%time interval
Init_theta = Init(1);Init_phi = Init(2); Init_psi = Init(3);
[ alp_init, mu_init, phibar_init ] = Initial_conditions( Init_theta, Init_phi, Init_psi, W_perp, W_par );
Init_asy = [alp_init, phibar_init, mu_init];

% Solve full ODE
[tps,y] = ode15s(@(t,y) ode_num(t,y,G,B,W_perp,W_par),[T_init T_max],Init,options);
theta = y(:,1);phi = y(:,2);psi = y(:,3);
% Solve the reduced ODE
[ t_asy,y_asy ] = Solve_reduced_ODE(G,B,w,Init_asy,tps);
alp = y_asy(:,1);phibar = y_asy(:,2);mu = y_asy(:,3);

nexttile;
plot_sphere(theta,phi,alp,phibar);

% c. Oblate object
B = B_oblate;
Init = Init_oblate;
T_max = T_max_slender(3);%time interval
Init_theta = Init(1);Init_phi = Init(2); Init_psi = Init(3);
[ alp_init, mu_init, phibar_init ] = Initial_conditions( Init_theta, Init_phi, Init_psi, W_perp, W_par );
Init_asy = [alp_init, phibar_init, mu_init];

% Solve full ODE
[tps,y] = ode15s(@(t,y) ode_num(t,y,G,B,W_perp,W_par),[T_init T_max],Init,options);
theta = y(:,1);phi = y(:,2);psi = y(:,3);
% Solve the reduced ODE
[ t_asy,y_asy ] = Solve_reduced_ODE(G,B,w,Init_asy,tps);
alp = y_asy(:,1);phibar = y_asy(:,2);mu = y_asy(:,3);

nexttile;
plot_sphere(theta,phi,alp,phibar);

toc
%% 4. General dynamics, critical omega (= sqrt(2)), the averaged dynamics are those of a sphere
tic

W_par = 15;
W_perp = 15*sqrt(2);

% Quantities used in the analysis.
w = W_par/W_perp; %spinning ratio
lambda = sqrt(1 + w.^2); 

% a. Slender object
B = B_slender;
Init = Init_slender;
T_max = T_max_slender(4);%time interval
Init_theta = Init(1);Init_phi = Init(2); Init_psi = Init(3);
[ alp_init, mu_init, phibar_init ] = Initial_conditions( Init_theta, Init_phi, Init_psi, W_perp, W_par );
Init_asy = [alp_init, phibar_init, mu_init];

% Solve full ODE
[tps,y] = ode15s(@(t,y) ode_num(t,y,G,B,W_perp,W_par),[T_init T_max],Init,options);
theta = y(:,1);phi = y(:,2);psi = y(:,3);
% Solve the reduced ODE
[ t_asy,y_asy ] = Solve_reduced_ODE(G,B,w,Init_asy,tps);
alp = y_asy(:,1);phibar = y_asy(:,2);mu = y_asy(:,3);

nexttile;
plot_sphere(theta,phi,alp,phibar);

% b. Spheroidal object
B = B_spheroidal;
Init = Init_spheroidal;
T_max = T_max_non_slender(4);%time interval
Init_theta = Init(1);Init_phi = Init(2); Init_psi = Init(3);
[ alp_init, mu_init, phibar_init ] = Initial_conditions( Init_theta, Init_phi, Init_psi, W_perp, W_par );
Init_asy = [alp_init, phibar_init, mu_init];

% Solve full ODE
[tps,y] = ode15s(@(t,y) ode_num(t,y,G,B,W_perp,W_par),[T_init T_max],Init,options);
theta = y(:,1);phi = y(:,2);psi = y(:,3);
% Solve the reduced ODE
[ t_asy,y_asy ] = Solve_reduced_ODE(G,B,w,Init_asy,tps);
alp = y_asy(:,1);phibar = y_asy(:,2);mu = y_asy(:,3);

nexttile;
plot_sphere(theta,phi,alp,phibar);

% c. Oblate object
B = B_oblate;
Init = Init_oblate;
T_max = T_max_slender(4);%time interval
Init_theta = Init(1);Init_phi = Init(2); Init_psi = Init(3);
[ alp_init, mu_init, phibar_init ] = Initial_conditions( Init_theta, Init_phi, Init_psi, W_perp, W_par );
Init_asy = [alp_init, phibar_init, mu_init];

% Solve full ODE
[tps,y] = ode15s(@(t,y) ode_num(t,y,G,B,W_perp,W_par),[T_init T_max],Init,options);
theta = y(:,1);phi = y(:,2);psi = y(:,3);
% Solve the reduced ODE
[ t_asy,y_asy ] = Solve_reduced_ODE(G,B,w,Init_asy,tps);
alp = y_asy(:,1);phibar = y_asy(:,2);mu = y_asy(:,3);

nexttile;
plot_sphere(theta,phi,alp,phibar);

toc
%% 5. General dynamics, high omega. 
% The equivalent B is of opposite sign, but it is not easy to see because the new initial condition is close to the equator anyway.
tic

W_par = 15;
W_perp = 100;

% Quantities used in the analysis.
w = W_par/W_perp; %spinning ratio
lambda = sqrt(1 + w.^2); 

% a. Slender object
B = B_slender;
Init = Init_slender;
T_max = T_max_slender(5);%time interval
Init_theta = Init(1);Init_phi = Init(2); Init_psi = Init(3);
[ alp_init, mu_init, phibar_init ] = Initial_conditions( Init_theta, Init_phi, Init_psi, W_perp, W_par );
Init_asy = [alp_init, phibar_init, mu_init];

% Solve full ODE
[tps,y] = ode15s(@(t,y) ode_num(t,y,G,B,W_perp,W_par),[T_init T_max],Init,options);
theta = y(:,1);phi = y(:,2);psi = y(:,3);
% Solve the reduced ODE
[ t_asy,y_asy ] = Solve_reduced_ODE(G,B,w,Init_asy,tps);
alp = y_asy(:,1);phibar = y_asy(:,2);mu = y_asy(:,3);

nexttile;
plot_sphere(theta,phi,alp,phibar);

% b. Spheroidal object
B = B_spheroidal;
Init = Init_spheroidal;
T_max = T_max_non_slender(5);%time interval
Init_theta = Init(1);Init_phi = Init(2); Init_psi = Init(3);
[ alp_init, mu_init, phibar_init ] = Initial_conditions( Init_theta, Init_phi, Init_psi, W_perp, W_par );
Init_asy = [alp_init, phibar_init, mu_init];

% Solve full ODE
[tps,y] = ode15s(@(t,y) ode_num(t,y,G,B,W_perp,W_par),[T_init T_max],Init,options);
theta = y(:,1);phi = y(:,2);psi = y(:,3);
% Solve the reduced ODE
[ t_asy,y_asy ] = Solve_reduced_ODE(G,B,w,Init_asy,tps);
alp = y_asy(:,1);phibar = y_asy(:,2);mu = y_asy(:,3);

nexttile;
plot_sphere(theta,phi,alp,phibar);

% c. Oblate object
B = B_oblate;
Init = Init_oblate;
T_max = T_max_slender(5);%time interval
Init_theta = Init(1);Init_phi = Init(2); Init_psi = Init(3);
[ alp_init, mu_init, phibar_init ] = Initial_conditions( Init_theta, Init_phi, Init_psi, W_perp, W_par );
Init_asy = [alp_init, phibar_init, mu_init];

% Solve full ODE
[tps,y] = ode15s(@(t,y) ode_num(t,y,G,B,W_perp,W_par),[T_init T_max],Init,options);
theta = y(:,1);phi = y(:,2);psi = y(:,3);
% Solve the reduced ODE
[ t_asy,y_asy ] = Solve_reduced_ODE(G,B,w,Init_asy,tps);
alp = y_asy(:,1);phibar = y_asy(:,2);mu = y_asy(:,3);

nexttile;
plot_sphere(theta,phi,alp,phibar);

toc

%% Post-processing and save

tl.Padding = 'compact';
tl.TileSpacing = 'none';
% delete(findall(gcf,'type','annotation'))
% annotation('textbox',[.15 .89 .5 .1],'String','$B = 0.99$','EdgeColor','none','Interpreter','latex','FontSize',20)
% annotation('textbox',[.45 .89 .5 .1],'String','$B = 0.5$','EdgeColor','none','Interpreter','latex','FontSize',20)
% annotation('textbox',[.74 .89 .5 .1],'String','$B = -0.99$','EdgeColor','none','Interpreter','latex','FontSize',20)
% annotation('textbox',[.3 .86 .5 .1],'String','Classical Jeffery orbits (no spinning)','EdgeColor','none','Interpreter','latex','FontSize',20)
% annotation('textbox',[.33 .68 .5 .1],'String','Bacterial sublimit ($\Omega_{\bot} \ll \Omega_{\parallel}$)','EdgeColor','none','Interpreter','latex','FontSize',20)
% annotation('textbox',[.33 .5 .5 .1],'String','General dynamics, $\omega = 0.5$','EdgeColor','none','Interpreter','latex','FontSize',20)
% annotation('textbox',[.28 .32 .5 .1],'String','General dynamics, critical case ($\omega = \sqrt{2}$)','EdgeColor','none','Interpreter','latex','FontSize',20)
% annotation('textbox',[.35 .14 .5 .1],'String','General dynamics, $\omega \gg 1$','EdgeColor','none','Interpreter','latex','FontSize',20)

% exportgraphics(gcf,'figure_rotational_general_dynamics.png','Resolution',450)
%% Auxiliary functions %%%%%%%%

function dy = ode_num(t,y,G,B,w1,w2)

dy = zeros(size(y));

theta = y(1);
phi = y(2);
psi = y(3);

f1 = -B*(sin(2*theta)*sin(2*phi))/4;
f2 = (1 - B*cos(2*phi))/2;
f3 = (B/2)*cos(theta)*cos(2*phi);

dy(1) = w1*cos(psi) + G*f1;
dy(2) = w1*sin(psi)./sin(theta) + G*f2;
dy(3) = w2 - w1*sin(psi).*cot(theta) + G*f3;

end


function [ t_asy,y_asy ] = Solve_reduced_ODE(G,B,a,Init_asy,t)

%SOLVE_REDUCED_ODE Solve the reduced (averaged) ODE for ther 3D achiral
%system.

om = sqrt(1 + a.^2);

B_eff = B*(2*a^2 - 1)./(2*(a^2 + 1));

options = odeset('RelTol',1e-10,'AbsTol',1e-10);

[t_asy,y_asy] = ode15s(@(t,y) ode_num_reduced(t,y,G,B_eff,a,om),t,Init_asy,options);


    function dy = ode_num_reduced(t,y,G,B_eff,a,om)

        dy = zeros(size(y));

        alpha = y(1);
        phi_bar = y(2);
        %         mu = y(3);

        f1 = -B_eff*(sin(2*alpha)*sin(2*phi_bar))/4 ;
        f2 = (1 - B_eff*cos(2*phi_bar))/2 ;
        f3 = (B_eff/2)*cos(alpha)*cos(2*phi_bar) ;

        dy(1) = G*(f1);
        dy(2) = G*(f2);
        dy(3) = G*(f3);

    end

end

function [] = plot_sphere(theta,phi,alp,phibar)

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

% Full traj.
plot3(sin(theta).*cos(phi),sin(theta).*sin(phi),cos(theta),'b')

% Average traj.
plot3(sin(alp).*cos(phibar),sin(alp).*sin(phibar),cos(alp),'r','LineWidth',4)

% Mark the initial conditions. 
plot3(sin(theta(1)).*cos(phi(1)),sin(theta(1)).*sin(phi(1)),cos(theta(1)),'k.','MarkerSize',30)
plot3(sin(alp(1)).*cos(phibar(1)),sin(alp(1)).*sin(phibar(1)),cos(alp(1)),'k.','MarkerSize',30)

view(120,25)
camlight
axis equal off
grid on
box on

drawnow

end

function [] = plot_sphere_nospin(theta,phi)

% Orbit in the orientation space
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

% Plot the trajectory
plot3(sin(theta).*cos(phi),sin(theta).*sin(phi),cos(theta),'Color',[0.4 0.8 1],'LineWidth',4)
view(120,25);camlight;axis equal off;
% Initial condition
plot3(sin(theta(1)).*cos(phi(1)),sin(theta(1)).*sin(phi(1)),cos(theta(1)),'k.','MarkerSize',30)

drawnow

end
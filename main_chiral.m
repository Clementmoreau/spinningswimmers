%%%%%%%%%% Figure for wobbling swimmer with translational dynamics %%%%%%%%%
% (with chirality)

% Solves the full dynamics and plots the translational dynamics.

close all;

%% Setup.
% G = gamma, the shear rate of the flow.
G = 1;
% B = Bretherton constant.
B = 0.5;
% Ishimoto constant.
C = 0;
% Extra constants to account for translational chiral drift.
a1 = 0;
a2 = 0;
a3 = 1.5;

% W_par is the intrinsic spin of the swimmer about its axis of helicoidal symmetry.
W_par = 100;
% W_perp is the other direction of spin (perpendicular to the axis of helicoidal symmetry).
W_perp = 0.1;

% Swimming velocity along axis of helicoidal symmetry, e_hat_1.
V1 = 1;
% Swimming velocity along axis e_hat_2.
V2 = 0;
% Swimming velocity along axis e_hat_3.
V3 = 0;

% Quantities used in the asymptotic analysis.
w = W_perp / W_par; % spinning ratio.
lambda = sqrt(1 + w^2); 
% Effective Bretherton constant.
B_eff = B * (2 - w^2) / (2 * (1 + w^2));
% Effective chirality parameters.
C_eff_3 = C./lambda^3;
% Effective speed.
V_hat = abs((V1 + w*V2) / lambda);

% ODE options.
options = odeset('RelTol', 1e-11, 'AbsTol', 1e-11);

% Generate specific IC for theta, phi, and psi.
alp_c = asin(w/lambda);

init_theta = pi/4;
init_phi = -0.5;
init_psi = 3*pi/4;

% Initial condition for swimmer position.
X0 = [0;0;0];

% Generate IC for the long-time variables, alpha_bar, mu_bar, and phi_bar.
[init_alpha_bar, init_mu_bar, init_phi_bar] = Initial_conditions_2(init_theta, init_phi, init_psi, W_perp, W_par)

%init_phi_bar = init_phi_bar-pi;



% The IC vector for the full simulations.
init_full = [init_theta; init_phi; init_psi; X0];

% The IC vector for the reduced simulations.
init_reduced = [init_alpha_bar; init_phi_bar; init_mu_bar; X0];

% Time interval
T_min = 0;
T_max = 20;
tps = linspace(T_min,T_max,1e5);

params = struct();
params.G = G;
params.B = B;
params.C = C;
params.a1 = a1;
params.a2 = a2;
params.a3 = a3;
params.W_perp = W_perp;
params.W_par = W_par;
params.w = w;
params.lambda = lambda;
params.V = [V1; V2; V3];
params.V_hat = V_hat;

params.B_eff = B_eff;
params.C_eff3 = C_eff_3;

%% Solve full ODE.
tic
[~, sol_full] = ode15s(@(t,state) ode_full(t,state,params),tps,init_full,options);
toc

theta_full = sol_full(:,1);
phi_full = sol_full(:,2);
psi_full = sol_full(:,3);
x_full = sol_full(:,4);
y_full = sol_full(:,5);
z_full = sol_full(:,6);

%% Solve the reduced ODE.
tic
[~, sol_reduced] = ode15s(@(t,state) ode_reduced(t,state,params),tps,init_reduced,options);
toc

alpha_bar = sol_reduced(:,1);
phi_bar = sol_reduced(:,2);
mu_bar = sol_reduced(:,3);
x_bar = sol_reduced(:,4);
y_bar = sol_reduced(:,5);
z_bar = sol_reduced(:,6);

%% Plots.

T = W_par*lambda*tps';
a = 1/w;
% Define useful intermediate quantities.
w_cos_theta_asy = a*cos(alpha_bar) - sin(alpha_bar).*cos(T+mu_bar);
w_sin_theta_sin_psi_asy = cos(alpha_bar) + a*sin(alpha_bar).*cos(T+mu_bar);
tan_phi_minus_phibar = sin(T+mu_bar)./(cos(alpha_bar).*cos(T+mu_bar) + a*sin(alpha_bar));
phi_minus_phibar = atan2(sin(T+mu_bar),cos(alpha_bar).*cos(T+mu_bar) + a*sin(alpha_bar));

%%%% Plot omega*cos(theta) and compare to asymptotics
figure(8);clf;
set(gcf, 'Position',  [1, 1, 1500, 400])
tiledlayout(3,1)
nexttile;
plot(tps,lambda*a*cos(theta_full),'LineWidth',1)
hold on
plot(tps,w_cos_theta_asy,'--','LineWidth',1)
% xlim([T_max-1 T_max])

xlabel('$t$','Interpreter','latex','FontSize',18);
ylabel('$\omega \cos \theta$','Interpreter','latex','FontSize',18);
set(gcf,'Color',[1,1,1])
set(gca,'TickLabelInterpreter','latex','FontSize',18);
%%%

%%%% Plot omega*sin(theta)*sin(psi) and compare to asymptotics
nexttile;
plot(tps,lambda*a*sin(psi_full).*sin(theta_full),'LineWidth',1)
hold on
plot(tps,w_sin_theta_sin_psi_asy,'--','LineWidth',1)
% xlim([T_max-1 T_max])

xlabel('$t$','Interpreter','latex','FontSize',18);
ylabel('$\omega \sin \theta \sin \psi$','Interpreter','latex','FontSize',18);
set(gcf,'Color',[1,1,1])
set(gca,'TickLabelInterpreter','latex','FontSize',18);
%%%

%%% Plot cos(phi) [to bound it] and compare to asymptotics
nexttile;
plot(tps,cos(phi_full),'LineWidth',1)
hold on
plot(tps,cos(phi_minus_phibar + phi_bar),'--','LineWidth',1)

ylim([-1 1])
xlabel('$t$','Interpreter','latex','FontSize',18);
ylabel('$\cos (\phi)$','Interpreter','latex','FontSize',18);
set(gcf,'Color',[1,1,1])
set(gca,'TickLabelInterpreter','latex','FontSize',18);
%%%

figure(7); clf
set(gcf, 'Position',  [1, 600, 700, 700])
plot_sphere(theta_full,phi_full,alpha_bar,phi_bar);

figure(9); clf

set(gcf, 'Position',  [700, 600, 700, 700])

% Plot the average trajectory.
tick=100; %trick to make the line look smoother
plot3(x_bar(1:tick:end),y_bar(1:tick:end),z_bar(1:tick:end),'r','LineWidth',4);
hold on

% Plot the full trajectory.
plot3(x_full,y_full,z_full,'k','LineWidth',2);

% Ribbon plot representing the psi angle.
vertices = {[x_full, y_full, z_full]};
psi_mod = mod(psi_full, 2*pi) - pi;
twistangle = {[psi_mod(1);diff(psi_mod)]};
str = streamribbon(vertices, twistangle ,0.1);

% Graphic setup.
grid on
xlabel('x'); ylabel('y'); zlabel('z');
view(60,25)
shading interp
lighting gouraud
material dull
camlight
axis equal
col = othercolor('Bu_10');
col = [col(end:-1:1,:);col];
colormap(col)
colorbar

%% Auxiliary functions 
function d_state = ode_full(t,state,params)
% The RHS of the full ODEs for both the angular and translational
% dynamics. Parameters are passed via the structure params.
    d_state = zeros(6,1);

    % Unpack the states.
    theta = state(1);
    phi = state(2);
    psi = state(3);
    x = state(4);
    y = state(5);
    z = state(6);

    % Unpack the params.
    B = params.B;
    C = params.C;
    a1 = params.a1;
    a2 = params.a2;
    a3 = params.a3;
    G = params.G;
    w1 = params.W_perp;
    w2 = params.W_par;
    V1 = params.V(1);
    V2 = params.V(2);
    V3 = params.V(3);

    % Angular dynamics.
    f1 = -B * (sin(2*theta) * sin(2*phi)) / 4 - C*(sin(theta)*cos(2*phi))/2;
    f2 = (1 - B * cos(2*phi)) / 2 + (C/2)*sin(2*phi)*cos(theta);
    f3 = (B/2) * cos(theta) * cos(2*phi) - (C/2)*(cos(theta)).^2 .* sin(2*phi);

    d_state(1) = w1*cos(psi) + G*f1;
    d_state(2) = w1*sin(psi) ./ sin(theta) + G*f2;
    d_state(3) = w2 - w1*sin(psi).*cot(theta) + G*f3;

    % Translational dynamics.
    d_state(4) = V1 * cos(theta) + ...
                 V2 * sin(theta)*sin(psi) + ...
                 V3 * sin(theta)*cos(psi);

    d_state(5) = V1 * sin(phi)*sin(theta) + ...
                 V2 * ( cos(phi)*cos(psi) - cos(theta)*sin(phi)*sin(psi)) + ...
                 V3 * (-cos(phi)*sin(psi) - cos(theta)*sin(phi)*cos(psi));

    d_state(6) = -V1 * cos(phi)*sin(theta) + G*y + ...
                  V2 * ( sin(phi)*cos(psi) + cos(theta)*cos(phi)*sin(psi)) + ...
                  V3 * (-sin(phi)*sin(psi) + cos(theta)*cos(phi)*cos(psi));

    % Extra drift terms.
    g1 = -(a1-a3)*sin(2*phi)*sin(theta)*sin(2*theta)/4 +...
         a2*sin(theta)^2*cos(2*phi)/2;
    g2 = -(a1-a3)*sin(2*phi)*sin(theta)*sin(phi)*sin(theta)^2/2 -...
         a3*cos(phi)*sin(theta)/2 + ...
         a2*sin(2*theta)*sin(phi)/4;
    g3 = (a1-a3)*sin(2*phi)*sin(theta)*cos(phi)*sin(theta)^2/2 +...
         a3*sin(phi)*sin(theta)/2 + ...
         a2*sin(2*theta)*cos(phi)/4;

    d_state(4) = d_state(4) + G*g1;
    d_state(5) = d_state(5) + G*g2;
    d_state(6) = d_state(6) + G*g3;
    
end

function d_state = ode_reduced(t,state,params)
% The RHS of the reduced ODEs for both the auxilliary angular functions and
% average translational dynamics. Parameters are passed via the structure
% params.
    d_state = zeros(6,1);

    % Unpack the states.
    alpha_bar = state(1);
    phi_bar = state(2);
    mu_bar = state(3);
    x_bar = state(4);
    y_bar = state(5);
    z_bar = state(6);

    % Unpack the params.
    B_eff = params.B_eff;
    C = params.C;
    C_eff3 = params.C_eff3;
    a1 = params.a1;
    a2 = params.a2;
    a3 = params.a3;
    G = params.G;
    V_hat = params.V_hat;
    w = params.w;
    lambda = params.lambda;
    a = 1/w; %old notation for use inside I_hat and J_hat. also in those w is equal to sqrt(1+a^2). may be nice to change
   
    % Auxiliary angular dynamics.
    f1 = -B_eff * (sin(2*alpha_bar) * sin(2*phi_bar)) / 4 - C_eff3*(sin(alpha_bar)*cos(2*phi_bar))/2;
    f2 = (1 - B_eff * cos(2*phi_bar)) / 2 + C_eff3*sin(2*phi_bar)*cos(alpha_bar)/2;
    f3 = (B_eff/2) * cos(alpha_bar) * cos(2*phi_bar) - C_eff3*(cos(alpha_bar)).^2 .* sin(2*phi_bar)/2;

    % Extra terms.
    f3_extra = (- 3*w^2./2./lambda^3)*sin(alpha_bar)^2;
        
    d_state(1) = G*(f1);
    d_state(2) = G*(f2);
    d_state(3) = G*(f3 + (C/2)*sin(2*phi_bar)*f3_extra);

    % Average translational dynamics.
    d_state(4) = V_hat * cos(alpha_bar);
    d_state(5) = V_hat * sin(phi_bar) * sin(alpha_bar);
    d_state(6) = -V_hat * cos(phi_bar) * sin(alpha_bar) + G*y_bar;

    % Extra drift terms.
    A1 = G/4/lambda^2*(-(a1*w^2+a3)/lambda*2*sin(alpha_bar)*cos(2*phi_bar)+a2*(1-w^2/2)*sin(2*phi_bar)*sin(2*alpha_bar));
    A2 = G/4/lambda^2*((a1*w^2+a3)/lambda*sin(2*phi_bar)*sin(2*alpha_bar)+a2*(1-w^2/2)*2*sin(alpha_bar)*cos(2*phi_bar));
    A3 = G/4/lambda^3*((a1-3*a3)*w^2-2*a1)*sin(2*phi_bar)*sin(alpha_bar)^2;

    d_state(4) = d_state(4) + A2*sin(alpha_bar) + A3*cos(alpha_bar);
    d_state(5) = d_state(5) + A1*cos(phi_bar) - A2*sin(phi_bar)*cos(alpha_bar) + A3*sin(phi_bar)*sin(alpha_bar);
    d_state(6) = d_state(6) + A1*sin(phi_bar) + A2*cos(phi_bar)*cos(alpha_bar) - A3*cos(phi_bar)*sin(alpha_bar);
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
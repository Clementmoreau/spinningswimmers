%%%%%%%%%% Movie for rotational dynamics for simple shape %%%%%%%%%
% Solves the dynamics and saves the data for the movie

clear all

% Setup.
% G = gamma, the shear rate of the flow.
G = 1;

% B = Bretherton parameter.
B = 0.99;

% Ishimoto parameter.
C = 0;

% Extra chiral parameter.
D=0;

% Generate specific IC for theta, phi, and psi.
init_theta = pi/18;
init_phi = pi/2;
init_psi = -pi/4;

% The IC vector for the full simulations.
init_full = [init_theta; init_phi; init_psi];

% ODE options.
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);

% Time interval
T_min = 0;
T_max = 100;
tps = linspace(T_min,T_max,3e4);

% W_perp is the other direction of spin (perpendicular to the axis of helicoidal symmetry)
% and takes variable values here.
W_perp_var = [0.5,1.5,7.5,15*sqrt(2),100];

data_movie_rot_achiral = cell(length(W_perp_var)+1,9);
% the data contains:
% w_par
% w_perp
% B_hat
% theta_full
% phi_full
% psi_full
% alpha_bar
% phi_bar
% mu_bar

% first of all, the orbit without spinningz
W_par = 0;
W_perp = 0;

params = struct();
params.W_perp = W_perp;
params.W_par = W_par;
params.G = G;
params.B = B;
params.C = C;
params.D = D;

[tps,y] = ode15s(@(t,state) ode_full(t,state,params),tps,init_full,options);% Solve full ODE
theta_nospin= y(:,1);phi_nospin= y(:,2);psi_nospin= y(:,3);

data_movie_rot_achiral{1,4} = theta_nospin;
data_movie_rot_achiral{1,5} = phi_nospin;
data_movie_rot_achiral{1,6} = psi_nospin;

W_par = 15;

for i = 1:length(W_perp_var)

    W_perp = W_perp_var(i);

    % Quantities used in the asymptotic analysis.
    w = W_perp / W_par; % spinning ratio.
    lambda = sqrt(1 + w^2);

    % Effective Bretherton constant.
    B_eff = B * (2 - w^2) / (2 * (1 + w^2));
    % Effective chirality parameters.
    C_eff = (C+w^2*D)./lambda^3;
    D_eff = (3*w^2*C+(2-w^2)*D)/2/lambda^3;


    params.W_perp = W_perp;
    params.W_par = W_par;
    params.w = w;
    params.lambda = lambda;

    params.B_eff = B_eff;
    params.C_eff = C_eff;
    params.D_eff = D_eff;

    % Generate IC for the long-time variables, alpha_bar, mu_bar, and phi_bar.
[init_alpha_bar, init_mu_bar, init_phi_bar] = Initial_conditions_2(init_theta, init_phi, init_psi, W_perp, W_par);

    % Solve full ODE.
    tic
    [~, sol_full] = ode15s(@(t,state) ode_full(t,state,params),tps,init_full,options);
    toc

    theta_full = sol_full(:,1);
    phi_full = sol_full(:,2);
    psi_full = sol_full(:,3);

    % The IC vector for the reduced simulations.
init_reduced = [init_alpha_bar; init_phi_bar; init_mu_bar];

    % Solve the reduced ODE.
    tic
    [~, sol_reduced] = ode15s(@(t,state) ode_reduced(t,state,params),tps,init_reduced,options);
    toc

    alpha_bar = sol_reduced(:,1);
    phi_bar = sol_reduced(:,2);
    mu_bar = sol_reduced(:,3);

    % save the data
    data_movie_rot_achiral{i+1,1} = W_par;
    data_movie_rot_achiral{i+1,2} = W_perp;
    data_movie_rot_achiral{i+1,3} = B_eff;
    data_movie_rot_achiral{i+1,4} = theta_full;
    data_movie_rot_achiral{i+1,5} = phi_full;
    data_movie_rot_achiral{i+1,6} = psi_full;
    data_movie_rot_achiral{i+1,7} = alpha_bar;
    data_movie_rot_achiral{i+1,8} = phi_bar;
    data_movie_rot_achiral{i+1,9} = mu_bar;

end

% Save the data.
save('data_movie_1.mat')

%% Auxiliary functions
function d_state = ode_full(t,state,params)
% The RHS of the full ODEs for both the angular and translational
% dynamics. Parameters are passed via the structure params.
d_state = zeros(3,1);

% Unpack the states.
theta = state(1);
phi = state(2);
psi = state(3);

% Unpack the params.
B = params.B;
C = params.C;
D = params.D;
G = params.G;
w1 = params.W_perp;
w2 = params.W_par;


% Angular dynamics.
f1 = -B * (sin(2*theta) * sin(2*phi)) / 4 - C*(sin(theta)*cos(2*phi))/2;
f2 = (1 - B * cos(2*phi)) / 2 + (C/2)*sin(2*phi)*cos(theta);
f3 = (B/2) * cos(theta) * cos(2*phi) - (C/2)*(cos(theta)).^2 .* sin(2*phi) - (D/2)*sin(theta).^2*sin(2*phi);

d_state(1) = w1*cos(psi) + G*f1;
d_state(2) = w1*sin(psi) ./ sin(theta) + G*f2;
d_state(3) = w2 - w1*sin(psi).*cot(theta) + G*f3;

end

function d_state = ode_reduced(t,state,params)
% The RHS of the reduced ODEs for both the auxilliary angular functions and
% average translational dynamics. Parameters are passed via the structure
% params.
d_state = zeros(3,1);

% Unpack the states.
alpha_bar = state(1);
phi_bar = state(2);
mu_bar = state(3);

% Unpack the params.
B_eff = params.B_eff;
C_eff = params.C_eff;
D_eff = params.D_eff;
G = params.G;

% Auxiliary angular dynamics.
f1 = -B_eff * (sin(2*alpha_bar) * sin(2*phi_bar)) / 4 - C_eff*(sin(alpha_bar)*cos(2*phi_bar))/2;
f2 = (1 - B_eff * cos(2*phi_bar)) / 2 + C_eff*sin(2*phi_bar)*cos(alpha_bar)/2;
f3 = (B_eff/2) * cos(alpha_bar) * cos(2*phi_bar) - C_eff*(cos(alpha_bar)).^2 .* sin(2*phi_bar)/2 - (D_eff/2)*sin(alpha_bar).^2*sin(2*phi_bar);

d_state(1) = G*(f1);
d_state(2) = G*(f2);
d_state(3) = G*(f3);

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
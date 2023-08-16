% Emergent three-dimensional dynamics of rapidly spinning microswimmers in shear flow. Part I: Achiral bodies

% Figure for translational dynamics, general dynamics

close all

%% Setup.
% G = gamma, the shear rate of the flow.
G = 1;
% B = Bretherton constant.
B = 0.5;

% W_par is the intrinsic spin of the swimmer about its axis of helicoidal symmetry.
W_par = 10;
% W_perp is the other direction of spin (perpendicular to the axis of helicoidal symmetry).
W_perp = 0.01;

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
% Effective speed.
V_hat = (V1 + w*V2) / lambda;

% ODE options.
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);

% Generate specific IC for theta, phi, and psi.
init_theta = 2*pi/5;
init_phi = -7*pi/15;
init_psi = pi/2;

% Initial condition for swimmer position.
X0 = [0;0;0];

% Generate IC for the long-time variables, alpha_bar, mu_bar, and phi_bar.
[init_alpha_bar, init_mu_bar, init_phi_bar] = compute_initial_conditions(init_theta, init_phi, init_psi, W_perp, W_par);

% The IC vector for the full simulations.
init_full = [init_theta; init_phi; init_psi; X0];

% The IC vector for the reduced simulations.
init_reduced = [init_alpha_bar; init_phi_bar; init_mu_bar; X0];

% Time interval
T_min = 0;
T_max = 9;
tps = linspace(T_min,T_max,1e5);

params = struct();
params.G = G;
params.B = B;
params.W_perp = W_perp;
params.W_par = W_par;
params.w = w;
params.B_eff = B_eff;
params.lambda = lambda;
params.V = [V1; V2; V3];
params.V_hat = V_hat;

%Prepare the figure.
F1 = figure(1);clf
set(gcf, 'Position',  [1, 640, 1200, 680])
tl = tiledlayout(2,3);
graphic_params = struct();
ribw = 1/2;
graphic_params.ribbonwidth = ribw;
graphic_params.v = [-287,17];
v=[-287,17];

V2_var = [0 0.25 sqrt(2)/2 1 2];
W_perp_var = [-40 -sqrt(2)*10 -5 5 sqrt(2)*10 40];

i_plot = [1 2 3 6 5 4];

%%
for i_w = 1:length(W_perp_var)

    i_p = i_plot(i_w);
    
    W_perp = W_perp_var(i_w);
    w = W_perp / W_par; % spinning ratio.
    lambda = sqrt(1 + w^2);
    B_eff = B * (2 - w^2) / (2 * (1 + w^2));% Effective Bretherton constant.
    
    ribbonwidth_var = ribw*ones(1,length(V2_var));
    
    nexttile(i_p);
    
    for i = 1:length(V2_var)
        
        params.V = [V1; V2_var(i); V3];
        V_hat = (V1 + w*V2_var(i)) / lambda;% Effective speed.
        
        % Generate IC for the long-time variables, alpha_bar, mu_bar, and phi_bar.
        [init_alpha_bar, init_mu_bar, init_phi_bar] = compute_initial_conditions(init_theta, init_phi, init_psi, W_perp, W_par);
        
        % The IC vector for the reduced simulations.
        init_reduced = [init_alpha_bar; init_phi_bar; init_mu_bar; X0];
        
        params.w = w;
        params.W_perp = W_perp;
        params.W_par = W_par;
        params.B_eff = B_eff;
        params.lambda = lambda;
        params.V_hat = V_hat;
        
        % Solve full ODE.
        tic
        [~, sol_full] = ode15s(@(t,state) ode_full(t,state,params),tps,init_full,options);
        toc
        
        theta_full = sol_full(:,1);
        phi_full = sol_full(:,2);
        psi_full = sol_full(:,3);
        x_full = sol_full(:,4);
        y_full = sol_full(:,5);
        z_full = sol_full(:,6);
        
        % Solve the reduced ODE.
        tic
        [~, sol_reduced] = ode15s(@(t,state) ode_reduced(t,state,params),tps,init_reduced,options);
        toc
        
        alpha_bar = sol_reduced(:,1);
        phi_bar = sol_reduced(:,2);
        mu_bar = sol_reduced(:,3);
        
        x_bar = sol_reduced(:,4);
        y_bar = sol_reduced(:,5);
        z_bar = sol_reduced(:,6);
        
        % Ribbon plot.
        graphic_params.ribbonwidth = ribbonwidth_var(i);
        plot_traj_ribbon(x_bar,z_bar,y_bar,x_full,z_full,y_full,psi_full,graphic_params,i)
        hold on
        
    end
    
    % Post-processing
    box on
    % camlight(-30,15)
    camlight('headlight')
    camlight(-20,20)
    axis equal
    set(gca, 'YDir', 'reverse')
    %axis([-1 4.1 -3 4 -12 3])
    col = colormap(gray(75));
    col = [col(50:end,:);col(end:-1:50,:)];
    colormap(col)
    % colorbar
    set(gca,'FontSize',20)
    set(gca,'TickLabelInterpreter','latex')
    %zticks([0 8])
    xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex'); zlabel('$y$','Interpreter','latex');
    view(v)
    drawnow
    if i_w == 2
        title(['$\omega = -\sqrt{2}$'],'Interpreter','latex')
    elseif i_w == 5
        title(['$\omega = \sqrt{2}$'],'Interpreter','latex')
    else
        title(['$\omega = ',num2str(w),'$'],'Interpreter','latex')
    end
    
end

%%

tl.Padding = 'tight';
tl.TileSpacing = 'tight';
exportgraphics(gcf,'figure_translational_general_dynamics_V.png','Resolution',500)
% exportgraphics(gcf,'figure_translational_general dynamics_3.eps','ContentType','vector')

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
G = params.G;
w1 = params.W_perp;
w2 = params.W_par;
V1 = params.V(1);
V2 = params.V(2);
V3 = params.V(3);

% Angular dynamics.
f1 = -B * (sin(2*theta) * sin(2*phi)) / 4;
f2 = (1 - B * cos(2*phi)) / 2;
f3 = (B/2) * cos(theta) * cos(2*phi);

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
G = params.G;
V_hat = params.V_hat;

% Auxiliary angular dynamics.
f1 = -B_eff * (sin(2*alpha_bar) * sin(2*phi_bar)) / 4;
f2 = (1 - B_eff * cos(2*phi_bar)) / 2;
f3 = (B_eff/2) * cos(alpha_bar) * cos(2*phi_bar);

d_state(1) = G * f1;
d_state(2) = G * f2;
d_state(3) = G * f3;

% Average translational dynamics.
d_state(4) = V_hat * cos(alpha_bar);
d_state(5) = V_hat * sin(phi_bar) * sin(alpha_bar);
d_state(6) = -V_hat * cos(phi_bar) * sin(alpha_bar) + G*y_bar;
end

function [] = plot_traj_ribbon(x_bar,y_bar,z_bar,x_full,y_full,z_full,psi_full,graphic_params,i)

% Plot the avg trajectory.
tick=100;
col_v = colormap(autumn(5));
plot3(x_bar(1:tick:end),y_bar(1:tick:end),z_bar(1:tick:end),'Color',col_v(i,:),'LineWidth',3);

hold on

ribbonwidth = graphic_params.ribbonwidth;
v = graphic_params.v;

% Plot the full trajectory.
tick2=10;
plot3(x_full(1:tick2:end),y_full(1:tick2:end),z_full(1:tick2:end),'k','LineWidth',2);

% Plot the ribbon representing the evolution pf psi.
vertices = {[x_full, y_full, z_full]};
psi_mod = mod(psi_full, 2*pi) - pi;
twistangle = {[psi_mod(1);diff(psi_mod)]};
str = streamribbon(vertices, twistangle, ribbonwidth);
str.FaceAlpha = 0.3;

grid on
shading interp
lighting gouraud
material dull


end
%%%%%%%%%% Figure for wobbling swimmer with translational dynamics %%%%%%%%%
% (without chirality)

% Figure for translational dynamics, bacterial limit

close all;

%% Setup.
% G = gamma, the shear rate of the flow.
G = 1;
% B = Bretherton constant.
B = 0.5;

% W_par is the intrinsic spin of the swimmer about its axis of helicoidal symmetry.
W_par = 15;
% W_perp is the other direction of spin (perpendicular to the axis of helicoidal symmetry).
W_perp = 0.0001;

% Swimming velocity along axis of helicoidal symmetry, e_hat_1.
V1 = 1;
% Swimming velocity along axis e_hat_2.
V2 = 1;
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
[init_alpha_bar, init_mu_bar, init_phi_bar] = Initial_conditions(init_theta, init_phi, init_psi, W_perp, W_par);

% The IC vector for the full simulations.
init_full = [init_theta; init_phi; init_psi; X0];

% The IC vector for the reduced simulations.
init_reduced = [init_alpha_bar; init_phi_bar; init_mu_bar; X0];

% Time interval
T_min = 0;
T_max = 25;
tps = linspace(T_min,T_max,3e4);

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
figure(2);clf
set(gcf, 'Position',  [100, 100, 1000, 1000])
tiledlayout(2,3)
graphic_params = struct();
graphic_params.ribbonwidth = 0.8;
graphic_params.v = [60,25];

% In this figure, we want to show that V_2 and V_3 have no influence on the
% average trajectory.

%% First, let's plot a trajectory with no spinning.

% Redefine the parameters.
W_par = 0;
W_perp = 0;
V1 = 1; V2 = 0 ; V3 = 0;
params.V = [V1; V2; V3];
params.W_perp = W_perp;
params.W_par = W_par;
w = W_perp / W_par; % spinning ratio.
lambda = sqrt(1 + w^2); 
% Effective Bretherton constant.
B_eff = B * (2 - w^2) / (2 * (1 + w^2));
% Effective speed.
V_hat = (V1 + w*V2) / lambda;
params.w = w;
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

% Ribbon plot.
nexttile(1);
plot_traj_ribbon_nospin(x_full,y_full,z_full,psi_full,graphic_params)

%% Now, some spinning. 
% We see that "nothing changes" in terms of full trajectory or average dynamics.
% The only thing that changes is \psi : the ribbon is now twisted.

W_par = 15;
W_perp = 0;
params.W_perp = W_perp;
params.W_par = W_par;
w = W_perp / W_par; % spinning ratio.
lambda = sqrt(1 + w^2); 
% Effective Bretherton constant.
B_eff = B * (2 - w^2) / (2 * (1 + w^2));
% Effective speed.
V_hat = (V1 + w*V2) / lambda;
params.w = w;
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
nexttile(2);
plot_traj_ribbon(x_bar,y_bar,z_bar,x_full,y_full,z_full,psi_full,graphic_params)

%% And now we add V2.
% The full dynamics trajectory changes a little, but the avg dyamics stay
% the same. 

V1 = 1; V2 = 1; V3 = 0;
params.V = [V1; V2; V3];
V_hat = (V1 + w*V2) / lambda;
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
nexttile(3);
plot_traj_ribbon(x_bar,y_bar,z_bar,x_full,y_full,z_full,psi_full,graphic_params)

%% And a small V3.
% The full dynamics trajectory changes a little, but the avg dyamics stay
% the same. 

V1 = 1; V2 = 1; V3 = 1;
params.V = [V1; V2; V3];
V_hat = (V1 + w*V2) / lambda;
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
nexttile(4);
plot_traj_ribbon(x_bar,y_bar,z_bar,x_full,y_full,z_full,psi_full,graphic_params)

%% Large V2.
% The full dynamics trajectory changes a lot, but the avg dyamics stay the same. 

V1 = 1; V2 = 5; V3 = 0;
params.V = [V1; V2; V3];
V_hat = (V1 + w*V2) / lambda;
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
nexttile(5);
plot_traj_ribbon(x_bar,y_bar,z_bar,x_full,y_full,z_full,psi_full,graphic_params)

%% Large V2 and V3.
% The full dynamics trajectory changes a lot, but the avg dyamics still stay the same! 

V1 = 1; V2 = 5; V3 = 5;
params.V = [V1; V2; V3];
V_hat = (V1 + w*V2) / lambda;
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
nexttile(6);
plot_traj_ribbon(x_bar,y_bar,z_bar,x_full,y_full,z_full,psi_full,graphic_params)

%% Some graphical setup.

c=colorbar;
c.Location='manual';
c.Position=[0.92 0.35 0.01 0.35];
c.TickLabelInterpreter='latex';
c.FontSize=24;
c.TickLength=0.01;
c.Ticks=[-3.141, 0, 3.141];
c.TickLabels={'$-\pi$','0','$\pi$'}; 
c.Label.String='$\psi$';
c.Label.Interpreter='latex';
c.Label.Rotation=0;
c.Label.Position=[0.5 3.7 0];
c.Label.FontSize=24;

%% Save the figure. 

%exportgraphics(gcf,'figure_translational_bacterial_limit.png','Resolution',500)
%exportgraphics(gcf,'figure_translational_bacterial_limit.eps','ContentType','vector')

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

function [] = plot_traj_ribbon(x_bar,y_bar,z_bar,x_full,y_full,z_full,psi_full,graphic_params)
    
    % Plot the averaged and full trajectories.
    tick=100; % trick to make the line look better
    plot3(x_bar(1:tick:end),y_bar(1:tick:end),z_bar(1:tick:end),'r','LineWidth',4)
    hold on
    
    tick2=10;
    plot3(x_full(1:tick2:end),y_full(1:tick2:end),z_full(1:tick2:end),'k','LineWidth',2);
    
    % Prepare the ribbon data.
    vertices = {[x_full, y_full, z_full]};
    psi_mod = mod(psi_full, 2*pi) - pi;
    twistangle = {[psi_mod(1);diff(psi_mod)]};
    % Unpack the parameters and plot the ribbon.
    ribbonwidth = graphic_params.ribbonwidth;
    v = graphic_params.v;
    str = streamribbon(vertices, twistangle, ribbonwidth);
    
    % Post-processing.
    xlabel('x','Interpreter','latex'); ylabel('y','Interpreter','latex'); zlabel('z','Interpreter','latex');
    grid on
    shading interp
    lighting gouraud
    material dull
    view(v)
    camlight(-30,15)
    camlight('headlight')
    % camlight(10,10)
    axis equal
    axis([-1 6.5 -3 3.5 -14 3])
    col = othercolor('Bu_10');
    col = colormap(gray(100));
    col = [col(50:end,:);col(end:-1:50,:)];
    colormap(col)
    % colorbar
    set(gca,'FontSize',20)
    set(gca,'TickLabelInterpreter','latex')
    %zticks([0 8])

end

function [] = plot_traj_ribbon_nospin(x_full,y_full,z_full,psi_full,graphic_params)
    tick=100; % BW: not sure this is needed here?
    plot3(x_full,y_full,z_full,'k','LineWidth',2);
    grid on
    xlabel('x','Interpreter','latex'); ylabel('y','Interpreter','latex'); zlabel('z','Interpreter','latex');
    hold on

    ribbonwidth = graphic_params.ribbonwidth;
    v = graphic_params.v;
    
    vertices = {[x_full, y_full, z_full]};
    psi_mod = mod(psi_full, 2*pi) - pi;
    twistangle = {[psi_mod(1);diff(psi_mod)]};
    str = streamribbon(vertices, twistangle, ribbonwidth);
    
    shading interp
    lighting gouraud
    material dull
    view(v)
    camlight(20,45)
    camlight('headlight')
    % camlight(10,10)
    axis equal
    axis([-1 6.5 -3 3.5 -14 3])
    %col = othercolor('Bu_10');
    col = colormap(gray(100));
    col = [col(end:-1:1,:);col];
    colormap(col)
    % colorbar

    set(gca,'FontSize',20)
    set(gca,'TickLabelInterpreter','latex')
    %zticks([0 8])
end
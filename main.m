%%%%%%%%%% Figure for wobbling swimmer with translational dynamics %%%%%%%%%
% (without chirality)

% Solves the full dynamics and plots the translational dynamics.

close all;

%% Setup.
% G = gamma, the shear rate of the flow.
G = 1;
% B = Bretherton constant.
B = 0.5;

% W_par is the intrinsic spin of the swimmer about its axis of helicoidal symmetry.
W_par = 10;
% W_perp is the other direction of spin (perpendicular to the axis of helicoidal symmetry).
W_perp = -100;

% Swimming velocity along axis of helicoidal symmetry, e_hat_1.
V1 = 1;
% Swimming velocity along axis e_hat_2.
V2 = 2;
% Swimming velocity along axis e_hat_3.
V3 = 0;

% Quantities used in the asymptotic analysis.
w = W_perp / W_par; % spinning ratio.
lambda = sqrt(1 + w^2); 
% Effective Bretherton constant.
B_eff = B * (2 - w^2) / (2 * (1 + w^2));
% Effective speed.
V_hat = abs((V1 + w*V2) / lambda);

% ODE options.
options = odeset('RelTol', 1e-11, 'AbsTol', 1e-11);

% Generate specific IC for theta, phi, and psi.
init_theta = pi/2+0.5;
init_phi = pi/2;
init_psi = pi/2;

% Initial condition for swimmer position.
X0 = [0;0;0];

% Generate IC for the long-time variables, alpha_bar, mu_bar, and phi_bar.
[init_alpha_bar, init_mu_bar, init_phi_bar] = Initial_conditions_2(init_theta, init_phi, init_psi, W_perp, W_par);

if W_perp<0
    init_alpha_bar = pi-init_alpha_bar;
    init_phi_bar = pi-init_phi_bar;
    init_mu_bar = pi+init_mu_bar;
end


% The IC vector for the full simulations.
init_full = [init_theta; init_phi; init_psi; X0];

% The IC vector for the reduced simulations.
init_reduced = [init_alpha_bar; init_phi_bar; init_mu_bar; X0];

% Time interval
T_min = 0;
T_max = 1;
tps = linspace(T_min,T_max,1e4);

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

figure(7); clf

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

% % Components of position.
figure(2); clf
plot(tps,x_full,tps,x_bar,tps,y_full,tps,y_bar,tps,z_full,tps,z_bar)
legend('x','x avg','y','y avg','z','z avg')

T = W_par*lambda*tps';
a = 1/w;
% Define useful intermediate quantities.
w_cos_theta_asy = a*cos(alpha_bar) - sin(alpha_bar).*cos(T+mu_bar);
w_sin_theta_sin_psi_asy = cos(alpha_bar) + a*sin(alpha_bar).*cos(T+mu_bar);
tan_phi_minus_phibar = sin(T+mu_bar)./(cos(alpha_bar).*cos(T+mu_bar) + a*sin(alpha_bar));
phi_minus_phibar = atan2(sin(T+mu_bar),cos(alpha_bar).*cos(T+mu_bar) + a*sin(alpha_bar));

%%%% Plot omega*cos(theta) and compare to asymptotics
figure(8);clf;
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
figure(9);clf;
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
figure(10);clf;
plot(tps,cos(phi_full),'LineWidth',1)
hold on
plot(tps,cos(phi_minus_phibar + phi_bar),'--','LineWidth',1)

ylim([-1 1])
xlabel('$t$','Interpreter','latex','FontSize',18);
ylabel('$\cos (\phi)$','Interpreter','latex','FontSize',18);
set(gcf,'Color',[1,1,1])
set(gca,'TickLabelInterpreter','latex','FontSize',18);
%%%


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
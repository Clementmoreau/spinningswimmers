%%%%%%%%%% Figure for swimmer with translational dynamics %%%%%%%%%
% (with chirality, without spinning)

% Solves the full dynamics and plots the translational dynamics.
% A possible figure showing the influence of each parameter

close all;

% Introduction

% G = gamma, the shear rate of the flow.
G = 1;
% B = Bretherton constant.
B = 0.9;
% Swimming velocity along axis of helicoidal symmetry, e_hat_1.
V1 = 1;
% Swimming velocity along axis e_hat_2.
V2 = 0;
% Swimming velocity along axis e_hat_3.
V3 = 0;
% Time interval
T_min = 0;
T_max = 100;
tps = linspace(T_min,T_max,2.5e2);
% Initial conditions.
init_theta = pi/3;
init_phi = pi/6;
init_psi = 2*pi/3;
% Initial condition for swimmer position.
X0 = [0;0;0];
% The IC vector for the simulations.
init_full = [init_theta; init_phi; init_psi; X0];
% Figure
figure(9); clf
set(gcf, 'Position',  [1, 500, 1050, 900])
tl=tiledlayout(4,4);

% Then, let's show the influence of each parameter.
% Note that D does not have any influence on the translational dynamics --
% it would only have one indirectly if we add spinning, by changing the
% effective C.

% ----- beta_hat -----
beta_var = -0.5:0.2:0.5;
col_beta = colormap(parula(256));
N_col = 256;

for i = 1:length(beta_var)
[theta,phi,psi,x,y,z] = traj_full(init_full,G,B,0,0,beta_var(i),0,0,V1,V2,V3,tps);

% Setup colour gradient.
i_col = floor(abs(beta_var(i))/beta_var(end)*(N_col-1))+1;

% Plot the trajectory.
nexttile(1);
plot3(x,y,z,'LineWidth',1,'Color',col_beta(i_col,:));
hold on
nexttile(2);
plot(tps,x,'r','LineWidth',1,'Color',col_beta(i_col,:));
hold on
nexttile(3);
plot(tps,y,'r','LineWidth',1,'Color',col_beta(i_col,:));
hold on
nexttile(4);
plot(tps,z,'r','LineWidth',1,'Color',col_beta(i_col,:));
hold on
colormap(gca,col_beta)
c_beta = colorbar;
drawnow
end

% ----- gamma_hat ----- 
gamma_var = [-1.5:0.3:-0.3,0.3:0.3:1.5];
col_gamma = colormap(parula(256));
N_col = 256;

for i = 1:length(gamma_var)
[theta,phi,psi,x,y,z] = traj_full(init_full,G,B,0,0, 0,gamma_var(i),0,V1,V2,V3,tps);

% Setup colour gradient.
i_col = floor(abs(gamma_var(i))/gamma_var(end)*(N_col-1))+1;

% Plot the trajectory.
nexttile(5);
plot3(x,y,z,'LineWidth',1,'Color',col_gamma(i_col,:));hold on
nexttile(6);
plot(tps,x,'LineWidth',1,'Color',col_gamma(i_col,:));hold on
nexttile(7);
plot(tps,y,'LineWidth',1,'Color',col_gamma(i_col,:));hold on
nexttile(8);
plot(tps,z,'LineWidth',1,'Color',col_gamma(i_col,:));hold on
colormap(gca,col_gamma)
c_gamma = colorbar;
drawnow
end

% ----- delta_hat -----
delta_var = [-7:-1,1:7];
col_delta = colormap(parula(256));
N_col = 256;

for i = 1:length(delta_var)
[theta,phi,psi,x,y,z] = traj_full(init_full,G,B,0,0,0,0,delta_var(i),V1,V2,V3,tps);

% Setup colour gradient.
i_col = floor(abs(delta_var(i))/delta_var(end)*(N_col-1))+1;

% Plot the trajectory.
nexttile(9);
plot3(x,y,z,'LineWidth',1,'Color',col_delta(i_col,:));hold on
nexttile(10);
plot(tps,x,'LineWidth',1,'Color',col_delta(i_col,:));hold on
nexttile(11);
plot(tps,y,'LineWidth',1,'Color',col_delta(i_col,:));hold on
nexttile(12);
plot(tps,z,'LineWidth',1,'Color',col_delta(i_col,:));hold on
colormap(gca,col_delta)
c_delta = colorbar;

drawnow
end

% ----- C_hat -----
C_var = [-0.15:0.03:-0.03,0.03:0.03:0.15];
col_C = colormap(parula(256));
N_col = 256;

for i = 1:length(C_var)
[theta,phi,psi,x,y,z] = traj_full(init_full,G,B,C_var(i),0,0,0,0,V1,V2,V3,tps);

% Setup colour gradient.
i_col = floor(abs(C_var(i))/C_var(end)*(N_col-1))+1;

% Plot the trajectory.
nexttile(13);
plot3(x,y,z,'LineWidth',1,'Color',col_C(i_col,:));hold on
nexttile(14);
plot(tps,x,'LineWidth',1,'Color',col_C(i_col,:));hold on
nexttile(15);
plot(tps,y,'LineWidth',1,'Color',col_C(i_col,:));hold on
nexttile(16);
plot(tps,z,'LineWidth',1,'Color',col_C(i_col,:));hold on
colormap(gca,col_C)
c_C = colorbar;
drawnow
end

% Finally, the traj without any chirality or asymmetry parameters.
[theta,phi,psi,x,y,z] = traj_full(init_full,G,B,0,0,0,0,0,V1,V2,V3,tps);
for i = 1:4
% Plot the average trajectory.
nexttile(4*(i-1)+1);
plot3(x,y,z,'k','LineWidth',2);
nexttile(4*(i-1)+2);
plot(tps,x,'k','LineWidth',2);
nexttile(4*(i-1)+3);
plot(tps,y,'k','LineWidth',2);
nexttile(4*(i-1)+4);
plot(tps,z,'k','LineWidth',2);
end

% Plenty of graphical setup.

for i = 1:4
nexttile(4*(i-1)+1);
xlabel('$X$','Interpreter','latex')
ylabel('$Y$','Interpreter','latex')
zlabel('$Z$','Interpreter','latex')
grid on
set(gca,'FontSize',18)
set(gca,'TickLabelInterpreter','latex')
view(30,10)
nexttile(4*(i-1)+2);
xlabel('$t$','Interpreter','latex')
ylabel('$X$','Interpreter','latex')
grid on
set(gca,'FontSize',18)
set(gca,'TickLabelInterpreter','latex')
nexttile(4*(i-1)+3);
xlabel('$t$','Interpreter','latex')
ylabel('$Y$','Interpreter','latex')
grid on
set(gca,'FontSize',18)
set(gca,'TickLabelInterpreter','latex')
nexttile(4*(i-1)+4);
xlabel('$t$','Interpreter','latex')
ylabel('$Z$','Interpreter','latex')
grid on
set(gca,'FontSize',18)
set(gca,'TickLabelInterpreter','latex')
end

c_beta.Ticks = [0 1];
c_beta.TickLabels = {'0',num2str(beta_var(end))};
c_beta.TickLabelInterpreter = 'latex';
c_beta.Label.String='$| \hat{\beta} |$';
c_beta.Label.Interpreter='latex';
c_beta.Label.Rotation=0;
c_beta.Label.Position=[2.5  0.55];
c_beta.Label.FontSize=20;

c_delta.Ticks = [0 1];
c_delta.TickLabels = {'0',num2str(delta_var(end))};
c_delta.TickLabelInterpreter = 'latex';
c_delta.Label.String='$| \hat{\delta} |$';
c_delta.Label.Interpreter='latex';
c_delta.Label.FontSize=20;
c_delta.Label.Rotation=0;
c_delta.Label.Position=[2.5  0.55];

c_gamma.Ticks = [0 1];
c_gamma.TickLabels = {'0',num2str(gamma_var(end))};
c_gamma.TickLabelInterpreter = 'latex';
c_gamma.Label.String='$| \hat{\gamma} |$';
c_gamma.Label.Interpreter='latex';
c_gamma.Label.FontSize=20;
c_gamma.Label.Rotation=0;
c_gamma.Label.Position=[2.5  0.55];

c_C.Ticks = [0 1];
c_C.TickLabels = {'0',num2str(C_var(end))};
c_C.TickLabelInterpreter = 'latex';
c_C.Label.String='$| \hat{C} |$';
c_C.Label.Interpreter='latex';
c_C.Label.FontSize=20;
c_C.Label.Rotation=0;
c_C.Label.Position=[2.5  0.58];

tl.Padding = 'none';
tl.TileSpacing = 'compact';

% exportgraphics(gcf,'beta_gamma_delta.eps','ContentType','vector')


%% Auxiliary functions 

function [theta,phi,psi,x,y,z] = traj_full(init_full,G,B,C,D,beta,gamma,delta,V1,V2,V3,tps)

% pack the parameters.
params = struct();
params.G = G;
params.B = B;
params.C = C;
params.D = D;
params.beta = beta;
params.gamma = gamma;
params.delta = delta;
params.V = [V1; V2; V3];

% solve the full dynamics.
options = odeset('RelTol', 1e-11, 'AbsTol', 1e-11);
[~, sol_full] = ode15s(@(t,state) ode_full(t,state,params),tps,init_full,options);
theta = sol_full(:,1);
phi = sol_full(:,2);
psi = sol_full(:,3);
x= sol_full(:,4);
y= sol_full(:,5);
z= sol_full(:,6);

end


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
    D = params.D;
    beta = params.beta;
    gamma = params.gamma;
    delta = params.delta;
    G = params.G;
    V1 = params.V(1);
    V2 = params.V(2);
    V3 = params.V(3);

  % Angular dynamics.
    f1 = -B * (sin(2*theta) * sin(2*phi)) / 4 - C*(sin(theta)*cos(2*phi))/2;
    f2 = (1 - B * cos(2*phi)) / 2 + (C/2)*sin(2*phi)*cos(theta);
    f3 = (B/2) * cos(theta) * cos(2*phi) - (C/2)*(cos(theta)).^2 .* sin(2*phi) - (D/2)*sin(theta).^2*sin(2*phi);

    d_state(1) = G*f1;
    d_state(2) = G*f2;
    d_state(3) = G*f3;

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

    g1 = -(delta-gamma)*sin(2*phi)*sin(theta)*sin(2*theta)/4 +...
         beta*sin(theta)^2*cos(2*phi)/2;

    g2 = -(delta-gamma)*sin(2*phi)*sin(theta)*sin(phi)*sin(theta)^2/2 -...
         gamma*cos(phi)*sin(theta)/2 + ...
         beta*sin(2*theta)*sin(phi)/4;

    g3 = (delta-gamma)*sin(2*phi)*sin(theta)*cos(phi)*sin(theta)^2/2 +...
         gamma*sin(phi)*sin(theta)/2 + ...
         beta*sin(2*theta)*cos(phi)/4;

    d_state(4) = d_state(4) + G*g1;
    d_state(5) = d_state(5) + G*g2;
    d_state(6) = d_state(6) + G*g3;
    
end




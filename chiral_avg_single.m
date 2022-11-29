%%%%%%%%%% Figure for wobbling swimmer with translational dynamics %%%%%%%%%
% (with chirality)

% Study of the averaged system.

% close all;

%% Setup.
% G = gamma, the shear rate of the flow.
G = 1;
% B = Bretherton constant.
B = 0;
% Ishimoto constant.
C = 1.5;

% Quantities used in the asymptotic analysis.
w = 1e-2; % spinning ratio.

lambda = sqrt(1 + w^2); 
% Effective Bretherton constant.
B_eff = B * (2 - w^2) / (2 * (1 + w^2));
% Effective chirality parameters.
C_eff_3 = C./lambda^3;

% ODE options.
options = odeset('RelTol', 1e-11, 'AbsTol', 1e-11);

% Time interval
T_min = 0;
T_max = 1e3;
tps = linspace(T_min,T_max,1e5);

params = struct();
params.G = G;
params.B = B;
params.C = C;
params.W_perp = W_perp;
params.W_par = W_par;
params.w = w;
params.lambda = lambda;

params.B_eff = B_eff;
params.C_eff3 = C_eff_3;

%% plot the sphere

figure(9); clf
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

%% Solve the reduced ODE.

            
init_theta = 0.8;
init_phi = -1.6;
init_psi = 1;
            
            [alpha,mu,varphi] = Initial_conditions_2(init_theta,init_phi,init_psi,w,1);
 
            % The IC vector for the reduced simulations.
            init_reduced = [alpha; varphi; mu];
            
            tic
            [~, sol_reduced] = ode15s(@(t,state) ode_reduced(t,state,params),tps,init_reduced,options);
            toc
            
            alpha_bar = sol_reduced(:,1);
            phi_bar = sol_reduced(:,2);
            mu_bar = sol_reduced(:,3);
            

            plot3(sin(alpha_bar).*cos(phi_bar),sin(alpha_bar).*sin(phi_bar),cos(alpha_bar),'b','LineWidth',1)
            plot3(sin(alpha_bar(1)).*cos(phi_bar(1)),sin(alpha_bar(1)).*sin(phi_bar(1)),cos(alpha_bar(1)),'ok','LineWidth',3)
            
            axis equal off
            
            drawnow
            

%% Auxiliary functions 

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
    C = params.C;
    C_eff3 = params.C_eff3;
    G = params.G;
    w = params.w;
    lambda = params.lambda;
    a = 1/w; %old notation for use inside I_hat and J_hat. also in those w is equal to sqrt(1+a^2). may be nice to change
   
    % Auxiliary angular dynamics.
    f1 = -B_eff * (sin(2*alpha_bar) * sin(2*phi_bar)) / 4 - C_eff3*(sin(alpha_bar)*cos(2*phi_bar))/2;
    f2 = (1 - B_eff * cos(2*phi_bar)) / 2 + (C_eff3/2)*sin(2*phi_bar)*cos(alpha_bar);
    f3 = (B_eff/2) * cos(alpha_bar) * cos(2*phi_bar) - (C_eff3/2)*(cos(alpha_bar)).^2 .* sin(2*phi_bar);

    % Extra terms.
    f3_extra = (- 3*w^2./2./lambda^3)*sin(2*alpha_bar)^2;
        
    d_state(1) = G*(f1);
    d_state(2) = G*(f2);
    d_state(3) = G*(f3 + (C/2)*sin(2*phi_bar)*f3_extra);

end

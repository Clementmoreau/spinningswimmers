%%%%%%%%%% Figure for wobbling swimmer with translational dynamics %%%%%%%%%
% (with chirality)

% Generated the data for translational, chiral movie.

close all;

% Data file
data_movie_trans_chiral = cell(5,4);
% data structure:
% W_par, W_perp
% B,C,D,beta,gamma,delta
% sol_full
% sol_reduced

% --- Parameters unchanged through the five simulations
% G = gamma, the shear rate of the flow.
G = 1;
% B = Bretherton parameter.
B = 0.7;
% Swimming velocity along axis of helicoidal symmetry, e_hat_1.
V1 = 1;
% Swimming velocity along axis e_hat_2.
V2 = 0;
% Swimming velocity along axis e_hat_3.
V3 = 0;

% ODE options.
options = odeset('RelTol', 1e-11, 'AbsTol', 1e-11);

% Generate specific IC for theta, phi, and psi.
init_theta = pi/4;
init_phi = pi/2;
init_psi = pi/4;

% Initial condition for swimmer position.
X0 = [0;0;0];

% Time interval
T_min = 0;
T_max = 50;
tps = linspace(T_min,T_max,1e5);

% --- Variable parameters
% Ishimoto parameter.
C_var = [0,0,-0.5,-0.5,-0.5];
% Extra chiral parameter.
D_var = [0,0,1,1,1];
% Extra constants to account for translational chiral drift.
beta_var = [0,0,0.1,0.1,0.1];
gamma_var = [0,0.5,0.5,0.5,0.5];
delta_var = [0,1,1,1,1];
% W_par is the intrinsic spin of the swimmer about its axis of helicoidal symmetry.
W_par_var = [1e-5,1e-5,1e-5,20,20];
% W_perp is the other direction of spin (perpendicular to the axis of helicoidal symmetry).
W_perp_var = [1e-5,1e-5,1e-5,4,20];

% main loop

for i = 1:length(C_var)
    C = C_var(i);
    D = D_var(i);
    beta = beta_var(i);
    gamma = gamma_var(i);
    delta = delta_var(i);
    W_par = W_par_var(i);
    W_perp = W_perp_var(i);

% Quantities used in the asymptotic analysis.
w = W_perp / W_par; % spinning ratio.
lambda = sqrt(1 + w^2); 

% Effective Bretherton constant.
B_eff = B * (2 - w^2) / (2 * (1 + w^2));
% Effective chirality parameters.
C_eff = (C+w^2*D)./lambda^3;
D_eff = (3*w^2*C+(2-w^2)*D)/2/lambda^3;
% Effective speed.
V_hat = (V1 + w*V2) / lambda;
% Effective translational parameters.
beta_eff = beta * (2 - w^2) / (2 * (1 + w^2));
gamma_eff = (gamma+w^2*delta)./lambda^3;
delta_eff = (3*w^2*gamma+(2-w^2)*delta)/2/lambda^3;

% Generate IC for the long-time variables, alpha_bar, mu_bar, and phi_bar.
[init_alpha_bar, init_mu_bar, init_phi_bar] = Initial_conditions_2(init_theta, init_phi, init_psi, W_perp, W_par);

% The IC vector for the full simulations.
init_full = [init_theta; init_phi; init_psi; X0];

% The IC vector for the reduced simulations.
init_reduced = [init_alpha_bar; init_phi_bar; init_mu_bar; X0];

params = struct();
params.G = G;
params.B = B;
params.C = C;
params.D = D;
params.beta = beta;
params.gamma = gamma;
params.delta = delta;
params.W_perp = W_perp;
params.W_par = W_par;
params.w = w;
params.lambda = lambda;
params.V = [V1; V2; V3];
params.V_hat = V_hat;

params.B_eff = B_eff;
params.C_eff = C_eff;
params.D_eff = D_eff;
params.beta_eff = beta_eff;
params.gamma_eff = gamma_eff;
params.delta_eff = delta_eff;

% Solve full ODE.
tic
[~, sol_full] = ode15s(@(t,state) ode_full(t,state,params),tps,init_full,options);
toc

% Solve the reduced ODE.
tic
[~, sol_reduced] = ode15s(@(t,state) ode_reduced(t,state,params),tps,init_reduced,options);
toc

% Save data

data_movie_trans_chiral{i,1} = [W_par W_perp];
data_movie_trans_chiral{i,2} = [B C D beta gamma delta];
data_movie_trans_chiral{i,3} = sol_full;
data_movie_trans_chiral{i,4} = sol_reduced;

end

save('data_movie_trans_chiral.mat','data_movie_trans_chiral')


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
    D = params.D;
    beta = params.beta;
    gamma = params.gamma;
    delta = params.delta;
    G = params.G;
    w1 = params.W_perp;
    w2 = params.W_par;
    V1 = params.V(1);
    V2 = params.V(2);
    V3 = params.V(3);

    % Angular dynamics.
    f1 = -B * (sin(2*theta) * sin(2*phi)) / 4 - C*(sin(theta)*cos(2*phi))/2;
    f2 = (1 - B * cos(2*phi)) / 2 + (C/2)*sin(2*phi)*cos(theta);
    f3 = (B/2) * cos(theta) * cos(2*phi) - (C/2)*(cos(theta)).^2 .* sin(2*phi) - (D/2)*sin(theta).^2*sin(2*phi);

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
    C_eff = params.C_eff;
    D_eff = params.D_eff;
    beta_eff = params.beta_eff;
    gamma_eff = params.gamma_eff;
    delta_eff = params.delta_eff;
    G = params.G;
    V_hat = params.V_hat;

    % Auxiliary angular dynamics.
    f1 = -B_eff * (sin(2*alpha_bar) * sin(2*phi_bar)) / 4 - C_eff*(sin(alpha_bar)*cos(2*phi_bar))/2;
    f2 = (1 - B_eff * cos(2*phi_bar)) / 2 + C_eff*sin(2*phi_bar)*cos(alpha_bar)/2;
    f3 = (B_eff/2) * cos(alpha_bar) * cos(2*phi_bar) - C_eff*(cos(alpha_bar)).^2 .* sin(2*phi_bar)/2 - (D_eff/2)*sin(alpha_bar).^2*sin(2*phi_bar);
        
    d_state(1) = G*(f1);
    d_state(2) = G*(f2);
    d_state(3) = G*(f3);

    % Average translational dynamics.
    d_state(4) = V_hat * cos(alpha_bar);
    d_state(5) = V_hat * sin(phi_bar) * sin(alpha_bar);
    d_state(6) = -V_hat * cos(phi_bar) * sin(alpha_bar) + G*y_bar;

    % Extra drift terms.
    g1 = -(delta_eff-gamma_eff)*sin(2*phi_bar)*sin(alpha_bar)*sin(2*alpha_bar)/4 +...
         beta_eff*sin(alpha_bar)^2*cos(2*phi_bar)/2;

    g2 = -(delta_eff-gamma_eff)*sin(2*phi_bar)*sin(alpha_bar)*sin(phi_bar)*sin(alpha_bar)^2/2 -...
         gamma_eff*cos(phi_bar)*sin(alpha_bar)/2 + ...
         beta_eff*sin(2*alpha_bar)*sin(phi_bar)/4;

    g3 = (delta_eff-gamma_eff)*sin(2*phi_bar)*sin(alpha_bar)*cos(phi_bar)*sin(alpha_bar)^2/2 +...
         gamma_eff*sin(phi_bar)*sin(alpha_bar)/2 + ...
         beta_eff*sin(2*alpha_bar)*cos(phi_bar)/4;

    % Extra drift terms.

    d_state(4) = d_state(4) + G*g1;
    d_state(5) = d_state(5) + G*g2;
    d_state(6) = d_state(6) + G*g3;
end
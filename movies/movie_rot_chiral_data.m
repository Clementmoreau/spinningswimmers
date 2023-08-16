%%%%%%%%%% Figure for wobbling swimmer with translational dynamics %%%%%%%%%
% (with chirality)

% Generates the data for the movies.

close all;



%% Setup.
% G = gamma, the shear rate of the flow.
G = 1;

% B = Bretherton parameter.
B_var = [0.7,0.7,0,0];

% Ishimoto parameter.
C_var = [0.7,1.5,0.7,1.5];

% Extra chiral parameter (mostly irrelevant here).
D = 0;

% W_par is the intrinsic spin of the swimmer about its axis of helicoidal symmetry.
W_par = 15;
% W_perp is the other direction of spin (perpendicular to the axis of helicoidal symmetry).
W_perp_var = W_par * [1e-8,0.5,1,2];

data_movie_rot_chiral = cell(length(B_var),length(W_perp_var),9);

for i = 1:length(B_var)
    for j = 1:length(W_perp_var)

        B = B_var(i);
        C = C_var(i);
        W_perp = W_perp_var(j);

        % Quantities used in the asymptotic analysis.
        w = W_perp / W_par; % spinning ratio.
        lambda = sqrt(1 + w^2);

        % Effective Bretherton constant.
        B_eff = B * (2 - w^2) / (2 * (1 + w^2));
        % Effective chirality parameters.
        C_eff = (C+w^2*D)./lambda^3;
        D_eff = (3*w^2*C+(2-w^2)*D)/2/lambda^3;

        % ODE options.
        options = odeset('RelTol', 1e-11, 'AbsTol', 1e-11);

        % Generate specific IC for theta, phi, and psi.

        init_theta = pi/2;
        init_phi = -pi-eps;
        init_psi = 0;

        % Generate IC for the long-time variables, alpha_bar, mu_bar, and phi_bar.
        [init_alpha_bar, init_mu_bar, init_phi_bar] = Initial_conditions_2(init_theta, init_phi, init_psi, W_perp, W_par);

        % The IC vector for the full simulations.
        init_full = [init_theta; init_phi; init_psi];

        % The IC vector for the reduced simulations.
        init_reduced = [init_alpha_bar; init_phi_bar; init_mu_bar];

        % Time interval
        T_min = 0;
        T_max = 30;
        tps = linspace(T_min,T_max,7200);

        params = struct();
        params.G = G;
        params.B = B;
        params.C = C;
        params.D = D;
        params.W_perp = W_perp;
        params.W_par = W_par;
        params.w = w;
        params.lambda = lambda;

        params.B_eff = B_eff;
        params.C_eff = C_eff;
        params.D_eff = D_eff;

        % Solve full ODE.
        tic
        [~, sol_full] = ode15s(@(t,state) ode_full(t,state,params),tps,init_full,options);
        toc

        theta = sol_full(:,1);
        phi = sol_full(:,2);
        psi = sol_full(:,3);

        % Solve the reduced ODE.
        tic
        [~, sol_reduced] = ode15s(@(t,state) ode_reduced(t,state,params),tps,init_reduced,options);
        toc

        alpha = sol_reduced(:,1);
        phir = sol_reduced(:,2);
        mu = sol_reduced(:,3);

        % Save the data.

        % the data contains:
        % w_par
        % w_perp
        % [B,C,B_eff,C_eff];
        % theta_full
        % phi_full
        % psi_full
        % alpha_bar
        % phi_bar
        % mu_bar
        data_movie_rot_chiral{i,j,1} = W_par;
        data_movie_rot_chiral{i,j,2} = W_perp;
        data_movie_rot_chiral{i,j,3} = [B,C,B_eff,C_eff];
        data_movie_rot_chiral{i,j,4} = theta;
        data_movie_rot_chiral{i,j,5} = phi;
        data_movie_rot_chiral{i,j,6} = psi;
        data_movie_rot_chiral{i,j,7} = alpha;
        data_movie_rot_chiral{i,j,8} = phir;
        data_movie_rot_chiral{i,j,9} = mu;

    end
end

save('data_movie_rot_chiral.mat')

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
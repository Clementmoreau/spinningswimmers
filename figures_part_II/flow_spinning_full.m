%%%%% Figure showing the rotational dynamics of the spinning swimmer in the
% (phi,theta) plane, for different values of B, C and omega.

% for each row, we set C,B,phi_init.
FlowDat=[ 0.0, 0.7, -1;...
    0.7, 0.7, -1;...
    1.6, 0.7, -1;...
    0.7, 0, -1;...
    1.6, 0, -1;...
    ];

% foe each column, a value of omega.
w_var = [-1e-8 0.5 1 2];

% setup the figure.
figure(5);clf;
set(gcf, 'Position',  [100, 100, 760, 1000])
tl = tiledlayout(size(FlowDat,1),length(w_var));


% loop over both (B,C) and omega.
for label=1:size(FlowDat,1)
    C =-FlowDat(label,1);
    B =-FlowDat(label,2);
    Yinit =FlowDat(label,3)*pi;
    for i = 1:length(w_var)
        nexttile((label-1)*4+i);
        fig_flow_draw(label,w_var(i),C,B,Yinit);
    end
end

% graphical setup.
tl.Padding = 'none';

% export the figure.
exportgraphics(gcf,'fig_flow_spinning.eps','ContentType','vector')
exportgraphics(gcf,'fig_flow_spinning.png','Resolution',300)

function out = fig_flow_draw(label,w,C,B,Yinit)

% resolution for the streamlines.
Ntheta=400;
Nphi=800;
% resolution for the vector field.
NQ=40;
% some graphical parameters.
col = [0.2,0.2,0.2]; %vector field color
str_width = 1; %streamline width
linejoin = 'chamfer'; %line drawing
str_precision = 1.5; %streamline precision

% effective values of the parameters.
D = 0;
B_hat = (2-w^2)/2/(1+w^2)*B;
C_hat = 1/(1+w^2)^(3/2)*C;
D_hat = (3*w^2*C+(2-w^2)*D)/2/(1+w^2)^(3/2);

% first plot the aspect of one full trajectory

% record the parameters.
params.B = B;
params.C = C;
params.D = 0;
params.B_hat = B_hat;
params.C_hat = C_hat;
params.D_hat = D_hat;
params.G = 1;
params.w2 = 10;
params.w1 = w*params.w2;

% initial condition.
init_full = [pi/2,-pi-eps,0];

% time interval.
T_min = 0;
T_max = 7 + 5*(w==0.5)*(label==3);
tps = linspace(T_min,T_max,1e3);

% solve the ode.
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
[~, sol_full] = ode15s(@(t,state) ode_full(t,state,params),tps,init_full,options);

% Generate IC for the long-time variables, alpha_bar, mu_bar, and phi_bar.
[init_alpha_bar, init_mu_bar, init_phi_bar] = Initial_conditions_2(init_full(1), init_full(2), init_full(3), 10*w, 10);
init_reduced = [init_alpha_bar; init_phi_bar; init_mu_bar];

% solve the reduced ode.
[~, sol_reduced] = ode15s(@(t,state) ode_reduced(t,state,params),tps,init_reduced,options);

% plot the solutions.

if (w==0.5)*(label==3)==1
    % this is to avoid a bad plotting in a single case when phi crosses the
    % -pi,pi interval, manually removing the three lines we dont want
    I1 =1:359;I2=360:407;I3=408:455;I4=456:1000;
    plot((mod(sol_full(I1,2)+pi,2*pi)-pi)/pi,sol_full(I1,1)/pi,'Color',[0.4,0.4,1],'LineWidth',1)
    hold on
    plot((mod(sol_full(I2,2)+pi,2*pi)-pi)/pi,sol_full(I2,1)/pi,'Color',[0.4,0.4,1],'LineWidth',1)
    plot((mod(sol_full(I3,2)+pi,2*pi)-pi)/pi,sol_full(I3,1)/pi,'Color',[0.4,0.4,1],'LineWidth',1)
    plot((mod(sol_full(I4,2)+pi,2*pi)-pi)/pi,sol_full(I4,1)/pi,'Color',[0.4,0.4,1],'LineWidth',1)
    drawnow
else
    plot((mod(sol_full(:,2)+pi,2*pi)-pi)/pi,sol_full(:,1)/pi,'Color',[0.4,0.4,1],'LineWidth',1)
    hold on
    drawnow
end

% prepare the vector field and streamlines.
delta_theta=pi/Ntheta;
delta_phi=2*pi/Nphi;
[phi,theta] = meshgrid(-2*pi:delta_phi:2*pi,0:delta_theta:pi);

U=0.5*(1-B_hat*cos(2*phi)) + C_hat*cos(theta).*cos(phi).*sin(phi);
V= -0.25*B_hat*sin(2*theta).*sin(2*phi)-0.5*C_hat*sin(theta).*cos(2*phi);

[phiQ,thetaQ] = meshgrid(-2*pi:delta_phi*NQ:2*pi,0:delta_theta*NQ:pi);

UQ=0.5*(1-B_hat*cos(2*phiQ)) + C_hat*cos(thetaQ).*cos(phiQ).*sin(phiQ);
VQ= -0.25*B_hat*sin(2*thetaQ).*sin(2*phiQ)-0.5*C_hat*sin(thetaQ).*cos(2*phiQ);

% starting positions of the streamlines.
starty = 0.1:0.1:0.9;

startx1 = (Yinit/pi)*ones(size(starty));
starty1= starty;

startx2 = (Yinit/pi)*ones(size(starty));
starty2= starty;

hold on
% plot the vector field.
quiver(phiQ/pi,thetaQ/pi,UQ/pi,VQ/pi,'LineWidth',1,'Color',[0.6 0.6 0.6]);

% plot the streamlines.
str1=streamline(phi/pi,theta/pi, U/pi, V/pi,startx1,starty1,[str_precision 2000]);
str2=streamline(phi/pi,theta/pi, -U/pi, -V/pi,startx2,starty2,[str_precision 2000]);

set(str1,'LineWidth',str_width,'Color',col,'LineJoin',linejoin)
set(str2,'LineWidth',str_width,'Color',col,'LineJoin',linejoin)

% extra streamlines manually adjusted in asymmetrical plots
if label==2
    if w < 0.1
        str=streamline(phi/pi,theta/pi, -U/pi, -V/pi,1*ones(1,9),0.1:0.1:0.9,str_precision);
        set(str,'LineWidth',str_width,'Color',col,'LineJoin',linejoin)
    end
end
if label==3
    if w < 0.1
        str=streamline(phi/pi,theta/pi, -U/pi, -V/pi,1*ones(1,9),0.1:0.1:0.9,[str_precision 2000]);
        set(str,'LineWidth',str_width,'Color',col,'LineJoin',linejoin)
    end
    if w ==0.5
        str=streamline(phi/pi,theta/pi, -U/pi, -V/pi,1*ones(1,4),[0.5 0.7 0.8 0.9],[str_precision 2000]);
        set(str,'LineWidth',str_width,'Color',col,'LineJoin',linejoin)
    end
end
if label==5
    if w < 0.1
        startx = [-0.25*ones(1,2),0.75*ones(1,2),0.25*ones(1,2),-0.75*ones(1,2)];
        starty = [0.9 0.8 0.9 0.8 0.1 0.2 0.1 0.2];
        str=streamline(phi/pi,theta/pi, U/pi, V/pi,startx,starty,[0.16 2000]);
        set(str,'LineWidth',str_width,'Color',col,'LineJoin',linejoin)
    end
    if w ==0.5
        startx = [-0.25,0.75,0.25,-0.75];
        starty = [0.95 0.95 0.05 0.05];
        str=streamline(phi/pi,theta/pi, U/pi, V/pi,startx,starty,[0.1 2000]);
        set(str,'LineWidth',str_width,'Color',col,'LineJoin',linejoin)
    end
end

% and finally, plot the averaged dynamics.
plot((mod(sol_reduced(3:end,2)+pi,2*pi)-pi)/pi,sol_reduced(3:end,1)/pi,'Color',[0.8,0.1,0.1],'LineWidth',3)

% graphical setup: axis labels, etc.
if label ==5
    xlabel('$\varphi$','Interpreter','latex')
end
if w < 0.1
    ylabel('$\theta$','Interpreter','latex')
end

set(gca,'FontSize',16)

if label==1
    if w < 0.1
    title(['$\omega = 0$'],'Interpreter','latex','FontSize',24)
    else
    title(['$\omega = ',num2str(w),'$'],'Interpreter','latex','FontSize',24)
    end
end

xlim([-1 1])
xticks([-1 1])
xticklabels({'$-\pi$','$\pi$'})
ylim([0 1])
yticks([0 1])
yticklabels({'$0$','$\pi$'})

set(gca,'TickLabelInterpreter','latex')

grid on
box on

end

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
w1 = params.w1;
w2 = params.w2;

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
    B_eff = params.B_hat;
    C_eff = params.C_hat;
    D_eff = params.D_hat;
    G = params.G;

    % Auxiliary angular dynamics.
    f1 = -B_eff * (sin(2*alpha_bar) * sin(2*phi_bar)) / 4 - C_eff*(sin(alpha_bar)*cos(2*phi_bar))/2;
    f2 = (1 - B_eff * cos(2*phi_bar)) / 2 + C_eff*sin(2*phi_bar)*cos(alpha_bar)/2;
    f3 = (B_eff/2) * cos(alpha_bar) * cos(2*phi_bar) - C_eff*(cos(alpha_bar)).^2 .* sin(2*phi_bar)/2 - (D_eff/2)*sin(alpha_bar).^2*sin(2*phi_bar);
        
    d_state(1) = G*(f1);
    d_state(2) = G*(f2);
    d_state(3) = G*(f3);

end



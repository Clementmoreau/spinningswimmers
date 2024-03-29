%%%%%%%%%% Movie for wobbling swimmer with translational dynamics %%%%%%%%%
% (without chirality)

% Movie associated to the figure for translational dynamics, bacterial limit

%Load the data
load('data_movie_BL.mat');

%Prepare the figure.
f=figure(2);clf
f.Position=[0 0 1000 1100];
tiledlayout(2,3)
ribbonwidth = 0.7;
v = [60,25];

% Do we want to save the movie ?
savemov = 0;

if savemov
    v1 = VideoWriter('movie_BL.mp4','MPEG-4');
    v1.Quality = 50;
    open(v1);
end

% Interval between frames.
tick = 200;

a = 1e-1*(max(z_full_6)-min(z_full_6)); %swimmer size

% Axis.
xmin = min(x_full_6)-1.2*a;
xmax = max(x_full_6)+1.2*a;
ymin = min(y_full_6)-1.2*a;
ymax = max(y_full_6)+1.2*a;
zmin = min(z_full_6)-1.2*a;
zmax = max(z_full_6)+1.2*a;

colormap(parula(256))

% Aspect ratio.
r = sqrt((1-B)/(1+B));

for i = 1:tick:length(tps)

    nexttile(1);
    % plot the full trajectory.
    plot3(x_full_1(1:i),y_full_1(1:i),z_full_1(1:i),'k','LineWidth',0.5)
    hold on
    plot_swimmer(r,a,x_full_1(i),y_full_1(i),z_full_1(i),phi_full_1(i),theta_full_1(i),psi_full_1(i))
    hold off
    % some graphical parameters
    axis equal
    set(gcf,'color','w');
    axis([xmin xmax ymin ymax zmin zmax])
    view(60,25)
    xlabel('x','interpreter','latex')
    ylabel('y','interpreter','latex')
    zlabel('z','interpreter','latex')
    grid on
    set(gca,'FontSize',20)
    set(gca,'TickLabelInterpreter','latex')

    nexttile(2);
    % Plot the average trajectory.
    plot3(x_bar_2(1:i),y_bar_2(1:i),z_bar_2(1:i),'r','LineWidth',1.5)
    hold on
    % plot the full trajectory.
    plot3(x_full_2(1:i),y_full_2(1:i),z_full_2(1:i),'k','LineWidth',0.5)
    plot_swimmer(r,a,x_full_2(i),y_full_2(i),z_full_2(i),phi_full_2(i),theta_full_2(i),psi_full_2(i))
    hold off
    % some graphical parameters
    axis equal
    set(gcf,'color','w');
    axis([xmin xmax ymin ymax zmin zmax])
    view(60,25)
    xlabel('x','interpreter','latex')
    ylabel('y','interpreter','latex')
    zlabel('z','interpreter','latex')
    grid on
    set(gca,'FontSize',20)
    set(gca,'TickLabelInterpreter','latex')

    nexttile(3);
    % Plot the average trajectory.
    plot3(x_bar_3(1:i),y_bar_3(1:i),z_bar_3(1:i),'r','LineWidth',1.5)
    hold on
    % plot the full trajectory.
    plot3(x_full_3(1:i),y_full_3(1:i),z_full_3(1:i),'k','LineWidth',0.5)
    plot_swimmer(r,a,x_full_3(i),y_full_3(i),z_full_3(i),phi_full_3(i),theta_full_3(i),psi_full_3(i))
    hold off
    % some graphical parameters
    axis equal
    set(gcf,'color','w');
    axis([xmin xmax ymin ymax zmin zmax])
    view(60,25)
    xlabel('x','interpreter','latex')
    ylabel('y','interpreter','latex')
    zlabel('z','interpreter','latex')
    grid on
    set(gca,'FontSize',20)
    set(gca,'TickLabelInterpreter','latex')

    nexttile(4);
    % Plot the average trajectory.
    plot3(x_bar_4(1:i),y_bar_4(1:i),z_bar_4(1:i),'r','LineWidth',1.5)
    hold on
    % plot the full trajectory.
    plot3(x_full_4(1:i),y_full_4(1:i),z_full_4(1:i),'k','LineWidth',0.5)
    plot_swimmer(r,a,x_full_4(i),y_full_4(i),z_full_4(i),phi_full_4(i),theta_full_4(i),psi_full_4(i))
    hold off
    % some graphical parameters
    axis equal
    set(gcf,'color','w');
    axis([xmin xmax ymin ymax zmin zmax])
    view(60,25)
    xlabel('x','interpreter','latex')
    ylabel('y','interpreter','latex')
    zlabel('z','interpreter','latex')
    grid on
    set(gca,'FontSize',20)
    set(gca,'TickLabelInterpreter','latex')

    nexttile(5);
    % Plot the average trajectory.
    plot3(x_bar_5(1:i),y_bar_5(1:i),z_bar_5(1:i),'r','LineWidth',1.5)
    hold on
    % plot the full trajectory.
    plot3(x_full_5(1:i),y_full_5(1:i),z_full_5(1:i),'k','LineWidth',0.5)
    plot_swimmer(r,a,x_full_5(i),y_full_5(i),z_full_5(i),phi_full_5(i),theta_full_5(i),psi_full_5(i))
    hold off
    % some graphical parameters
    axis equal
    set(gcf,'color','w');
    axis([xmin xmax ymin ymax zmin zmax])
    view(60,25)
    xlabel('x','interpreter','latex')
    ylabel('y','interpreter','latex')
    zlabel('z','interpreter','latex')
    grid on
    set(gca,'FontSize',20)
    set(gca,'TickLabelInterpreter','latex')

    nexttile(6);
    % Plot the average trajectory.
    plot3(x_bar_6(1:i),y_bar_6(1:i),z_bar_6(1:i),'r','LineWidth',1.5)
    hold on
    % plot the full trajectory.
    plot3(x_full_6(1:i),y_full_6(1:i),z_full_6(1:i),'k','LineWidth',0.5)
    plot_swimmer(r,a,x_full_6(i),y_full_6(i),z_full_6(i),phi_full_6(i),theta_full_6(i),psi_full_6(i))
    hold off
    % some graphical parameters
    axis equal
    set(gcf,'color','w');
    axis([xmin xmax ymin ymax zmin zmax])
    view(60,25)
    xlabel('x','interpreter','latex')
    ylabel('y','interpreter','latex')
    zlabel('z','interpreter','latex')
    grid on
    set(gca,'FontSize',20)
    set(gca,'TickLabelInterpreter','latex')

    if savemov
        frame = getframe(gcf);
        writeVideo(v1,frame);
    end

    drawnow

end

if savemov
        close(v1);
end

%% post-processing

nexttile(1);
title('(a)','Interpreter','latex')
nexttile(2);
title('(b)','Interpreter','latex')
nexttile(3);
title('(c)','Interpreter','latex')
nexttile(4);
title('(d)','Interpreter','latex')
nexttile(5);
title('(e)','Interpreter','latex')
nexttile(6);
title('(f)','Interpreter','latex')

%% auxiliary functions

function [] = plot_swimmer(aspect_ratio,swimmer_size,x,y,z,phi,theta,psi)
% Show the swimmer.
    r = aspect_ratio;
    a = swimmer_size;
    % First define a sphere.
    [sx0,sy0,sz0]=sphere(50);
    % Rescale to make it a spheroid.
    sx=min(r,1)*sx0;sy=min(r,1)*sy0;sz=min(1/r,1)*sz0;
    % Plot it and colour it so the spinning will be visible.
    s = surf(a*sx+x,a*sy+y,a*sz+z,3*sz0+cos(6*sx0)+cos(4*sy0),'LineStyle','none','FaceLighting','gouraud');
    % Finally, rotate it to the appropriate orientation.
    o = [x,y,z];
    % First, azimuthal angle phi.
    rotate(s,[0 0 1],rad2deg(phi),o)
    % Then, polar angle theta.
    rotate(s,[cos(phi),sin(phi),0],rad2deg(theta),o)
    % And finally proper rotation psi.
    rotate(s,[sin(theta)*sin(phi),-sin(theta)*cos(phi),cos(theta)],rad2deg(psi),o)
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
    axis([-1 4.1 -3 4 -12 3])
    colormap(parula(256));
    col = [col(end:-1:1,:);col];
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
    axis([-1 4.1 -3 4 -12 3])
    colormap(parula(256));
    col = [col(end:-1:1,:);col];
    colormap(col)
    % colorbar

    set(gca,'FontSize',20)
    set(gca,'TickLabelInterpreter','latex')
    %zticks([0 8])
end
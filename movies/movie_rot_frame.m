%% Prepare the figure.

%Load the data
load('data_movie_1.mat');

save_movie = 1;

IF = [133,9,9,7,3,2];
I_end_movie = [I_end,I_end/3,I_end/3,I_end/5,I_end/10,I_end/12];
RS = [1,1,0.9,0.8,0.65,0.4];

for i_movie = 1:6

%Prepare the figure.
f=figure(1);clf
set(gcf, 'Position',  [1, 640, 800, 800])
set(gcf,'color','w');
%tl = tiledlayout(2,2);
%v=[-105,40];

theta = data_movie_rot_achiral{i_movie,4};
phi = data_movie_rot_achiral{i_movie,5};
psi = data_movie_rot_achiral{i_movie,6};
alpha = data_movie_rot_achiral{i_movie,7};
phir = data_movie_rot_achiral{i_movie,8};
mu = data_movie_rot_achiral{i_movie,9};

I_end = length(theta);

if save_movie
    v1 = VideoWriter(['movie_rot_achiral',num2str(i_movie),'.mp4'],'MPEG-4');
    % v1.Quality = 50;
    open(v1);
end

% number of frames depending on the case
i_frame = IF(i_movie);

for i = 1:i_frame:I_end_movie(i_movie)
    hold off

    % Plot the swimmer
    [SX,SY,SZ]=sphere(60);
    % Deform the swimmer
    B = 0.99;
    r = sqrt((1+B)/(1-B));
    if r>1
        SX = SX/r;
        SY = SY/r;
    else
        SZ = SZ*r;
    end
    o = [0,0,0];
    % Rotate the swimmer
    s = surf(SX,SY,SZ,3*sin(pi/2*SZ)+2*cos(4*r*SX)+2*cos(4*r*SY),'LineStyle','none','FaceLighting','gouraud');
    rotate(s,[0 1 0],rad2deg(theta(i)),o)
    rotate(s,[0 0 1],rad2deg(phi(i)),o)
    rotate(s,[sin(theta(i))*cos(phi(i)),sin(theta(i))*sin(phi(i)),cos(theta(i))],rad2deg(psi(i)),o)

    colormap(parula(256));
   

    shading interp
    material dull
    hold on

    % Plot the equivalent swimmer (except in the first case)

    if i_movie > 1
        [SX,SY,SZ]=sphere(30);
        % Deform the swimmer
        Bhat = data_movie_rot_achiral{i_movie,3};
        r = sqrt((1+Bhat)/(1-Bhat));
        SX = SX/r;
        SY = SY/r;
        % resize to be more visually pleasing
        SX = SX*RS(i_movie);
        SY = SY*RS(i_movie);
        SZ = SZ*RS(i_movie);
        o = [0,0,0];
        % Rotate the swimmer
        s = surf(SX,SY,SZ,'FaceLighting','gouraud','FaceAlpha',0.35,'FaceColor',[1 1 1],'EdgeColor',[0.8 0.8 0.8],'EdgeAlpha',0.4);
        rotate(s,[0 1 0],rad2deg(alpha(i)),o)
        rotate(s,[0 0 1],rad2deg(phir(i)),o)
        rotate(s,[sin(alpha(i))*cos(phir(i)),sin(alpha(i))*sin(phir(i)),cos(alpha(i))],rad2deg(mu(i)),o)
    end

    hold on

    % Plot the axis and (xy) unit circle
    plot3([0 1.5],[0 0],[0 0],'k','LineWidth',1)
    plot3([0 0],[0 1.5],[0 0],'k','LineWidth',1)
    plot3([0 0],[0 0],[0 1.5],'k','LineWidth',1)
    tt=0:0.1:2*pi+0.1;
    plot3(cos(tt),sin(tt),zeros(1,length(tt)),'k','LineWidth',1)

    % Full traj so far.
    plot3(sin(theta(1:i)).*cos(phi(1:i)),sin(theta(1:i)).*sin(phi(1:i)),cos(theta(1:i)),'k')

    % Average traj so far (except in the first case).
    if i_movie > 1
        plot3(sin(alpha(1:i)).*cos(phir(1:i)),sin(alpha(1:i)).*sin(phir(1:i)),cos(alpha(1:i)),'Color',[0.7 0.7 0.7],'LineWidth',4)
        plot3(sin(alpha(1)).*cos(phir(1)),sin(alpha(1)).*sin(phir(1)),cos(alpha(1)),'k.','MarkerSize',30)
    end

    % Mark the initial conditions.
    plot3(sin(theta(1)).*cos(phi(1)),sin(theta(1)).*sin(phi(1)),cos(theta(1)),'k.','MarkerSize',30)


    view(120,25)
    camlight
    grid on
    box on

    lim = 1.15;
    
    axis off equal
    axis([-lim lim -lim lim -lim lim])

    % save the frame
    if save_movie
        frame = getframe(gcf);
        writeVideo(v1,frame);
    end

    drawnow
end

if save_movie
    close(v1);
end

end

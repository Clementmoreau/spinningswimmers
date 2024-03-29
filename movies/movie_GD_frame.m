%%%%%%%%%% Movie for wobbling swimmer with translational dynamics %%%%%%%%%
% (without chirality)

% Movie associated to the figure for translational dynamics, general
% dynamics

%% Prepare the figure.

%Load the data
load('data_movie_GD.mat');

%Prepare the figure.
f=figure(2);clf
set(gcf, 'Position',  [1, 640, 1100, 700])
tl = tiledlayout(2,2);
v=[-105,40];

% Do we want to save the movie ?
savemov = 1;

if savemov
    v1 = VideoWriter('movie_GD.mp4','MPEG-4');
    v1.Quality = 50;
    open(v1);
end

% Interval between frames.
tick = 10;

colormap(parula(256))

% Aspect ratio.
r = sqrt((1-B)/(1+B));

% Axis for each subplot.

xlim=zeros(2,length(V2_var));
ylim=zeros(2,length(V2_var));
zlim=zeros(2,length(V2_var));

for j=1:length(V2_var)
    % Look for the extreme points.
    for i = 1:length(W_par_var)
        mx = min(data_full{i,j,4});
        Mx = max(data_full{i,j,4});
        my = min(data_full{i,j,5});
        My = max(data_full{i,j,5});
        mz = min(data_full{i,j,6});
        Mz = max(data_full{i,j,6});

        if mx < xlim(1,j)
            xlim(1,j) = mx;
        end
        if Mx > xlim(2,j)
            xlim(2,j) = Mx;
        end
        if my < ylim(1,j)
            ylim(1,j) = my;
        end
        if My > ylim(2,j)
            ylim(2,j) = My;
        end
        if mz < zlim(1,j)
            zlim(1,j) = mz;
        end
        if Mz > zlim(2,j)
            zlim(2,j) = Mz;
        end

    end
end

% Small padding.

a = 5e-2*(max(data_full{1,1,6})-min(data_full{1,1,6})); %swimmer size

xlim(1,:) = xlim(1,:)-1.2*a;
xlim(2,:) = xlim(2,:)+1.2*a;
ylim(1,:) = ylim(1,:)-1.2*a;
ylim(2,:) = ylim(2,:)+1.2*a;
zlim(1,:) = zlim(1,:)-1.2*a;
zlim(2,:) = zlim(2,:)+1.2*a;


%% Frames.

titles = {'(a)','(b)','(c)','(d)'};

for f = 1:tick:length(tps)

    for j = 1:length(V2_var)

        nexttile(j);
        

        for i = 1:length(W_par_var)

            % Access the data.

            theta_full = data_full{i,j,1};
            phi_full = data_full{i,j,2};
            psi_full = data_full{i,j,3};
            x_full = data_full{i,j,4};
            y_full = data_full{i,j,5};
            z_full = data_full{i,j,6};

            x_bar = data_bar{i,j,4};
            y_bar = data_bar{i,j,5};
            z_bar = data_bar{i,j,6};

            % Plot the average trajectory.
            plot3(x_bar(1:f),z_bar(1:f),y_bar(1:f),'r','LineWidth',1.5)
            hold on

            % plot the full trajectory
            plot3(x_full(1:f),z_full(1:f),y_full(1:f),'k','LineWidth',0.5)

            % plot the swimmer
            plot_swimmer(r,a,x_full(f),z_full(f),y_full(f),phi_full(f),theta_full(f),psi_full(f))

        end

        % some graphical parameters
        axis equal
        set(gcf,'color','w');
        axis([xlim(1,j) xlim(2,j) zlim(1,j) zlim(2,j) ylim(1,j) ylim(2,j)])
        set(gca, 'YDir', 'reverse')
        view(v)
        xlabel('x','interpreter','latex')
        ylabel('z','interpreter','latex')
        zlabel('y','interpreter','latex')
        title(titles{j},'Interpreter','latex')
        grid on
        set(gca,'FontSize',20)
        set(gca,'TickLabelInterpreter','latex')

        hold off

    end

    if savemov
        frame = getframe(gcf);
        writeVideo(v1,frame);
    end

    drawnow

end

if savemov
        close(v1);
end

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



   
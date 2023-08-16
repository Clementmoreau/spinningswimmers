% Load the data
load('data_movie_trans_chiral.mat');

% Prepare the figure

B=0.8;

f=figure(1);clf
set(gcf, 'Position',  [1, 640, 2400, 800])
set(gcf,'color','w');
tl = tiledlayout(1,3);

% Set the axis limits for the left panel (through all trajectories)
xm = 0;xM = 0;ym=0;yM=0;zm=0;zM=0;
for i = 1:5
    x = data_movie_trans_chiral{i,3}(:,4);
    y = data_movie_trans_chiral{i,3}(:,5);
    z = data_movie_trans_chiral{i,3}(:,6);
    xm = min(xm,min(x));
    xM = max(xM,max(x));
    ym = min(ym,min(y));
    yM = max(yM,max(y));
    zm = min(zm,min(z));
    zM = max(zM,max(z));
end
pad = 3;
xm = xm-pad;
ym = ym-pad;
zm = zm-pad;
xM = xM+pad;
yM = yM+pad;
zM = zM+pad;

% Color scales (to be changed by the user)
blues = colormap(parula(256));
reds = colormap(parula(256));
greens = colormap(parula(256));
purples = colormap(parula(256));
oranges = colormap(parula(256));

cols = cell(5);
cols{1} = blues;cols{2} = reds;cols{3}=greens;cols{4}=purples;cols{5}=oranges;

save_movie = 1;

for i_movie = 1:5

    % Frame rate

    if save_movie
        v1 = VideoWriter(['movie_trans_chiral',num2str(i_movie),'.mp4'],'MPEG-4');
        v1.FrameRate = 50;
        v1.Quality = 50;
        open(v1);
    end

    if i_movie > 3
        i_frame = 25;
    else
        i_frame = 100;
    end

    % Extract the data
    x = data_movie_trans_chiral{i_movie,3}(:,4);
    y = data_movie_trans_chiral{i_movie,3}(:,5);
    z = data_movie_trans_chiral{i_movie,3}(:,6);

    theta = data_movie_trans_chiral{i_movie,3}(:,1);
    phi = data_movie_trans_chiral{i_movie,3}(:,2);
    psi = data_movie_trans_chiral{i_movie,3}(:,3);

    xr = data_movie_trans_chiral{i_movie,4}(:,4);
    yr = data_movie_trans_chiral{i_movie,4}(:,5);
    zr = data_movie_trans_chiral{i_movie,4}(:,6);

    thr = data_movie_trans_chiral{i_movie,4}(:,1);
    phr = data_movie_trans_chiral{i_movie,4}(:,2);
    psr = data_movie_trans_chiral{i_movie,4}(:,3);

    i_light = 50;

    % Generate each frame.

    for i = 1:i_frame:length(x)

        nexttile(1);

        hold off

        % Plot each trajectory in light colour
        for j = 1:5
            i_traj = 3 + 1*(j>3);
            xa = data_movie_trans_chiral{j,i_traj}(:,4);
            ya = data_movie_trans_chiral{j,i_traj}(:,5);
            za = data_movie_trans_chiral{j,i_traj}(:,6);
            plot3(xa,ya,za,'LineWidth',2,'Color',cols{j}(i_light,:));
            hold on
        end

        if i_movie > 3

            % Full traj so far.
            plot3(x(1:i),y(1:i),z(1:i),'LineWidth',0.5,'Color','k')
            hold on
            % Avg traj so far.
            plot3(xr(1:i),yr(1:i),zr(1:i),'LineWidth',4,'Color',cols{i_movie}(end,:))

        else
            % Full traj so far.
            plot3(x(1:i),y(1:i),z(1:i),'LineWidth',4,'Color',cols{i_movie}(end,:))

        end

        % vector e1hat
        quiver3([x(i)],[y(i)],[z(i)],cos(theta(i)),sin(phi(i))*sin(theta(i)),-cos(phi(i))*sin(theta(i)),2,'Color','k','LineWidth',3)

        % vector e2hat
        quiver3([x(i)],[y(i)],[z(i)],sin(psi(i))*sin(theta(i)),...
            ( cos(phi(i))*cos(psi(i)) - cos(theta(i))*sin(phi(i))*sin(psi(i))),...
            ( sin(phi(i))*cos(psi(i)) + cos(theta(i))*cos(phi(i))*sin(psi(i))),2,'Color','k','LineWidth',3)

        % vector e3hat
        quiver3([x(i)],[y(i)],[z(i)],cos(psi(i))*sin(theta(i)),...
            ( -cos(phi(i))*sin(psi(i)) - cos(theta(i))*sin(phi(i))*cos(psi(i))),...
            ( - sin(phi(i))*sin(psi(i)) + cos(theta(i))*cos(phi(i))*cos(psi(i))),2,'Color','k','LineWidth',3)

        axis equal
        axis([xm xM ym yM zm zM])
        grid on
        box on
        hold off

        xlabel('$x$','Interpreter','latex')
        ylabel('$y$','Interpreter','latex')
        zlabel('$z$','Interpreter','latex')
        set(gca,'FontSize',24)
        set(gca,'TickLabelInterpreter','latex')

        nexttile(2);

        hold off

        % Plot the swimmer
        n_sphere = 100;
        [SX,SY,SZ]=sphere(n_sphere);
        % aspect ratio
        r = sqrt((1+B)/(1-B));
        if r>1
            SX = SX/r;
            SY = SY/r;
        else
            SZ = SZ*r;
        end
        o = [0,0,0];
        % Add wings
        n_wings = 5;
        b = 0.4;
        SX = SX.*(1+b*cos(n_wings*atan2(SY,SX)));
        SY = SY.*(1+b*cos(n_wings*atan2(SY,SX)));
        % Make the swimmer chiral
        if i_movie >2
            p = 0.8;
            th = SZ;
            SXt = SX;
            SX = SX.*cos(p*th) - SY.*sin(p*th);
            SY = SXt.*sin(p*th) + SY.*cos(p*th);
        end

        % Break fore-aft symmetry
        if i_movie >1
            zmax = 0.5;
            zmin = -0.5;
            azm = 0.2;
            azM = 0.6;
            % PZ = azm*sin(pi*SZ).*(SZ<0) + azM*sin(pi*SZ).*(SZ>0);
            PZ = zmin + (zmax-zmin)*(SZ+1)/2;
            SX = SX + SX.*PZ;
            SY = SY + SY.*PZ;
        end

        % Move the swimmer to its position
        SX = SX + x(i);
        SY = SY + y(i);
        SZ = SZ + z(i);

        % Rotate the swimmer
        o = [x(i),y(i),z(i)];

        s = surf(SX,SY,SZ,sqrt(r^2*(SX.^2+SY.^2))+1*SZ,'LineStyle','none','FaceLighting','gouraud');

        % Rotate the swimmer by pi/2 to make it aligned with x axis
        rotate(s,[0 1 0],90,o)

        % First, azimuthal angle phi.
        rotate(s,[1 0 0],rad2deg(phi(i)),o)
        % Then, polar angle theta.
        rotate(s,[0,cos(phi(i)),sin(phi(i))],rad2deg(theta(i)),o)
        % And finally proper rotation psi.
        rotate(s,[cos(theta(i)),sin(theta(i))*sin(phi(i)),-sin(theta(i))*cos(phi(i))],rad2deg(psi(i)),o)

        colormap(cols{i_movie});

        shading interp
        material shiny
        axis equal
        box on
        light
        camlight

        xlabel('$x$','Interpreter','latex')
        ylabel('$y$','Interpreter','latex')
        zlabel('$z$','Interpreter','latex')
        set(gca,'FontSize',24)
        set(gca,'TickLabelInterpreter','latex')

        axlim = 1.5;
        axis([x(i)-axlim x(i)+axlim y(i)-axlim y(i)+axlim z(i)-axlim z(i)+axlim]);

        %view(135,40)


        hold on

        if i_movie > 3

            % Full traj so far.
            plot3(x(1:i),y(1:i),z(1:i),'LineWidth',0.5,'Color','k')
            hold on
            % Avg traj so far.
            plot3(xr(1:i),yr(1:i),zr(1:i),'LineWidth',4,'Color',cols{i_movie}(end,:))

        else

            % Full traj so far.
            plot3(x(1:i),y(1:i),z(1:i),'LineWidth',4,'Color',cols{i_movie}(end,:))

        end


        nexttile(3);

        hold off

        % Plot the axis and (xy) unit circle
        plot3([0 1.5],[0 0],[0 0],'k','LineWidth',1)
        hold on
        plot3([0 0],[0 1.5],[0 0],'k','LineWidth',1)
        plot3([0 0],[0 0],[0 1.5],'k','LineWidth',1)
        tt=0:0.1:2*pi+0.1;
        plot3(cos(tt),sin(tt),zeros(1,length(tt)),'k','LineWidth',1)

        % Plot each trajectory in light colour
        for j = 1:5
            i_traj = 3 + 1*(j>3);
            tha = data_movie_trans_chiral{j,i_traj}(:,1);
            pha = data_movie_trans_chiral{j,i_traj}(:,2);
            if j==2
                plot3(sin(tha).*cos(pha),sin(tha).*sin(pha),cos(tha),'--','LineWidth',2,'Color',cols{j}(i_light,:));
            else
                plot3(sin(tha).*cos(pha),sin(tha).*sin(pha),cos(tha),'LineWidth',2,'Color',cols{j}(i_light,:));
            end
            hold on
        end

        if i_movie>3
            % Full traj so far.
            plot3(sin(theta(1:i)).*cos(phi(1:i)),sin(theta(1:i)).*sin(phi(1:i)),cos(theta(1:i)),'k')

            % Avg traj so far
            plot3(sin(thr(1:i)).*cos(phr(1:i)),sin(thr(1:i)).*sin(phr(1:i)),cos(thr(1:i)),'Color',cols{i_movie}(end,:),'LineWidth',3)
            plot3(sin(thr(1)).*cos(phr(1)),sin(thr(1)).*sin(phr(1)),cos(thr(1)),'k.','MarkerSize',30)

        else

            % Full traj so far.
            plot3(sin(theta(1:i)).*cos(phi(1:i)),sin(theta(1:i)).*sin(phi(1:i)),cos(theta(1:i)),'Color',cols{i_movie}(end,:),'LineWidth',3)

        end

        % Mark the initial conditions.
        plot3(sin(theta(1)).*cos(phi(1)),sin(theta(1)).*sin(phi(1)),cos(theta(1)),'k.','MarkerSize',30)

        view(120,25)
        % camlight
        grid on
        box on

        lim = 1.5;

        axis off equal
        axis([-lim lim -lim lim -lim lim])

        %title(['$\omega =',num2str(w(j_w)),'$'],'Interpreter','latex')
        set(gca,'FontSize',24)

        % save the frame
        %     if save_movie
        %         frame = getframe(gcf);
        %         writeVideo(v1,frame);
        %     end

        % tl.Padding = 'tight';
        %tl.TileSpacing = 'none';

        drawnow

        if save_movie
            frame = getframe(gcf);
            writeVideo(v1,frame);
        end
    end
    close(v1);
end
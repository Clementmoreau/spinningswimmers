%%%% Diagram showing the qualitative regimes of the rotational dynamics of
%%%% the spinning swimmer, depending on B, C, D and omega.

% This is a schematic diagram so a lot is done 'by hand'.

% setup the figure.
figure(1);clf;
set(gcf, 'Position',  [1, 640, 860, 470])
tl = tiledlayout(3,6);

% predefined colours.
red = [1,0.15,0];
blue = [0,0.4,0.8];
purple = [0.35,0,0.95];

% the number of sides for the ellipses we will draw.
Nsides = 72;

% three different values of the chirality parameter D.
D_var = [0;0.3;-0.3];

% loop on the values of D.
for i = 1:3

    D = D_var(i);

    % now, we will loop over four values of omega, ignoring the special
    % cases omega = sqrt(2) and omega = infinity for now.

    w_var = [0,1,sqrt(3),2];
    w_var_str = {'0','1','\sqrt{3}','2'};
    j_w = [1,2,4,5];

    % loop over omega.
    for j = 1:length(j_w)

        w = w_var(j);

        % go to the right tile.
        nexttile(6*(i-1)+j_w(j));

        % define a circle and resize it to an ellipse with respect to
        % B_hat, C_hat, D_hat.
        pgon = nsidedpoly(Nsides);
        pgon = rotate(pgon,180/Nsides); % small rotation to allow accurate division into quarters
        Binv = abs(2*(1+w^2)/(2-w^2)); % effective B
        Cinv = abs((1+w^2)^(3/2)*(1+D*w^2)); % effective C (depending on D as well)
        pgon.Vertices(:,1) = pgon.Vertices(:,1)*Binv;
        pgon.Vertices(:,2) = pgon.Vertices(:,2)*Cinv;

        % divide the ellipse into four quarters.
        quarter1 = polyshape([pgon.Vertices(1:(Nsides/4+1),1);0],[pgon.Vertices(1:(Nsides/4+1),2);0]);
        quarter2 = polyshape([pgon.Vertices((Nsides/4+1):(Nsides/2+1),1);0],[pgon.Vertices((Nsides/4+1):(Nsides/2+1),2);0]);
        quarter3 = polyshape([pgon.Vertices((Nsides/2+1):(3*Nsides/4+1),1);0],[pgon.Vertices((Nsides/2+1):(3*Nsides/4+1),2);0]);
        quarter4 = polyshape([pgon.Vertices((3*Nsides/4+1):end,1);pgon.Vertices(1,1);0],[pgon.Vertices((3*Nsides/4+1):end,2);pgon.Vertices(1,2);0]);

        % define the four quadrants as well.
        square1 = polyshape([-15 -15.5 0 0],[0 15 15 0]);
        square2 = polyshape([0 0 15 15],[0 15 15 0]);
        square3 = polyshape([-15 -15 0 0],[0 -15 -15 0]);
        square4 = polyshape([0 0 15 15],[0 -15 -15 0]);

        % define the proper colours. red means that the attractors are on
        % the south pole or hemisphere and blue that they are on the north
        % hemisphere. when omega > sqrt(2) (j > 2 in the code), B_hat has the
        % opposite sign as B, so the colours are reversed. 
        % In the case omega = 2 and D = -0.3 (j==4 && i==3), C_hat also has 
        % the opposite sign as C because of the negative value of D, 
        % so the colours are back to initial layout.
        if j <= 2 || (j==4 && i==3)
            col1 = red;
            col2 = blue;
        else
            col1 = blue;
            col2 = red;
        end

        % plot all the shapes.
        plot(pgon,'LineWidth',0.5,'EdgeColor','k','FaceAlpha',0)
        hold on
        plot(square1,'FaceColor',col1,'FaceAlpha',0.2)
        plot(square4,'FaceColor',col1,'FaceAlpha',0.2)

        plot(square2,'FaceColor',col2,'FaceAlpha',0.2)
        plot(square3,'FaceColor',col2,'FaceAlpha',0.2)

        plot(quarter1,'FaceColor',col2,'FaceAlpha',0.4)
        plot(quarter3,'FaceColor',col2,'FaceAlpha',0.4)

        plot(quarter2,'FaceColor',col1,'FaceAlpha',0.4)
        plot(quarter4,'FaceColor',col1,'FaceAlpha',0.4)

        % plot lines and dashed lines indicating periodical trajectories.
        plot([0 0],[-Cinv Cinv],'LineWidth',4,'Color',col3)
        plot([-Binv Binv],[0 0],'LineWidth',4,'Color',col3)

        plot([0 0],[-15 -Cinv],':','LineWidth',3,'Color',col3)
        plot([0 0],[Cinv 15],':','LineWidth',3,'Color',col3)
        plot([-15 -Binv],[0 0],':','LineWidth',3,'Color',col3)
        plot([Binv 15],[0 0],':','LineWidth',3,'Color',col3)

        % add value of omega on the first row.
        if i==1
            title(['$\omega =',w_var_str{j},'$'],'Interpreter','latex')
        end
    end

    % special case w = sqrt(2)
    nexttile(6*(i-1)+3);

    % here we divide the (B,C) plane into three horizontal bands. The
    % height of the central band depends on D.
    square1 = polyshape([-15 -15 15 15],[5*(1+2*D) 15 15 5*(1+2*D)]);
    square2 = polyshape([-15 -15 15 15],[-5*(1+2*D) 5*(1+2*D) 5*(1+2*D) -5*(1+2*D)]);
    square3 = polyshape([-15 -15 15 15],[-15 -5*(1+2*D) -5*(1+2*D) -15]);

    hold on
    plot(square1,'FaceColor',col3,'FaceAlpha',0.2)
    plot(square3,'FaceColor',col3,'FaceAlpha',0.2)
    plot(square2,'FaceColor',col3)

    if i==1
        title('$\omega = \sqrt{2}$','Interpreter','latex')
    end

    % special case w = inf

    nexttile(6*(i-1)+6);

    % here we divide the (B,C) plane into three vertical bands, independent
    % of D. 
    square1 = polyshape([-15 -2 -2 -15],[-15 -15 15 15]);
    square2 = polyshape([-2 2 2 -2],[-15 -15 15 15]);
    square3 = polyshape([2 15 15 2],[-15 -15 15 15]);

    hold on

    plot(square1,'FaceColor',col3,'FaceAlpha',0.2)
    plot(square3,'FaceColor',col3,'FaceAlpha',0.2)
    plot(square2,'FaceColor',col3)

    if i==1
        title('$\omega \rightarrow \infty$','Interpreter','latex')
    end

    % graphical setup.
    ax = [-12 12 -12 12];
    for j = 1:6
        nexttile((i-1)*6+j);
        box on
        axis equal
        axis(ax)
        set(gca,'FontSize',18)
        set(gca,'TickLabelInterpreter','latex')
        if i == 3
            xlabel('$B$','Interpreter','latex');
        else
            xticks([]);
        end
        if j == 1
            ylabel('$C$','Interpreter','latex');
        else
            yticks([]);
        end
    end

end

tl.Padding = 'none';
tl.TileSpacing = 'compact';

% export the figure.
% exportgraphics(gcf,'bifurcation_diagram.eps','ContentType','vector')



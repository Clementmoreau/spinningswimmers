%%%%%% Figure showing the B_hat, C_hat, D_hat  with respect to omega. %%%%%%

% setup the figure
figure(1);clf;
set(gcf, 'Position',  [1, 640, 1100, 250])
tl = tiledlayout(1,3);

% range of omega
wmax = 10;
w = -wmax:0.1:wmax;

% first, we plot a figure for B_hat.
nexttile(1);

B_hat = (2-w.^2)./(2*(1+w.^2));

% for the colours, matching the ones of the translational dynamics figure.
col_beta = colormap(parula(256));

plot(w,B_hat,'k','LineWidth',1.5,'Color',col_beta(end,:))

% graphical setup.
grid on
box on
set(gca,'FontSize',18)
set(gca,'TickLabelInterpreter','latex')
xlabel('$\omega$','Interpreter','latex')
ylabel('$\hat{B}/B, \hat{\beta}/\beta$','Interpreter','latex')
axis([-wmax wmax -0.5 1])

% second panel is C_hat, for different values of D/C.
nexttile(2);hold on

% values of D/C and associated legend names
D_var = [-1 0 sqrt(2) 3];
D_str = {'-1','0','\sqrt{2}','3'};

col_gamma = colormap(parula(256));
N_col = 256;
i_col = [60,120,200,250];

% plot
for i = 1:length(D_var)

    C_hat = (1 + w.^2*D_var(i))./(1 + w.^2).^(3/2);
    plot(w,C_hat,'LineWidth',1.5,'Color',col_gamma(i_col(i),:),'DisplayName',['$D/C = ',D_str{i},'$'])
end

% also, plot the particular case C=0.
% plot(w,w.^2./(1 + w.^2).^(3/2),'k:','LineWidth',1.5,'DisplayName','$C = 0, D = 1$');

% graphical setup.
grid on
box on
set(gca,'FontSize',18)
set(gca,'TickLabelInterpreter','latex')
xlabel('$\omega$','Interpreter','latex')
ylabel('$\hat{C}/C, \hat{\gamma}/\gamma$','Interpreter','latex')
axis([-wmax wmax -0.5 1.5])
lgd = legend;
lgd.Interpreter = 'latex';
lgd.FontSize = 10;
lgd.Location = 'northwest';

% third panel is D_hat for different values of C/D.
nexttile(3);
hold on

% values ofC/D and associated legend names
C_var = [-1 1/3 4/3 3];
C_str = {'-1','1/3','4/3','3'};

col_delta = colormap(parula(256));
N_col = 256;
i_col = [60,120,200,250];

% plot
for i = 1:length(C_var)
    D_hat = (3*w.^2*C_var(i) + (2-w.^2))./2./(1 + w.^2).^(3/2);
    plot(w,D_hat,'LineWidth',1.5,'Color',col_delta(i_col(i),:),'DisplayName',['$C/D = ',C_str{i},'$'])
end

% also, plot the particular case D=0.
% plot(w,(3*w.^2)./2./(1 + w.^2).^(3/2),'k:','LineWidth',1.5,'DisplayName','$D = 0, C = 1$');

% graphical setup.
grid on
box on
set(gca,'FontSize',18)
set(gca,'TickLabelInterpreter','latex')
xlabel('$\omega$','Interpreter','latex')
ylabel('$\hat{D}/D, \hat{\delta}/\delta$','Interpreter','latex')
axis([-wmax wmax -1 2.2])
lgd = legend;
lgd.Interpreter = 'latex';
lgd.FontSize = 10;
lgd.Location = 'northwest'

tl.Padding = 'none';
tl.TileSpacing = 'compact';

% export the figure.
exportgraphics(gcf,'bhatchatdhat.eps','ContentType','vector')




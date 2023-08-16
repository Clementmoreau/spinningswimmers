function [r_effective,r,w] = effective_aspect_ratio
%EFFECTIVE_ASPECT_RATIO 
%Plots effective aspect ratio for spinning swimmer with respect to the
%aspect ratio and spinning angle alpha.

% LineWidth
L = 6;
% Font size
F = 50;
F2 = 40;
% F2 = 40;

w = linspace(0,10,1e3);
r_list = [0.1, 0.5, 1, 2, 10];

f1 = figure(1);
% Resize figure window
f1Pos = f1.Position;
f1Pos(3) = 982;
f1Pos(4) = 982;
f1.Position = f1Pos;

semilogy(2*sqrt(2),1e1,'HandleVisibility','off');
hold on

grey = [0.9, 0.9, 0.9];

area([sqrt(2),4],[2e1 2e1],'edgecolor','none','facecolor',grey,'HandleVisibility','off');
% area([-3*sqrt(2)+0.01,-sqrt(2)],[2e1 2e1],'edgecolor','none','facecolor',grey,'HandleVisibility','off');
set(gca, "Layer", "top")

C_divide = floor(linspace(1,256,length(r_list)));
colour_list = ametrine(256);
C_list = colour_list(C_divide,:);


for j = 1:length(r_list)
    r = r_list(j);
    r_effective = sqrt((4*r.^2 + w.^2.*(3 + r.^2))./(4 + w.^2.*(1 + 3.*r.^2)));
    semilogy(w,r_effective,'LineWidth',L,'Color',C_list(j,:));
    hold on
end
% legend('$r = 0.1$','$r = 0.5$','$r = 1$','$r = 2$','$r = 10$','Interpreter','latex')

% set(gca, "Layer", "top")

% plot([-sqrt(2),-sqrt(2)],[0 2e1],'k--','LineWidth',L-1,'HandleVisibility','off');
plot([sqrt(2),sqrt(2)],[1e-2 2e1],'k--','LineWidth',L-1,'HandleVisibility','off');

annotation(f1,'textbox',...
    [0.193464358452136 0.846451612903226 0.147676171079432 0.0683870967741942],...
    'String','$\omega < \sqrt{2}$',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',F2,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');


annotation(f1,'textbox',...
    [0.593668024439911 0.859354838709677 0.147676171079436 0.0529032258064522],...
    'String','$\omega > \sqrt{2}$',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',F2,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

annotation(f1,'textbox',...
    [0.549879837067198 0.667096774193548 0.23830753564155 0.04],...
    'String','Increasing $r$',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',F2,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

annotation(f1,'arrow',[0.667006109979626 0.667006109979633],...
    [0.646451612903226 0.438709677419355],'LineWidth',4,'HeadWidth',20,...
    'HeadLength',20);


axis([0 4 0.8e-1 1.2e1])
xticks([0;2;4]);
xticklabels({'$0$', '$2$', '$4$'});
yticks([0.1, 0.5, 1, 2, 10]);
% yticklabels({'$-5$', '$0$', '$5$'})
set(gcf, 'color', 'w'); 
ax = gca;
ax.FontSize = F;
ax.TickLabelInterpreter = 'latex';
xlabel('$\omega$','FontSize',F,'Interpreter','latex');%,'Position',[0,-3.6,0]);
ylabel('$\hat{r}$','FontSize',F,'Interpreter','latex','Rotation',0);%,'Position',[-5.8,0,0]);


w_list_pre = [0.1, 0.5, 1];%, sqrt(2), 4];
w_list_post = [2, 4];
r = logspace(-2,2,1e3);

f2 = figure(2);
% Resize figure window
f2Pos = f2.Position;
f2Pos(3) = 982;
f2Pos(4) = 982;
f2.Position = f2Pos;


loglog(2*sqrt(2),1e1,'HandleVisibility','off');
hold on

r_effective_inf = sqrt(((3 + r.^2))./((1 + 3.*r.^2)));

area([1e0,2e1],[1e0 1e0],'edgecolor','none','facecolor',grey,'HandleVisibility','off');
area(r(r>9e-1),r_effective_inf(r>9e-1),'edgecolor','none','facecolor',[1,1,1],'HandleVisibility','off');

area(r(r<1e0),r_effective_inf(r<1e0),'edgecolor','none','facecolor',grey,'HandleVisibility','off');
area([5e-2,1e0],[1e0 1e0],'edgecolor','none','facecolor',[1,1,1],'HandleVisibility','off');

% area([1e0,1e2],[1e0 1e2],'edgecolor','none','facecolor',grey,'HandleVisibility','off');
% area([1e-2,1e0],[1e0 1e0],'edgecolor','none','facecolor',grey,'HandleVisibility','off');
% area([1e-2,1e0,1e2],[1e-2 1e0 1e0],'edgecolor','none','facecolor',[1,1,1],'HandleVisibility','off');
% area([-3*sqrt(2)+0.01,-sqrt(2)],[2e1 2e1],'edgecolor','none','facecolor',grey,'HandleVisibility','off');
set(gca, "Layer", "top")



C_divide = floor(linspace(1,256,length(w_list_pre)+length(w_list_post)));
colour_list = (isolum(256));
C_list = colour_list(C_divide,:);

% w = 0;
% r_effective = sqrt((4*r.^2 + w.^2.*(3 + r.^2))./(4 + w.^2.*(1 + 3.*r.^2)));
%     loglog(r,r_effective,'k:','LineWidth',L-1);

for j = 1:length(w_list_pre)
    w = w_list_pre(j);
    r_effective = sqrt((4*r.^2 + w.^2.*(3 + r.^2))./(4 + w.^2.*(1 + 3.*r.^2)));
    loglog(r,r_effective,'LineWidth',L,'Color',C_list(j,:));
    hold on
end

w = sqrt(2);
r_effective = sqrt((4*r.^2 + w.^2.*(3 + r.^2))./(4 + w.^2.*(1 + 3.*r.^2)));
    loglog(r,r_effective,'k--','LineWidth',L-1);

for j = length(w_list_pre)+1:length(w_list_pre) + length(w_list_post)
    w = w_list_post(j-length(w_list_pre));
    r_effective = sqrt((4*r.^2 + w.^2.*(3 + r.^2))./(4 + w.^2.*(1 + 3.*r.^2)));
    loglog(r,r_effective,'LineWidth',L,'Color',C_list(j,:));
    hold on
end

w = 0;
r_effective = sqrt((4*r.^2 + w.^2.*(3 + r.^2))./(4 + w.^2.*(1 + 3.*r.^2)));
    loglog(r,r_effective,'k:','LineWidth',L-1,'HandleVisibility','off');

r_effective = sqrt(((3 + r.^2))./((1 + 3.*r.^2)));
loglog(r,r_effective,'k:','LineWidth',L-1);


% legend('$\omega = 0$','$\omega = 0.1$','$\omega = 0.2$','$\omega = \sqrt{2}$','$\omega = 10$','$\omega = \infty$','Interpreter','latex','Location','northwest')


axis([5e-2 2e1 5e-2 2e1])
xticks([1e-1;1e0;1e1]);
xticklabels({'$10^{-1}$', '$10^{0}$',  '$10^{1}$'});
yticks([1e-1;1e0;1e1]);
yticklabels({'$10^{-1}$',  '$10^{0}$', '$10^{1}$'});
set(gcf, 'color', 'w'); 
ax = gca;
ax.FontSize = F;
ax.TickLabelInterpreter = 'latex';
xlabel('$r$','FontSize',F,'Interpreter','latex');%,'Position',[0,-3.6,0]);
ylabel('$\hat{r}$','FontSize',F,'Interpreter','latex','Rotation',0);%,'Position',[-5.8,0,0]);

annotation(f2,'arrow',[0.450101832993888 0.289205702647655],...
    [0.353548387096774 0.705806451612903],'LineWidth',4,'HeadWidth',20,...
    'HeadLength',20);

annotation(f2,'textbox',...
    [0.349268839103864 0.318709677419354 0.23830753564155 0.04],...
    'String',{'Increasing $\omega$'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',F2,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

annotation(f2,'arrow',[0.601832993890008 0.762729124236241],...
    [0.745806451612903 0.393548387096774],'LineWidth',4,'HeadWidth',20,...
    'HeadLength',20);

annotation(f2,'textbox',...
    [0.474523421588586 0.769032258064515 0.238307535641549 0.04],...
    'String',{'Increasing $\omega$'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',F2,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

end


om = -5:0.01:5;
lam = sqrt(1+om.^2);
f = (2-om.^2)./(2*(1+om.^2));

figure(1);clf;
set(gcf, 'Position',  [1, 200, 900, 225])
tiledlayout(1,2)
nexttile;

plot(om,f,'k','LineWidth',2)
hold on
grid on
plot([-5 5],[0 0],'k');
plot([0 0],[-0.5 1],'k');
set(gca,'FontSize',18)
set(gca,'TickLabelInterpreter','latex')
%zticks([0 8])
xlabel('$\omega$','Interpreter','latex'); ylabel('$\hat{B}/B$','Interpreter','latex');

nexttile;
hold on
v2 = [0 0.5 1 3];
for i = 1:length(v2)
plot(om,(1+om*v2(i))./lam/sqrt(1+v2(i)^2),'LineWidth',2,'Color',[i i i]/6)
end

grid on
box on
set(gca,'FontSize',18)
set(gca,'TickLabelInterpreter','latex')
%zticks([0 8])
xlabel('$\omega$','Interpreter','latex'); ylabel('$\hat{V}/|V|$','Interpreter','latex');
legend('$V_2 = 0$','$V_2 = 0.5$','$V_2 = 1$','$V_2 = 3$','Location','southeast','Interpreter','latex','FontSize',14)


exportgraphics(gcf,'B_hat_V_hat.eps','ContentType','vector')
function out = fig_flow

global Yinit Yinit2 
global Ntheta Nphi 
global alpha beta

Ntheta=400;
Nphi=800;
load('contour_zero');
ContMat=Mat;

% ----- alpha,beta,Yinit,Yinit2 ---
FlowDat=[ 0.0, 0.8, -1.0,  -1.5;...
          0.2, 0.8, -1.21, -1.21;...
          0.5, 0.0, -1.0,  -1.25;...
          0.5,-0.5, -1.1,  -1.1;...
          1.5, 0.5, -1.0,  -1.0;...
          1.5, 0.0, -1.0,   0.0;...
          0.2, 1.5, -1.0,   0.0;];  

figure;
dpi = '-r400';

hold on
%viscircles([0,0],1,'Color','b','LineWidth',3);
LineX=linspace(-1,1,2);
LineY=linspace(-2,2,2);

plot (LineX,zeros(size(LineX)),'Color','k','LineWidth',3)
plot (zeros(size(LineY)),LineY,'Color','k','LineWidth',3)
plot (ContMat(1,2:size(ContMat,2)), ContMat(2,2:size(ContMat,2)),...
      'g--','LineWidth',2)
plot (-ContMat(1,2:size(ContMat,2)), ContMat(2,2:size(ContMat,2)),...
      'g--','LineWidth',2)
%viscircles([ 1,0],0.05,'Color','k','LineWidth',2);
%viscircles([-1,0],0.05,'Color','k','LineWidth',2);
s=scatter(LineX,zeros(size(LineX)),100,'w','filled');
    s.LineWidth = 1.0;
    s.MarkerEdgeColor = 'k';

s2=scatter(FlowDat(:,2),FlowDat(:,1),100,'r','filled');
    s2.LineWidth = 1.0;
    s2.MarkerEdgeColor = 'k';  
    
xlabel('$\beta$','Interpreter','latex')
ylabel('$\alpha$','Interpreter','latex')
xlim([-2 2])
xticks([-2:0.5:2])
ylim([0 2])
yticks([-2:0.5:2])

set(gca,'FontSize',24)
set(gca,'TickLabelInterpreter','latex')

hold off
grid on
daspect([1,1,1]);  

print('fig_flow_param','-depsc',dpi) 
           

% ----- draw flow ---      
for label=1:size(FlowDat,1)
    alpha =FlowDat(label,1);
    beta  =FlowDat(label,2);
    Yinit =FlowDat(label,3)*pi;
    Yinit2=FlowDat(label,4)*pi;
    fig_flow_draw(label);
end



end
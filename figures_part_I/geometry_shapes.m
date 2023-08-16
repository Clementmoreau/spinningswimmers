% Use to plot various shapes

figure(1);clf;
set(gcf, 'Position',  [1, 640, 800, 800])
set(gcf,'color','w');

% chirality ?
chiral = 1;
% fore-aft asymmetry ?
foreaft = 0;
% heterochirality ?
heterochirality = 0;
homochirality = 0;

savemovie = 0;

nb = 60;
p_var = linspace(-5,5,nb);
b_var = linspace(-0.5,1,nb);
zmax_var = linspace(0.5,1,nb);
zmin_var = linspace(-0.8,-0.2,nb);
r_var = linspace(0.2,4,nb);

fr = [1:length(p_var),length(p_var)-1:-1:1];

if savemovie
v1 = VideoWriter(['prettymovie.mp4'],'MPEG-4');
        v1.FrameRate = 30;
        v1.Quality = 50;
        open(v1);
end

for j = 1:length(fr)

    i = fr(j);

n_sphere = 200;
[SX,SY,SZ]=sphere(n_sphere);
% aspect ratio
r = r_var(i);
if r>1
    SX = SX/r;
    SY = SY/r;
else
    SZ = SZ*r;
end
o = [0,0,0];
% Add wings
n_wings = 5;
b = b_var(i);
SX = SX.*(1+b*cos(n_wings*atan2(SY,SX)));
SY = SY.*(1+b*cos(n_wings*atan2(SY,SX)));
% Make the swimmer chiral
if chiral

    p = p_var(i);
    th = SZ;
    SXt = SX;
    SX = SX.*cos(p*th) - SY.*sin(p*th);
    SY = SXt.*sin(p*th) + SY.*cos(p*th);
    % Break fore-aft symmetry

end

if foreaft

    zmax = zmax_var(i);
    zmin = zmin_var(i);
    azm = 0.2;
    azM = 0.6;
    % PZ = azm*sin(pi*SZ).*(SZ<0) + azM*sin(pi*SZ).*(SZ>0);
    PZ = zmin + (zmax-zmin)*(SZ+1)/2;
    SX = SX + SX.*PZ;
    SY = SY + SY.*PZ;

end

if heterochirality

    SX(1:n_sphere/2,:) = - SX(n_sphere+1:-1:n_sphere/2+2,:);
    SY(1:n_sphere/2,:) = SY(n_sphere+1:-1:n_sphere/2+2,:);

end

if homochirality

    SX(1:n_sphere/2,:) = SX(n_sphere+1:-1:n_sphere/2+2,:);
    SY(1:n_sphere/2,:) = SY(n_sphere+1:-1:n_sphere/2+2,:);

end

% Rotate the swimmer
theta = pi/2;
phi = 0;
psi = 0.2;
s = surf(SX,SY,SZ,sqrt(r^2*(SX.^2+SY.^2))+1.5*SZ,'LineStyle','none','FaceLighting','gouraud','FaceAlpha',1);
rotate(s,[0 1 0],rad2deg(theta),o)
rotate(s,[0 0 1],rad2deg(phi),o)
rotate(s,[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)],rad2deg(psi),o)

colormap parula;

shading interp
material shiny
axis off equal
camproj('perspective')
view(135+i/3,40+i/3)
camlight

drawnow
hold off

if savemovie
frame = getframe(gcf);
            writeVideo(v1,frame);
end

end

if savemovie
close(v1);
end



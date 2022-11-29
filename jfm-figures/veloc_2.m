function out = veloc_2(t,Xvar)

global alpha beta

out=zeros(2,1);
theta=Xvar(1);
phi=Xvar(2);


v_theta= 0.25*beta*sin(2*theta)*sin(2*phi)-0.5*alpha*sin(theta)*cos(2*phi);
v_phi=beta*cos(phi)*cos(phi) + 0.5*(1-beta)+alpha*cos(theta)*cos(phi)*sin(phi);

out=[v_theta; v_phi];

end


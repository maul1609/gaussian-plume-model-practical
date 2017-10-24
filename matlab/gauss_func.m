function C=gauss_func(Q,u,dir1,x,y,z,xs,ys,H,Dy,Dz,STABILITY)
u1=u;
x1=x-xs; % shift the coordinates so that stack is centre point
y1=y-ys; 

% components of u in x and y directions
wx=u1.*sin((dir1-180).*pi./180);
wy=u1.*cos((dir1-180).*pi./180);

% Need angle between point x, y and the wind direction, so use scalar product:
dot_product=wx.*x1+wy.*y1;
% product of magnitude of vectors:
magnitudes=u1.*sqrt(x1.^2+y1.^2); 

% angle between wind and point (x,y)
subtended=acos(dot_product./magnitudes);
% distance to point x,y from stack
hypotenuse=sqrt(x1.^2+y1.^2);

% distance along the wind direction to perpendilcular line that intesects
% x,y
downwind=cos(subtended).*hypotenuse;

% Now calculate distance cross wind.
crosswind=sin(subtended).*hypotenuse;

ind=find(downwind>0);
C=zeros(size(downwind));

% sig_y=sqrt(2.*Dy.*downwind./u1);
% sig_z=sqrt(2.*Dz.*downwind./u1);

% calculate sigmas based on stability and distance downwind
[sig_y,sig_z]=calc_sigmas(STABILITY,downwind);
    
C(ind)=Q./(2.*pi.*u1.*sig_y(ind).*sig_z(ind)) ...
    .*exp(-crosswind(ind).^2./(2.*sig_y(ind).^2)).* ... % note that a previous version had (2.*sig_y(ind)).^2 thanks for John Marino for pointing out
    (exp(-(z(ind)-H).^2./(2.*sig_z(ind).^2))+ ...
    exp(-(z(ind)+H).^2./(2.*sig_z(ind).^2)) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GAUSSIAN PLUME MODEL FOR TEACHING PURPOSES                              %
% PAUL CONNOLLY (UNIVERSITY OF MANCHESTER, 2014)                          %
% THIS CODE IS PROVIDED `AS IS' WITH NO GUARANTEE OF ACCURACY             %
% IT IS USED TO DEMONSTRATE THE EFFECTS OF ATMOSPHERIC STABILITY,         %
% WINDSPEED AND DIRECTION AND MULTIPLE STACKS ON THE DISPERSION OF        %
% POLLUTANTS FROM POINT SOURCES                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not change these variables                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SECTION 0: Definitions (normally don't modify this section)
% view
PLAN_VIEW=1;
HEIGHT_SLICE=2;
SURFACE_TIME=3;
NO_PLOT=4;

% wind field
CONSTANT_WIND=1;
FLUCTUATING_WIND=2;
PREVAILING_WIND=3;

% number of stacks
ONE_STACK=1;
TWO_STACKS=2;
THREE_STACKS=3;

% stability of the atmosphere
CONSTANT_STABILITY=1;
ANNUAL_CYCLE=2;
stability_str={'Very unstable','Moderately unstable','Slightly unstable',...
    'Neutral','Moderately stable','Very stable'};
% Aerosol properties
HUMIDIFY=2;
DRY_AEROSOL=1;

SODIUM_CHLORIDE=1;
SULPHURIC_ACID=2;
ORGANIC_ACID=3;
AMMONIUM_NITRATE=4;
nu=[2 2.5 1 2];
rho_s=[2160 1840 1500 1725];
Ms=[58.44e-3 98e-3 200e-3 80e-3];
Mw=18e-3;


dxy=100;          % resolution of the model in both x and y directions
dz=10;
x=-2500:dxy:2500; % solve on a 5 km domain
y=x;              % x-grid is same as y-grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% SECTION 1: Configuration
% Variables can be changed by the user+++++++++++++++++++++++++++++++++++++
RH=0.90;
aerosol_type=SODIUM_CHLORIDE;

dry_size=60e-9;
humidify=DRY_AEROSOL;

stab1=1; % set from 1-6
stability_used=CONSTANT_STABILITY;


output=PLAN_VIEW;
x_slice=26; % position (1-50) to take the slice in the x-direction
y_slice=1;  % position (1-50) to plot concentrations vs time

wind=PREVAILING_WIND;
stacks=ONE_STACK;
stack_x=[0 1000 -200];
stack_y=[0 250 -500];

Q=[40 40 40]; % mass emitted per unit time
H=[50 50 50]; % stack height, m
days=50;          % run the model for 365 days
%--------------------------------------------------------------------------
times=[1:days.*24]./24;

Dy=10;
Dz=10;

% SECTION 2: Act on the configuration information

% Decide which stability profile to use
switch stability_used
    case CONSTANT_STABILITY
        stability=stab1.*ones(days.*24,1);
        stability_str=stability_str{stab1};
    case ANNUAL_CYCLE
        stability=round(2.5.*cos(times.*2.*pi/(365))+3.5);
        stability_str='Annual cycle';
    otherwise
       return
end

% decide what kind of run to do, plan view or y-z slice, or time series
switch output
    case {PLAN_VIEW,SURFACE_TIME,NO_PLOT}
        C1=zeros(length(x),length(y),days.*24); % array to store data, initialised to be zero

        [x,y]=meshgrid(x,y); % x and y defined at all positions on the grid
        z=zeros(size(x));    % z is defined to be at ground level.
    case HEIGHT_SLICE
        z=0:dz:500;       % z-grid

        C1=zeros(length(y),length(z),days.*24); % array to store data, initialised to be zero

        [y,z]=meshgrid(y,z); % y and z defined at all positions on the grid
        x=x(x_slice).*ones(size(y));    % x is defined to be x at x_slice       
    otherwise
        return;
end
        

% Set the wind based on input flags++++++++++++++++++++++++++++++++++++++++
wind_speed=5.*ones(days.*24,1); % m/s
switch wind
    case CONSTANT_WIND
        wind_dir=0.*ones(days.*24,1);
        wind_dir_str='Constant wind';
    case FLUCTUATING_WIND
        wind_dir=360.*rand(days.*24,1);
        wind_dir_str='Random wind';
    case PREVAILING_WIND
        wind_dir=-sqrt(2).*erfcinv(2.*rand(days.*24,1)).*40; %norminv(rand(days.*24,1),0,40);
        % note at this point you can add on the prevailing wind direction, i.e.
        % wind_dir=wind_dir+200;
        wind_dir(find(wind_dir>=360))=...
            mod(wind_dir(find(wind_dir>=360)),360);
        wind_dir_str='Prevailing wind';
    otherwise 
        return;
end
%--------------------------------------------------------------------------


% SECTION 3: Main loop
% For all times...
h = waitbar(0,'Please wait...');
warning off
for i=1:length(wind_dir)
    for j=1:stacks
        C=gauss_func(Q(j),wind_speed(i),wind_dir(i),x,y,z,...
            stack_x(j),stack_y(j),H(j),Dy,Dz,stability(i));
        C1(:,:,i)=C1(:,:,i)+C;
    end
    waitbar(i/length(wind_dir),h);
end
warning on;
close(h);



% SECTION 4: Post process / output

% decide whether to humidify the aerosol and hence increase the mass
switch humidify
    case DRY_AEROSOL
        disp('do nothing');
    case HUMIDIFY
        mass=pi./6.*rho_s(aerosol_type).*dry_size.^3;
        moles=mass./Ms(aerosol_type);
        
        nw=RH.*nu(aerosol_type)*moles./(1-RH);
        mass2=nw.*Mw+moles.*Ms(aerosol_type);
        C1=C1.*mass2./mass; 
        
    otherwise
        return
end

% output the plots
switch output
    case PLAN_VIEW
        figure;
        pcolor(x,y,mean(C1,3).*1e6);shading flat
        caxis([0 1e2]);
        title({stability_str,wind_dir_str});
        xlabel('x (metres)');
        ylabel('y (metres)');
        h=colorbar;
        ylabel(h,'\mu g m^{-3}');
    case HEIGHT_SLICE
        figure;
        pcolor(y,z,mean(C1,3).*1e6);shading flat        
        caxis([0 1e2]);
        xlabel('y (metres)');
        ylabel('z (metres)');
        title({stability_str,wind_dir_str});
        h=colorbar;
        ylabel(h,'\mu g m^{-3}');
    case SURFACE_TIME
        figure;
        subplot(211);
        plot(times,1e6.*squeeze(C1(y_slice,x_slice,:)));hold on;
        try
            plot(times,smooth(1e6.*squeeze(C1(y_slice,x_slice,:)),24),'r');
            legend('Hourly mean','Daily mean')
        catch
        end
        xlabel('time (days)');
        ylabel('Mass loading (\mu g m^{-3})');
        title({stability_str,wind_dir_str});
        
        subplot(212)
        plot(times,stability);
        xlabel('time (days)');
        ylabel('Stability parameter');
        
    case NO_PLOT
        disp('don''t plot');
    otherwise
        return;
end

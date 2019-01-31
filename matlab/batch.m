wind_speeds=1:1:10;
q_strength=logspace(-3,0,5);


nnn=length(wind_speeds);
mmm=length(q_strength);

mytable=zeros(nnn,mmm);

for iii=1:nnn
    for jjj=1:mmm
        wind_use=wind_speeds(iii);
        Q_use=q_strength(jjj);
        gaussian_plume_model;
        mytable(iii,jjj)=C1(1,51+7).*1e6; % this is in microgram per cubic metre
    end
end

function[op]=fn_bemt_isolated(ip)
%% decoding parameters

theta75_deg = ip.theta;

Nb			= ip.Nb;
rad_ft		= ip.rad;
rad_in_ft	= ip.rad_in;
sigma_eq	= ip.sigma;
thetatw_deg	= ip.twist;
taper		= ip.taper;
sweep1_deg	= ip.sweep1;
sweep2_deg	= ip.sweep2;
rs			= ip.rs;

vc_knots	= ip.vc;
rho			= ip.rho; 

tc			= ip.tc;
tw1			= ip.tw1;

data		= ip.data;
soln		= ip.soln; 
use_soln	= ip.use_soln;
tiploss		= ip.tiploss;
itermax		= ip.itermax_bemt;
iter_cutoff	= ip.cutoff_bemt;
rw			= ip.rw;

%% unit conversion

theta75 = theta75_deg*(pi/180);		% radian
thetatw = thetatw_deg*(pi/180);		% radian
rad		= rad_ft*12*0.0254;			% meter
rad_in	= rad_in_ft*12*0.0254;		% meter
vc		= vc_knots*0.514444;		% meter/sec
% vc		= vc_knots;		% meter/sec
sweep1	= sweep1_deg*(pi/180);		% radian
sweep2	= sweep2_deg*(pi/180);		% radian

%% rpm

disc_area	= pi*rad^2;		% disc area
yct = ip.yct;
if (yct == 0)
	rpm = ip.rpm;
end
if (yct == 1)
	ctbysigma = ip.ctbysigma;
	rpm = (60/(2*pi*rad))*sqrt(ip.treq*4.44822/(ip.ctbysigma*sigma_eq*rho*disc_area*2));
end

omega		= rpm*2*pi/60;	% rotational speed
vtip		= omega*rad;    % tip speed (without considering fluid speed)
Mtip		= vtip/343;		% Tip Mach number
freq		= rpm/60*Nb;	% blade passage frequency
lambdac		= vc/vtip;      % non-dimensional climb velocity

%% bemt discretization
% N = number of panels,r = position of mid-point of each panels, dr = length of each panels

[N,r,dr]= dicretization();  
cutoff	= round(rad_in/rad*N);    % cut-off percentage
if (cutoff == 0)
    cutoff = 1;
end
if (cutoff > 1)
	rc = (r(cutoff)+r(cutoff-1))*0.5;
end
if (cutoff == 1)
    rc = 0;
end

%% calculate sectional configurations   

thetatw2= thetatw/(1-r(cutoff));
chd_eq	= sigma_eq*disc_area/Nb/(rad);									% effective chord              
chd0	= chd_eq/ (1+ (1/taper-1)*(0.75-r(cutoff))/ (1.0-r(cutoff)) );	% root chord
chdtw	= chd0*(1/taper-1)/(1-r(cutoff));								% change in chord
chd1	= chd0+((1-r(cutoff)))*chdtw;									% 
AR		= (rad-rad_in)/(chd0+chd1)*2;									% aspect ratio

theta	= zeros(N,1);
chd		= zeros(N,1);
sigma	= zeros(N,1);

for i = cutoff:N	
	chd(i,1)	= chd0 + chdtw*(r(i)-r(cutoff));			% sectional chord
	sigma(i,1)	= Nb*chd(i)*rad/disc_area;					% sectional solidity
	theta(i,1)	= theta75 + thetatw2 * (r(i)-0.75);			% sectional pitch
	if(r(i) <= tc)
		theta(i,1) = theta75+ thetatw2 * (tc-0.75)+ tw1*thetatw2 *(r(i)-tc);
	end
end

% calculate tangential distance of blade elements
rd1 = rs*tan(sweep1);
rd	= zeros(N,1);
for i=1:N
	if(r(i) <= rs)
		rd(i,1) = r(i)*tan(sweep1);
	end
	if(r(i) > rs)
		rd(i,1) = rd1-(r(i)-rs)*tan(sweep2);
	end
end

%% BEMT loop

cl			= zeros(N,1);
lambda		= zeros(N,1);
lambda2		= zeros(N,1);
lambda3		= zeros(N,1);
lambdai		= zeros(N,1);
phi			= zeros(N,1);
alpha		= zeros(N,1);
cl2d		= zeros(N,1);
cd			= zeros(N,1);
cm			= zeros(N,1);
dct			= zeros(N,1);
dcp_induced	= zeros(N,1);
dcp_profile = zeros(N,1);
dct2		= zeros(N,1);
dct3		= zeros(N,1);
dcmp		= zeros(N,1);
dcmp2		= zeros(N,1);
f			= zeros(N,1);
froot		= zeros(N,1);
F			= zeros(N,1);
vres		= zeros(N,1);
phi2		= zeros(N,1);
pm			= zeros(N,1);
inplane		= zeros(N,1);
outplane	= zeros(N,1);
check		= zeros(N,1);

mat_error	= zeros(itermax,1);
mat_ct		= zeros(itermax,1);
 
for i = 1:N
    F(i,1) = 1;			% tiploss function
end

ct = 0;
iter = 0;
flag = 0;

% use solution file if required
if (use_soln == 1)
	cl		= soln(:,1);
	lambda	= soln(:,2);
	lambda2 = lambda;
	F		= soln(:,3);
end


%%

while ((flag<1)&& (iter<itermax) )
    iter = iter+1;
    for i = cutoff:N
        if (r(i)>rw)
            lambdac=0;
        end
        if (r(i)<=rw)
            lambdac=vc/vtip;
		end 
        c_lambda = sigma(i)*cl(i)*r(i)/8/F(i)+lambda(i)*lambdac;
        if (c_lambda<0)
            c_lambda = 0;
        end
        lambda3(i,1) = lambda2(i);
        lambda2(i,1) = lambda(i);
        lambda(i,1) = sqrt(c_lambda);
%         lambda(i,1) = sqrt(ct/2);
        lambda(i) = (lambda(i)+lambda2(i)+lambda3(i))/3;
        lambdai(i,1) = lambda(i)-lambdac;
		phi(i,1) = atan(lambda(i)/r(i));
        alpha(i,1) = theta(i)-phi(i);       % sectional angle of attack
		if (ip.af == 1)
			[cl2d(i,1),cd(i,1),cm(i,1)] = fn_clcdcm(alpha(i),data);
		end
		if(ip.af == 0)
			cl2d(i,1) = ip.cl_alpha*alpha(i);
			cd(i,1) = ip.cd(1)+ip.cd(2)*alpha(i)+ip.cd(3)*alpha(i)^2;
		end
        cl(i,1) = cl2d(i,1)/(1+1/AR);
		dct(i,1) = 0.5*sigma(i)*r(i)^2*dr*( cl(i)*cos(phi(i)) -cd(i)*sin(phi(i)) );	% sectional thrust coefficient    
        dcp_induced(i,1) = 0.5*sigma(i)*cl(i)*r(i)^3*dr*sin(phi(i));				% sectional induced power coefficient
        dcp_profile(i,1) = 0.5*sigma(i)*cd(i)*r(i)^3*dr*cos(phi(i));				% sectional profile power coefficient
		
 
        dct2(i,1) = 0.5*r(i)^2*dr;
        dct3(i,1) = 0.5*sigma(i)*cd(i)*r(i)^2*dr+lambda(i,1)*dct(i)/r(i);
        dcmp(i,1) = 0.5*sigma(i)*cm(i)*r(i)^2*dr*chd(i)/rad;
        dcmp2(i,1) = dcmp(i)+dct(i)*rd(i);
        % tiploss effect
        f(i,1) = Nb*(1-r(i))/lambda(i)/2;    
        froot(i,1) = 0.5*Nb*(r(i)-rc)*r(i)/(1-r(i))/lambda(i);

        if (tiploss ==1)
%         F(i,1) = 2*acos( exp(-f(i)) )/pi;
			F(i,1) = (2*acos( exp(-f(i)) )/pi);
			Froot(i,1) = (2*acos( exp(-froot(i)) )/pi);
			F(i,1) =F(i,1)*Froot(i,1);
        end
        % tiploss effect ends
        vres(i,1) = sqrt((r(i))^2+(lambda(i))^2);
        phi2(i,1) = theta(i)-alpha(i);
        pm(i,1) = dct(i)*rd(i);
		inplane(i,1) = 0.5*rho*vres(i)^2*dr*rad*chd(i)*(cl(i)*sin(phi2(i)) + cd(i)*cos(phi2(i)));
        outplane(i,1) = 0.5*rho*vres(i)^2*dr*rad*chd(i)*(cl(i)*cos(phi2(i)) - cd(i)*sin(phi2(i)))/1;
    end
    ct_old = ct;
    ct = sum(dct);  % net thrust coefficient
%     theta_tip = theta_tip + 4*(ct_req-ct)/(sigma_eq*cl_alpha) + sqrt(0.5)*( sqrt(ct_req) -sqrt(ct) );   % updated pitch/collective
    error = abs((ct_old-ct)/ct*100);     % check error
    error2 = norm(lambda-lambda2);
    if((error<iter_cutoff) &&(iter_cutoff))
        flag = 1;
    end
  
end

%% final performance calculation

% PM = sum(pm);
ctbysigma = ct/sigma_eq;
cp = sum(dcp_induced)+sum(dcp_profile);		% coefficient of power

Thrust_N = ct*rho*disc_area*(omega*rad)^2;
Power_W  = cp*rho*disc_area*(omega*rad)^3;
Torque_Nm = Power_W/omega;

Thrust_lbs = Thrust_N*0.224809;				% thrust in lbs
Power_hp = Power_W/746;						% power in hp
FM = ct^1.5/sqrt(2)/cp;						% Figure fo merit
eta = Thrust_N*vc/Power_W;					% propulsive efficiency

theta_deg = theta.*180/pi;					% sectional pitch in deg.
alpha_deg = alpha.*180/pi;					% sectional AOA
vi = vtip*(lambdai);						% inflow speed

kk = sqrt(2)*sum(dcp_induced)/ct^1.5;

dt(:,1) = dct*rho*disc_area*(omega*rad)^2;
dt(:,2) = dct3*rho*disc_area*(omega*rad)^2;
dt(:,3) = dcmp*rho*disc_area*(omega*rad)^2*rad;
dt(:,4) = dcmp2*rho*disc_area*(omega*rad)^2*rad;
xx=cutoff;


Mp = sum(dcmp)*rho*disc_area*(omega*rad)^2*rad;
for j=1:360
    for i=1:length(inplane)
        inplane2(i,j)=inplane(i)*(omega*rad)^2;
        outplane2(i,j)=outplane(i)*(omega*rad)^2;
    end
end

%%
% vi = vtip*(lambdai);   % inflow speed
theta = theta*180/pi; % collective/pitch
PM = sum(pm)*rho*disc_area*(omega*rad)^2*rad;


%% encapsulation of output parameters

op.thrust_lbs = Thrust_lbs;
op.power_hp = Power_hp;
op.torque_Nm = Torque_Nm;
op.FM = FM;
op.eta = eta;
op.Mtip = Mtip;
op.rpm = rpm;
op.ct = ct;
op.cp = cp;
op.ctbysigma = ctbysigma;
op.freq = freq;
op.AR = AR;
op.Mp = Mp;
op.itn = iter;
op.kk = kk;


op.chd = chd;
op.pitch = theta_deg;
op.aoa = alpha_deg;
op.cl = cl;
op.vi = vi;		% inflow in meter/sec

% sectional force
op.inplane = inplane2;
op.outplane = outplane2;
op.dt = dt;

op.error = mat_error;
op.check = check;

end
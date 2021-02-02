function [op] =fn_bemt_isolated_ff(ip)
%% decoding parameters
% Tf = ip.Tf;
theta75_deg = ip.theta;
thetac_deg = ip.thetac;
thetas_deg = ip.thetas;
tilt_deg   = ip.tilt;

Nb			= ip.Nb;
rad_ft		= ip.rad;
rad_in_ft	= ip.rad_in;
sigma_eq	= ip.sigma;
thetatw_deg	= ip.twist;
taper		= ip.taper;
% sweep1_deg	= ip.sweep1;
% sweep2_deg	= ip.sweep2;
% rs			= ip.rs;

vc_knots	= ip.vc;
vinf_knots	= ip.vinf;
rho			= ip.rho; 

tc			= ip.tc;
tw1			= ip.tw1;

data		= ip.data;
% soln		= ip.soln; 
% use_soln	= ip.use_soln;
% tiploss		= ip.tiploss;
itermax		= ip.itermax_bemt;
iter_cutoff	= ip.cutoff_bemt;
% rw			= ip.rw;

delpsi	= ip.delpsi;
kp_ff	= ip.kp_ff;
kt_ff	= ip.kt_ff;
%% unit conversion

theta75 = theta75_deg*(pi/180);		% radian
thetac	= thetac_deg*(pi/180);		% radian
thetas	= thetas_deg*(pi/180);		% radian
tilt	= tilt_deg*(pi/180);		% radian
thetatw = thetatw_deg*(pi/180);		% radian
rad		= rad_ft*12*0.0254;			% meter
rad_in	= rad_in_ft*12*0.0254;		% meter
vc		= vc_knots*0.514444;		% meter/sec
v_inf	= vinf_knots*0.514444;		% meter/sec
% sweep1	= sweep1_deg*(pi/180);		% radian
% sweep2	= sweep2_deg*(pi/180);		% radian

%% rpm & some parameters

disc_area	= pi*rad^2;		% disc area
rpm			= ip.rpm;
omega		= rpm*2*pi/60;	% rotational speed
vtip		= omega*rad;    % tip speed (without considering fluid speed)
Mtip		= vtip/343;		% Tip Mach number
freq		= rpm/60*Nb;	% blade passage frequency
lambdac		= vc/vtip;      % non-dimensional climb velocity
mu = v_inf/vtip;
mux = v_inf*cos(tilt)/vtip;
muz = v_inf*sin(tilt)/vtip;
% ctf = ( Tf/(rho*disc_area*vtip^2) );
% ct_by_sigmaf = ctf/sigma_eq;


%% BEMT

[N,r,dr] = dicretization();
root_cut = round(rad_in/rad*N);
if (root_cut == 0)
    root_cut = 1;
end
% if (cutoff > 1)
% 	rc = (r(cutoff)+r(cutoff-1))*0.5;
% end
% if (cutoff == 1)
%     rc = 0;
% end

%% calculate sectional configurations   

thetatw2= thetatw/(1-r(root_cut));
chd_eq	= sigma_eq*disc_area/Nb/(rad);									% effective chord              
chd0	= chd_eq/ (1+ (1/taper-1)*(0.75-r(root_cut))/ (1.0-r(root_cut)) );	% root chord
chdtw	= chd0*(1/taper-1)/(1-r(root_cut));								% change in chord
chd1	= chd0+((1-r(root_cut)))*chdtw;									% 
AR		= (rad-rad_in)/(chd0+chd1)*2;									% aspect ratio

thetacol	= zeros(N,1);
chd		= zeros(N,1);
% sigma	= zeros(N,1);

for i = root_cut:N	
	chd(i,1)	= chd0 + chdtw*(r(i)-r(root_cut));			% sectional chord
% 	sigma(i,1)	= Nb*chd(i)*rad/disc_area;					% sectional solidity
	thetacol(i,1)	= theta75 + thetatw2 * (r(i)-0.75);			% sectional pitch
	if(r(i) <= tc)
		thetacol(i,1) = theta75+ thetatw2 * (tc-0.75)+ tw1*thetatw2 *(r(i)-tc);
	end
end

% calculate tangential distance of blade elements
% rd1 = rs*tan(sweep1);
% rd	= zeros(N,1);
% for i=1:N
% 	if(r(i) <= rs)
% 		rd(i,1) = r(i)*tan(sweep1);
% 	end
% 	if(r(i) > rs)
% 		rd(i,1) = rd1-(r(i)-rs)*tan(sweep2);
% 	end
% end

%% BEMT Loop

% azimuthal discretization

psi		= 0:delpsi:(360-delpsi);
psi		= psi'*pi/180;
npsi	= length(psi);

% matrix initialization

F			= zeros(N,npsi);
theta		= zeros(N,npsi);
lambdai		= zeros(N,npsi);
lambdai_old	= zeros(N,npsi);
Ut			= zeros(N,npsi);
Up			= zeros(N,npsi);
U			= zeros(N,npsi);
phi			= zeros(N,npsi);
alpha		= zeros(N,npsi);
cl2d		= zeros(N,npsi);
cd			= zeros(N,npsi);
cm			= zeros(N,npsi);
cl			= zeros(N,npsi);
dL			= zeros(N,npsi);
dD			= zeros(N,npsi);        
dFz			= zeros(N,npsi); 
dFtan		= zeros(N,npsi);
dFtan2		= zeros(N,npsi);
dQ			= zeros(N,npsi);
dQ2			= zeros(N,npsi);
dMy			= zeros(N,npsi);
dMx			= zeros(N,npsi);
dMp			= zeros(N,npsi);	

Ftan		= zeros(npsi,1);
Fx			= zeros(npsi,1);
Fy			= zeros(npsi,1);
Fz			= zeros(npsi,1);
Q			= zeros(npsi,1);
Q2			= zeros(npsi,1);
Mx			= zeros(npsi,1);
My			= zeros(npsi,1);
Mp			= zeros(npsi,1);

for j=1:length(psi)
    for i = 1:N
        F(i,j) = 1;
    end
end

lambda0 = 0;
lambda  = lambdac;

for iter=1:itermax    
	kai = atan2 (mux,(muz+lambda0) );
	kai3 = atan2 (mux,(muz) );
	kx = 15*pi/23*tan( kai/2 );
% 	kx = 0.5;
%   kx(i,j) = (sin(kai(i,j)))^2;
%   kx(i,j) = 4/3*(1-cos(kai(i,j))-1.8*mu^2)/sin(kai(i,j));
	ky = -2*mu*0;   
	for j =1:length(psi)        
		for i = root_cut:N
			theta(i,j) = thetacol(i) + thetac*cos(psi(j))+thetas*sin(psi(j));
			% dynamic inflow calculations
			lambdai_old(i,j) = lambdai(i,j);     
			lambdai(i,j) = lambda0*(1 + kx*r(i)*cos(psi(j)) + ky*r(i)*sin(psi(j)) );
			Ut(i,j) = ( r(i) + sin(psi(j))*mux )*vtip;
			Up(i,j) = ( lambdai(i,j) + lambdac + muz )*vtip;
			U(i,j) = sqrt( (Up(i,j))^2 + (Ut(i,j))^2 );
			phi(i,j) =atan2(Up(i,j),Ut(i,j));
			alpha(i,j) = theta(i,j)-phi(i,j);
% 			cl(i,j) = cl_alpha*alpha(i,j);
% 			cd(i,j) = cd0 + d1*alpha(i,j) + d2*alpha(i,j)^2;
% 			if (alpha(i,j)<0)
% 				alpha(i,j)=0;
% 			end
			[cl2d(i,j),cd(i,j),cm(i,j)] = fn_clcdcm(alpha(i,j),data);
			cl(i,j) = cl2d(i,j)/(1+1/AR);
			dL(i,j) = 0.5*rho*U(i,j)^2*chd(i)*cl(i,j)*dr*rad;
			dD(i,j) = 0.5*rho*U(i,j)^2*chd(i)*cd(i,j)*dr*rad;        
			dFz(i,j) = dL(i,j)*cos(phi(i,j)) - dD(i,j)*sin(phi(i,j));  
			dFtan(i,j) = dL(i,j)*sin(phi(i,j)) + dD(i,j)*cos(phi(i,j));
			dFtan2(i,j) = kp_ff*dL(i,j)*sin(phi(i,j)) + dD(i,j)*cos(phi(i,j));
			dQ(i,j) = -dFtan(i,j)*r(i)*rad;
			dQ2(i,j) = -dFtan2(i,j)*r(i)*rad;
			dMy(i,j) = -dFz(i,j)*r(i)*cos(psi(j))*rad;
			dMx(i,j) = dFz(i,j)*r(i)*sin(psi(j))*rad;
			dMp(i,j) = 0.5*rho*U(i,j)^2*chd(i)^2*cm(i,j)*dr*rad;        
% 			f(i,j) = Nb*(1-r(i))/lambda(i,j)/2;
% 			if (tiploss ==1)
% 				F(i,j) = 2*acos( exp(-f(i)) )/pi;
% 			end
		end    
		Ftan(j,1) = sum(dFtan(:,j));
		Fx(j,1) = sum(dFtan(:,j))*sin(psi(j));
		Fy(j,1) = -sum(dFtan(:,j))*cos(psi(j));
		Fz(j,1) = sum(dFz(:,j));
		Q(j,1) = sum(dQ(:,j)); 
		Q2(j,1) = sum(dQ2(:,j)); 
		Mx(j,1) = sum(dMx(:,j));
		My(j,1) = sum(dMy(:,j));
		Mp(j,1) = sum(dMp(:,j));
	end
    ct = Nb*mean(Fz)/(rho*disc_area*vtip^2);
	
% 	steady part of dynamic inflow update    
% 	lambda0 = 0.5*ct/sqrt( lambda0^2 + mu^2);
	ct_m = 1*ct;
    lambda_old = lambda;
	if (iter == 1)
		lambda = sqrt(0.5*( sqrt(mu^4+ ct_m^2)-mu^2 ) ) + lambdac;   
	end
	if (iter>1)
		lambda = mu*tan(tilt)*0+0.5*ct/sqrt(mu^2+lambda^2) + lambdac;
	end
    
    lambda = (lambda+lambda_old)/2;
	lambda0 = lambda-lambdac;
	
    error = lambdai-lambdai_old;
    norm_error = 0;
    for j =1:length(psi)
        for i=root_cut:N
            norm_error = norm_error+error(i,j)^2;
        end
    end
    error2(iter,1)=norm_error;
    
    if((error2(iter)<iter_cutoff) && (iter >1))
        break
    end
end



% kk=inducedk(theta0,R, R_in,sigma,area,Nb,thetatw*180/pi,taper,data);


Fx = kt_ff*Fx;
Fy = kt_ff*Fy;
Fz  = kt_ff*Fz;

net_Mx = Nb*mean(Mx);
net_My = Nb*mean(My);
net_Mz = Nb*mean(Q2);
net_Mp = mean(Mp);
Forcex = Nb*mean(Fx);
Forcey = Nb*mean(Fy);
Forcez = Nb*mean(Fz);
torque = -Nb*mean(Q);
torque2 = -Nb*mean(Q2);
power = omega*torque;
power2 = omega*torque2;
power = abs(power);
power2 = abs(power2);
Force_total = sqrt((Nb*mean(Fz))^2+(Nb*mean(Fx))^2);

thrust_N = Force_total;
power_watt = power2;
torque_Nm = torque2;

thrust_lbs = thrust_N*0.224809;
power_hp = power_watt/746;

theta_deg = theta.*180/pi;		% sectional pitch in deg.
alpha_deg = alpha.*180/pi;		% sectional AOA

phi_deg = phi*180/pi;
alpha2 = atan(Forcex/Forcez);
alphaf2 = tilt - alpha2;

load(1) = Forcex;	% Fx
load(2) = Forcey;	% Fy
load(3) = Forcez;	% Fz
load(4) = net_Mx;	% Mx
load(5) = net_My;	% My
load(6) = net_Mz;	% Mz

%% encapsulation of output parameters

op.thrust_lbs = thrust_lbs;	
op.power_hp = power_hp;
op.torque_Nm = torque_Nm;
op.Mtip = Mtip;
op.rpm = rpm;
% op.ct = ct;
% op.cp = cp;
op.ctbysigma = ct/sigma_eq;
op.freq = freq;
op.AR = AR;
% op.Mp = Mp;
op.iter = iter;

op.chd = chd;
op.pitch = theta_deg;
op.aoa = alpha_deg;
op.cl = cl;
% op.vi = vi;		% inflow in meter/sec
% 
% % sectional force
% op.inplane = inplane2;
% op.outplane = outplane2;
% op.dt = dt;
% 
% op.error = mat_error;
% op.check = check;

op.mean_inflow = lambda0*vtip;
op.load = load;

end
function [ip]=fn_input_modifier(ip)
vinf = ip.vinf*0.514444; % forward speed in meter
rad = ip.rad*12*0.0254;	% radius in m
fp = ip.fp*(0.3048^2);	% flat plate area in m^2
disc_area = pi*(rad)^2;

%%
yct = ip.yct;
if (yct == 0)
	rpm = ip.rpm;
end
if (yct == 1)
	ctbysigma = ip.ctbysigma;
	rpm = (60/(2*pi*rad))*sqrt(ip.treq*4.44822/(ctbysigma*ip.sigma*ip.rho*disc_area*2));
	ip.rpm = rpm;
end

%%
weight_N = ip.treq*4.44822;
drag_N = 0.5*ip.rho*vinf^2*fp;
total_force_N = sqrt( weight_N^2 + drag_N^2 );
alphaf = atan (drag_N/weight_N);
tilt = alphaf*180/pi;

%%

vtip = (rpm*2*pi/60)*(rad);
ctf = ( total_force_N/(2*ip.rho*disc_area*vtip^2) );
mu = vinf/vtip;
lambda_m =  sqrt(0.5*( sqrt(mu^4+ ctf^2)-mu^2 ) );
mux = vinf*cos(alphaf)/vtip;
muz = vinf*sin(alphaf)/vtip;
kai = atan2 (mux,(muz+lambda_m) )*180/pi;

gamma_UL = [-0.19 46.6 -0.06853]*[1 mu kai]'; % 3-point
gamma_LU = [0.8211 8.94 -0.01952]*[1 mu kai]'; % 3-point

%%	

ip.total_drag_N =drag_N;
ip.gamma_UL = gamma_UL;
ip.gamma_LU = gamma_LU;
end
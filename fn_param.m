function [ip] = fn_param()
%% rotor configuration

	ip.Nb		= 4;				% number of blades
	ip.rad		= 4.1;				% rotor radius in ft
	ip.rad_in	= 0.34*ip.rad;		% root cut in ft
	ip.sigma	= 0.2;				% equivalent solidity
	ip.twist	= -9;				% blade twist in deg
	ip.taper	= 2;				% blade taper

% option for blade sweep
% for unswept blades, keep sweep1 = sweep2=0

	ip.sweep1	= 5;				% forward sweep in deg
	ip.sweep2	= 20;				% backward sweep in deg
	ip.rs		= 0.6209;			% point of sweep reversal
	
% option for dual twist, keep tw1 = 1 for uniform twist rate

	ip.tw1	= 1;					% inner twist (normalized by by thetatw)
	ip.tc	= 0.4;					% radial location where twist rate changes
	
%% flow conditions

	ip.vc	= 0;					% axial velocity in knots
	ip.vinf = 30/0.514444;			% forward velocity in knots, unused for hover
	ip.rho	= 1.2256;				% density kg/m^3;
	
%% control parameters

	ip.treq = 700;					% required thrust in lbs
	ip.theta = 8;					% collective pitch (at 75% span) in degree
	ip.thetac = 0;					% cosinusoidal pitch, unused for hover
	ip.thetas = 0;					% sinusoidal pitch, unused for hover
	ip.tilt = 0;					% tilt, unused for hover
	
	ip.yct		= 1;
	if (ip.yct == 0)
		ip.rpm	= 950;				% rpm
		ip.ctbysigmacheck = ip.treq*4.44822/(ip.rho*2*pi*(ip.rad*12*0.0254)^2*(ip.rad*12*0.0254*ip.rpm*2*pi/60)^2*ip.sigma);
	end
	if (ip.yct == 1)
		
		ip.ctbysigma = 0.11;		% blade loading
% 		ip.rpm = (2*pi/60)*sqrt(ip.th_req/(ip.ctbysigma*ip.sigma_eq*ip.rho*ip.rad))
	end	
	
	
%% airfoil data
	% load airfoil data : C_L and C_D as function of AOA for certain Reynolds and Mach number
	% Present data is generated using Xfoil
	% include more airfoil data if needed
% 	cd airfoil_data
	load airfoil_data\rc410.dat;
% 	load ch10.dat;
% 	load s1223.dat;
% 	load naca0012.dat;
% 	load cy.dat;
	
	ip.af = 1;
	if (ip.af == 1)
		ip.data = rc410;		% choose airfoil section
	end
	if(ip.af == 0)
		ip.cl_alpha = 2*pi;
		ip.cd = [0.01,0.5,1];
	end

%% initial solution, if required
% keep use_soln=1 if you want to use solution file, otherwise use any other number

	ip.use_soln = 0;	
	ip.soln = xlsread('init2.xlsx');	
	
%% some more parameters
% tiploss flag,tiploss =1 means effect of tiploss is considered, 
% tiploss = any other number to neglect tiploss effect
	ip.tiploss = 1;	
	ip.fp = 24*0.35;		% flat plate area in ft^2
	ip.cutoff_bemt = 1e-6;	% residual cutoff value, bemt loop
	ip.cutoff_coax = 1e-5;	% residual cutoff value, coax interference loop
	ip.cutoff_trim = 1e-4;	% residual cutoff value, trim loop
	ip.itermax_bemt	= 30;	% max. number of iterations, bemt loop
	ip.itermax_coax = 10;	% max. number of iterations, coax interference loop
	ip.itermax_trim = 10;	% max. number of iterations, trim loop
	ip.rw	= 1;			% denotes the radial location that sees climb vel., keep it 1
	ip.rot_sep = 0.25;	% vertical distance between coax rotors
	ip.gamma_UL = 0.6;		% coax interference parameters, keep it fixed
	ip.gamma_LU = 0.4;		% coax interference parameters, keep it fixed
	ip.delpsi = 10;		% delta azimuth (deg.) for ff, keep 5 to 10
	ip.kp_ff = 1.00;		% increase in power due to non-uniform inflow in ff
	ip.kt_ff = 1.00;		% reduction in thrust due to non-uniform inflow in ff
	
end
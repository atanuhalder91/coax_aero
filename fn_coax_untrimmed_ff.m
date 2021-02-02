function [op1,op2,op0] = fn_coax_untrimmed_ff(ip1,ip2)
%%
itermax	= ip1.itermax_coax;
iter_cutoff	= ip1.cutoff_coax;
tilt = ip1.tilt*pi/180;
total_drag_N = ip1.total_drag_N;
weight_N = ip1.treq*4.44822;

gamma_UL = ip1.gamma_UL;
gamma_LU = ip1.gamma_LU;

d=ip1.rot_sep;
k_UL = 1+(d/sqrt(1+d^2))^gamma_UL;
k_LU = 1-(d/sqrt(1+d^2))^gamma_LU;

%%
vc1 = 0;
vc2 = 0;

vc1_old = 0;
vc2_old = 0;

vc1_old2 = 0;
vc2_old2 = 0;

for iter=1:itermax
    ip1.vc = vc1/0.514444;
	ip2.vc = vc2/0.514444;
	
    [op1] = fn_bemt_isolated_ff(ip1);	
    [op2] = fn_bemt_isolated_ff(ip2);
    
    %% climb velocity update
    vc1_old2 = vc1_old;
    vc2_old2 = vc2_old;
    
    vc1_old = vc1;
    vc2_old = vc2;
    
    vc1 = op1.mean_inflow*k_LU;
    vc2 = op2.mean_inflow*k_UL;
    
    vc1 = (vc1+vc1_old+vc1_old2*0)/2;
    vc2 = (vc2+vc2_old+vc2_old2*0)/2;
    
    mat_vc1(iter,1) = vc1;
    mat_vc2(iter,1) = vc2;
    error_vc(iter,1) = sqrt( (vc1-vc1_old)^2+(vc2-vc2_old)^2 );
	
	if (error_vc(iter)<iter_cutoff)
		break
	end
    
end

load1 = op1.load;
load2 = op2.load;

%% trim correction

delta = zeros(7,1);

% delta(1) = -42.48;
% delta(2) = -20.09;
% delta(3) = 169.19;
% delta(4) = 132.7201;
% delta(5) = -78.4217;
% delta(6) = 447.9705;
% delta(7) = -5.8416;


%% get rotor force and moments

Fx1 = load1(1);
Fy1 = load1(2);
Fz1 = load1(3);
Mx1 = load1(4);
My1 = load1(5);
Mz1 = load1(6);

Fx2 =  load2(1);
Fy2 = -load2(2);
Fz2 =  load2(3);
Mx2 = -load2(4);
My2 =  load2(5);
Mz2 = -load2(6);

%% calculate residuals

residual(1,1) = (Fx1+Fx2) - weight_N*sin(tilt)+total_drag_N*cos(tilt) + delta(1);
residual(2,1) = Fy1 + Fy2 + delta(2);
residual(3,1) = (Fz1+Fz2) - weight_N*cos(tilt)-total_drag_N*sin(tilt) + delta(3);
residual(4,1) = Mx1 + delta(4);
residual(5,1) = Mx2 + delta(5);
residual(6,1) = My1 + My2 + delta(6) ;
residual(7,1) = Mz1 + Mz2 + delta(7);

%% encapsulate outputs

op0.residual = residual;
op0.iter = iter;
op0.error = error_vc;


	
end
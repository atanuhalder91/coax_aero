function[op1,op2,op0]=fn_coax_untrimmed(ip1,ip2)
%%
itermax	= ip1.itermax_coax;
iter_cutoff	= ip1.cutoff_coax;
treq = ip1.treq;		% required thrust in lbs
gamma_UL = ip1.gamma_UL;
gamma_LU = ip1.gamma_LU;

d=ip1.rot_sep;
k_UL = 1+(d/sqrt(1+d^2))^gamma_UL;
k_LU = 1-(d/sqrt(1+d^2))^gamma_LU;
rw = sqrt(1/k_UL);
% rwu = 1;
[N,r,dr]= dicretization();

%%
vc1=0;
vc2=0;
vc1_old=0;
vc2_old=0;
mat_vc1 = zeros(itermax,1);
mat_vc2 = zeros(itermax,1);
error_vc = zeros(itermax,1);


for iter=1:itermax
	ip1.vc = vc1/0.514444;		% climb speed in knots
	ip2.vc = vc2/0.514444;		% climb speed in knots
% 	ip1.rw = 1;
	ip2.rw = rw;
	[op1]=fn_bemt_isolated(ip1);
	[op2]=fn_bemt_isolated(ip2);
%  %% mean induced velocity calculation
	for j=1:N
		vi_h1b(j,1)=2*r(j)*dr*op1.vi(j,1);
		vi_h2b(j,1)=2*r(j)*dr*op2.vi(j,1);
	end
	viU_mean = sum(vi_h1b);
	viL_mean = sum(vi_h2b);
 %% climb velocity update
	vc1_old2 = vc1_old;
	vc2_old2 = vc2_old;    
	vc1_old = vc1;
	vc2_old = vc2;    
	vc1 = viL_mean*k_LU;
	vc2 = viU_mean*k_UL;    
%     vc1 = (vc1+vc1_old+vc1_old2*0)/2;
%     vc2 = (vc2+vc2_old+vc2_old2*0)/2;    
	mat_vc1(iter,1) = vc1;
	mat_vc2(iter,1) = vc2;
	error_vc(iter,1) = sqrt( (vc1-vc1_old)^2+(vc2-vc2_old)^2 );   
	if (error_vc(iter)<iter_cutoff)
		break
	end
end

residual(1,1) = op1.power_hp-op2.power_hp;
residual(2,1) = (treq-(op1.thrust_lbs+op2.thrust_lbs));

op0.matvc1 = mat_vc1;
op0.matvc2 = mat_vc2;
op0.errorvc = error_vc;
op0.residual = residual;
end
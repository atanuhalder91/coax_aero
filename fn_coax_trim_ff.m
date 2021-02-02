function[op1,op2,opt,op0]=fn_coax_trim_ff(ip1,ip2,init_ctrl)
%% subroutine to trim coax rotor
	
	itermax = ip1.itermax_trim;	% max number of iterations
	iter_cutoff = ip1.cutoff_trim;	% iteration stops when residual norm reaches this cutoff value
	
	% initial guess of control imputs
	ctrl = init_ctrl;
	
%% trim loop

	for iter=1:itermax
%         iter
		ip1.theta	= ctrl(1);
		ip1.thetac	= ctrl(2);
		ip1.thetas	= ctrl(3);
		ip2.theta	= ctrl(4);
		ip2.thetac	= ctrl(5);
		ip2.thetas	= ctrl(6);
		ip1.tilt	= ctrl(7);
		ip2.tilt	= ctrl(7);
		[op1,op2,op0] = fn_coax_untrimmed_ff(ip1,ip2);		% call coax rotor perofrmance subroutine
		residual = op0.residual;
    
		mat_ctrl(iter,:) = ctrl; 
		mat_info(iter,1) = op1.thrust_lbs;
		mat_info(iter,2) = op2.thrust_lbs;
		mat_info(iter,3) = op1.power_hp;
		mat_info(iter,4) = op2.power_hp;
		mat_info(iter,5) = op1.torque_Nm;
		mat_info(iter,6) = op2.torque_Nm;
		mat_info(iter,7) = norm(residual);
		
		norm_res(iter,1)=norm(residual);
		if (norm(residual)<iter_cutoff)
			break
		end
	
		[J] = fn_coax_jacob_ff(ip1,ip2,residual);
		del_resi = -J\residual;
		ctrl = ctrl + del_resi;
  
	end
	opt.total_thrust_lbs = op1.thrust_lbs + op2.thrust_lbs;
	opt.total_power_hp = op1.power_hp + op2.power_hp;
	opt.total_torque_Nm = op1.torque_Nm + op2.torque_Nm;
	opt.ctbysigma = 0.5*(op1.ctbysigma+op2.ctbysigma);
	opt.ctrl = ctrl;
	opt.ctrl_iter = mat_ctrl;
	opt.info = mat_info;
	opt.norm_res = norm_res;	
end
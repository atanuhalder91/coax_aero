function[op1,op2,opt,op0]=fn_coax_trim(ip1,ip2,init_ctrl)
%% subroutine to trim coax rotor
	
	itermax = ip1.itermax_trim;	% max number of iterations
	iter_cutoff = ip1.cutoff_trim;	% iteration stops when residual norm reaches this cutoff value
	
	theta1 = init_ctrl(1);			% initial guess of collective (upper rotor)
	theta2 = init_ctrl(2);			% initial guess of collective (lower rotor)
	
%% trim loop

	for iter=1:itermax
%         iter
		
		ip1.theta = theta1;
		ip2.theta = theta2;
		[op1,op2,op0]=fn_coax_untrimmed(ip1,ip2);		% call coax rotor perofrmance subroutine
		
% 		residual(1,1) = op1.power_hp-op2.power_hp;
% 		residual(2,1) = (treq-(op1.thrust_lbs+op2.thrust_lbs));
		residual = op0.residual;
    
		mat_theta_col(iter,1) = theta1; 
		mat_theta_col(iter,2) = theta2; 
		mat_info(iter,1) = op1.thrust_lbs;
		mat_info(iter,2) = op2.thrust_lbs;
		mat_info(iter,3) = op1.power_hp;
		mat_info(iter,4) = op2.power_hp;
		mat_info(iter,5) = op1.torque_Nm;
		mat_info(iter,6) = op2.torque_Nm;
		mat_info(iter,7) = norm(residual);
		if (norm(residual)<iter_cutoff)
			break
		end
	
		[J]=fn_coax_jacob(ip1,ip2,residual);			% copmutes jacobian
		del_resi = -J\residual;
		theta1 = theta1 + del_resi(1);
		theta2 = theta2 + del_resi(2);
  
	end
	opt.total_thrust_lbs = op1.thrust_lbs + op2.thrust_lbs;
	opt.total_power_hp = op1.power_hp + op2.power_hp;
	opt.total_torque_Nm = op1.torque_Nm + op2.torque_Nm;
	opt.ctbysigma = 0.5*(op1.ctbysigma+op2.ctbysigma);
	opt.pitch_iter = mat_theta_col;
	opt.pitch = [theta1, theta2];
	opt.info = mat_info;
	
	opt.power_ideal = (opt.total_thrust_lbs*4.44822)^(1.5)/sqrt(2*ip1.rho*pi*(ip1.rad*12*0.0254)^2)/746;
	opt.FM = opt.power_ideal/opt.total_power_hp; 
end
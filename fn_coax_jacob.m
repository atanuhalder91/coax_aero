function[J]=fn_coax_jacob(ip1,ip2,residual)
%% subroutine to compute jacobian

	ctrl = [ip1.theta,ip2.theta];
	dof = length(ctrl);
	delx = zeros(dof,1);
	mat_residual = zeros(dof,dof);
	J = zeros(dof,dof);

	for i=1:dof
		delx(i,1) = 0.1;
	end

for i=1:dof
	x = ctrl;
	x(i) = ctrl(i)+delx(i);
	
	ip1.theta = x(1);
	ip2.theta = x(2);

	[op1,op2,op0] = fn_coax_untrimmed(ip1,ip2);
	mat_residual(:,i) = op0.residual;
end

for i=1:dof
    for j=1:dof
        J(j,i) = (mat_residual(j,i)-residual(j))/delx(i);
    end
end
	

end
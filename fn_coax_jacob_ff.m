function[J] = fn_coax_jacob_ff(ip1,ip2,residual)

ctrl = [ip1.theta,ip1.thetac,ip1.thetas,ip2.theta,ip2.thetac,ip2.thetas,ip1.tilt]';

dof = length(ctrl);
delx = zeros(dof,1);
mat_residual = zeros(dof,dof);
J = zeros(dof,dof);

for i=1:dof
    delx(i,1) = 1;
end

for i=1:dof
	x = ctrl;
	x(i) = ctrl(i)+delx(i);
	
	ip1.theta = x(1);
	ip1.thetac = x(2);
	ip1.thetas = x(3);
	ip2.theta = x(4);
	ip2.thetac = x(5);
	ip2.thetas = x(6);
	ip1.tilt = x(7);
	ip2.tilt = x(7);

	[op1,op2,op0] = fn_coax_untrimmed_ff(ip1,ip2);
	mat_residual(:,i) = op0.residual;
end

for i=1:dof
    for j=1:dof
        J(j,i) = (mat_residual(j,i)-residual(j))/delx(i);
    end
end

end
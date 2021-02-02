clear all;
close all;
hold on
clc;

[ip] = fn_param();	% input parameters

var1 = 0.01:0.01:0.15;
% var1 = 0.08;
n1 = length(var1);

ctrl = [15, 15]; % initial guess of collectives for upper and lower rotor, in deg.

%% Calling coaxial trimmed function
% op1, op2 denote outputs related to individual rotor performance
% Here1 denotes upper rotor, 2 denotes lower rotor
% opt denotes output related to total rotor performance
%-----------------------------------------------------------------------
% initiate variables
% first column upper rotor, second column lower rotor,
% third colum is total rotor performance

thrust_lbs	= zeros(n1,3);			% thrust in lbs
power_hp	= zeros(n1,3);			% power in hp
torque_Nm	= zeros(n1,3);			% torque in Nm 
FM			= zeros(n1,3);			% figure of merit
ctbysigma	= zeros(n1,3);			% blade loading
pitch		= zeros(n1,2);			% collective pitch
rpm			= zeros(n1,1);			% rpm
Mtip		= zeros(n1,1);			% tip Mach number

for i=1:n1 
	i
	ip.ctbysigma = var1(i);	
	[op1,op2,opt,op0]=fn_coax_trim(ip,ip,ctrl);
	ctrl = opt.pitch;
	
	% extracting performance parameters
	
	thrust_lbs(i,1) = op1.thrust_lbs;
	thrust_lbs(i,2) = op2.thrust_lbs;
	thrust_lbs(i,3) = opt.total_thrust_lbs;
	
	power_hp(i,1) = op1.power_hp;
	power_hp(i,2) = op2.power_hp;
	power_hp(i,3) = opt.total_power_hp;
	
	torque_Nm(i,1) = op1.torque_Nm;
	torque_Nm(i,2) = op2.torque_Nm;
	torque_Nm(i,3) = opt.total_torque_Nm;
	
	FM(i,1) = op1.FM;
	FM(i,2) = op2.FM;
	FM(i,3) = opt.FM;
	
	ctbysigma (i,1) = op1.ctbysigma;
	ctbysigma (i,2) = op2.ctbysigma;
	ctbysigma (i,3) = opt.ctbysigma;
	
	pitch(i,1) = opt.pitch(1);
	pitch(i,2) = opt.pitch(2);
	
	rpm(i,1) = op1.rpm;
	Mtip(i,1) = op1.Mtip;	
	
end

%% saving results
% cd C:\Users\halde\Downloads\psuwopwop\VSP2WOPWOP-master
save('aero_op.mat')

%% plotting

fs = 16;

figure(1)
hold on; grid on; box on;
set(gca,'FontName','Times','FontSize',fs);
plot(var1,thrust_lbs(:,1),'-sr','LineWidth',2)
plot(var1,thrust_lbs(:,2),'-sb','LineWidth',2)
plot(var1,thrust_lbs(:,3),'-sk','LineWidth',2)
legend('Upper Rotor','Lower Rotor','Total')
ylabel('Thrust (lbs.) ')
xlabel('Blade Loading, C_T/\sigma')

figure(2)
hold on; grid on; box on;
set(gca,'FontName','Times','FontSize',fs);
plot(var1,power_hp(:,1),'-sr','LineWidth',2)
plot(var1,power_hp(:,2),'--sb','LineWidth',2)
plot(var1,power_hp(:,3),'-sk','LineWidth',2)
legend('Upper Rotor','Lower Rotor','Total')
ylabel('Power (hp) ')
xlabel('Blade Loading, C_T/\sigma')

figure(3)
hold on; grid on; box on;
set(gca,'FontName','Times','FontSize',fs);
plot(var1,torque_Nm(:,1),'-sr','LineWidth',2)
plot(var1,torque_Nm(:,2),'--sb','LineWidth',2)
plot(var1,torque_Nm(:,3),'-sk','LineWidth',2)
legend('Upper Rotor','Lower Rotor','Total')
ylabel('Torque(N-m) ')
xlabel('Blade Loading, C_T/\sigma')

figure(4)
hold on; grid on; box on;
set(gca,'FontName','Times','FontSize',fs);
plot(var1,FM(:,1),'-sr','LineWidth',2)
plot(var1,FM(:,2),'-sb','LineWidth',2)
plot(var1,FM(:,3),'-sk','LineWidth',2)
legend('Upper Rotor','Lower Rotor','Total')
ylabel('Figure of Merit ')
xlabel('Blade Loading, C_T/\sigma')

figure(5)
hold on; grid on; box on;
set(gca,'FontName','Times','FontSize',fs);
plot(var1,ctbysigma(:,1),'-sr','LineWidth',2)
plot(var1,ctbysigma(:,2),'-sb','LineWidth',2)
plot(var1,ctbysigma(:,3),'-sk','LineWidth',2)
legend('Upper Rotor','Lower Rotor','Total')
ylabel('Blade Loading, C_T/\sigma')
xlabel('Blade Loading, C_T/\sigma')

figure(6)
hold on; grid on; box on;
set(gca,'FontName','Times','FontSize',fs);
plot(var1,pitch(:,1),'-sr','LineWidth',2)
plot(var1,pitch(:,2),'-sb','LineWidth',2)
legend('Upper Rotor','Lower Rotor')
ylabel('Collective Pitch (deg.) ')
xlabel('Blade Loading, C_T/\sigma')

figure(7)
hold on; grid on; box on;
set(gca,'FontName','Times','FontSize',fs);
plot(var1,rpm,'-sr','LineWidth',2)
ylabel('RPM ')
xlabel('Blade Loading, C_T/\sigma')

figure(8)
hold on; grid on; box on;
set(gca,'FontName','Times','FontSize',fs);
plot(var1,Mtip,'-sr','LineWidth',2)
ylabel('Tip Mach Number ')
xlabel('Blade Loading, C_T/\sigma')



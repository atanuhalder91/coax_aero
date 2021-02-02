clear all;
close all;
hold on;
clc;

[ip] = fn_param();	% input parameters

% variable 1 is forward speed in knots, don't change it
% choose variable 2 (i.e. rpm,sigma,ctbysigma etc.)

var1 = 0:10:100;			
var2 = 0.07:0.01:0.08;

% var1 = 0;
% var2 = 0.08;

n1 = length(var1);
n2 = length(var2);

% initial value of control parameters, in deg.
init_ctrl = [15,0,0,15,0,0,0]';

%% Calling coaxial trimmed function
% op1, op2 denote outputs related to individual rotor performance
% Here1 denotes upper rotor, 2 denotes lower rotor
% opt denotes output related to total rotor performance
%-----------------------------------------------------------------------
% initiate variables
% first column upper rotor, second column lower rotor,
% third colum is total rotor performance
thrust_lbs	= zeros(n1,3,n2);			% thrust in lbs
power_hp	= zeros(n1,3,n2);			% power in hp
torque_Nm	= zeros(n1,3,n2);			% torque in Nm 
FM			= zeros(n1,3,n2);			% figure of merit
ctbysigma	= zeros(n1,3,n2);			% blade loading
pitch		= zeros(n1,2,n2);			% collective pitch
rpm			= zeros(n1,n2);			% rpm
Mtip		= zeros(n1,n2);			% tip Mach number
mat_ctrl	= zeros(n1,7,n2);			% tip Mach number
		
for j=1:n2
	j
	ctrl = init_ctrl;
% second variable is prescribed to rpm, change it to prescribe any other parameters
    ip.ctbysigma = var2(j);		
	for i=1:n1
		i
		ip.vinf = var1(i);
		[ip]=fn_input_modifier(ip);
		[op1,op2,opt,op0]=fn_coax_trim_ff(ip,ip,ctrl);
		ctrl = opt.ctrl;
		
		mat_ctrl(i,:,j) = ctrl; 
		thrust_lbs(i,1,j) = op1.thrust_lbs;
		thrust_lbs(i,2,j) = op2.thrust_lbs;
		thrust_lbs(i,3,j) = opt.total_thrust_lbs;
	
		power_hp(i,1,j) = op1.power_hp;
		power_hp(i,2,j) = op2.power_hp;
		power_hp(i,3,j) = opt.total_power_hp;
		
		torque_Nm(i,1,j) = op1.torque_Nm;
		torque_Nm(i,2,j) = op2.torque_Nm;
		torque_Nm(i,3,j) = opt.total_torque_Nm;
	
		ctbysigma (i,1,j) = op1.ctbysigma;
		ctbysigma (i,2,j) = op2.ctbysigma;
		ctbysigma (i,3,j) = opt.ctbysigma;
	
		rpm(i,j) = op1.rpm;
		Mtip(i,j) = op1.Mtip;
	end
end

%% saving results
% cd C:\Users\halde\Downloads\psuwopwop\VSP2WOPWOP-master
save('aero_op_ff.mat')

%% plotting

fs = 16;
cc1 = [0 100 0]*(1/255);	% custom color (cc) 1
j=1;

% plot control inputs

figure(1)
hold on; grid on; box on;
set(gca,'FontName','Times','FontSize',fs);
plot(var1,mat_ctrl(:,1,j),'-r','LineWidth',2)
plot(var1,mat_ctrl(:,4,j),'--sb','LineWidth',2)
legend('Upper Rotor','Lower Rotor')
ylabel('Collective Pitch (deg.) ')
xlabel('Forward Speed (knots)')

figure(2)
hold on; grid on; box on;
set(gca,'FontName','Times','FontSize',fs);
plot(var1,mat_ctrl(:,2,j),'-r','LineWidth',2)
plot(var1,mat_ctrl(:,5,j),'--sr','LineWidth',2)
plot(var1,mat_ctrl(:,3,j),'-b','LineWidth',2)
plot(var1,mat_ctrl(:,6,j),'--sb','LineWidth',2)
legend('\theta_{1c},Upper Rotor','\theta_{1c},Lower Rotor','\theta_{1s},Upper Rotor','\theta_{1s},Lower Rotor')
ylabel('Cyclic Pitch (deg.) ')
xlabel('Forward Speed (knots)')

figure(3)
hold on; grid on; box on;
set(gca,'FontName','Times','FontSize',fs);
plot(var1,mat_ctrl(:,7,j),'-sr','LineWidth',2)
ylabel('Shaft Tilt (deg.) ')
xlabel('Forward Speed (knots)')

% plot performance

figure(4)
hold on; grid on; box on;
set(gca,'FontName','Times','FontSize',fs);
plot(var1,thrust_lbs(:,1,j),'-sr','LineWidth',2)
plot(var1,thrust_lbs(:,2,j),'-sb','LineWidth',2)
plot(var1,thrust_lbs(:,3,j),'-sk','LineWidth',2)
legend('Upper Rotor','Lower Rotor','Total')
ylabel('Thrust (lbs.) ')
xlabel('Forward Speed (knots)')

figure(5)
hold on; grid on; box on;
set(gca,'FontName','Times','FontSize',fs);
plot(var1,power_hp(:,1,j),'-sr','LineWidth',2)
plot(var1,power_hp(:,2,j),'--sb','LineWidth',2)
plot(var1,power_hp(:,3,j),'-sk','LineWidth',2)
legend('Upper Rotor','Lower Rotor','Total')
ylabel('Power (hp) ')
xlabel('Forward Speed (knots)')
ylim([0 150])

figure(6)
hold on; grid on; box on;
set(gca,'FontName','Times','FontSize',fs);
plot(var1,torque_Nm(:,1,j),'-sr','LineWidth',2)
plot(var1,torque_Nm(:,2,j),'--sb','LineWidth',2)
plot(var1,torque_Nm(:,3,j),'-sk','LineWidth',2)
legend('Upper Rotor','Lower Rotor','Total')
ylabel('Torque(N-m) ')
xlabel('Forward Speed (knots)')

figure(7)
hold on; grid on; box on;
set(gca,'FontName','Times','FontSize',fs);
plot(var1,ctbysigma(:,1,j),'-sr','LineWidth',2)
plot(var1,ctbysigma(:,2,j),'-sb','LineWidth',2)
plot(var1,ctbysigma(:,3,j),'-sk','LineWidth',2)
legend('Upper Rotor','Lower Rotor','Total')
ylabel('Blade Loading, C_T/\sigma')
xlabel('Forward Speed (knots)')

%% plot variation with second variable
figure(8)
hold on; grid on; box on;
set(gca,'FontName','Times','FontSize',fs);
plot(var1,power_hp(:,3,1),'-sr','LineWidth',2)
plot(var1,power_hp(:,3,2),'-sb','LineWidth',2)
% plot(var1,power_hp(:,3,3),'-sk','LineWidth',2)
% plot(var1,power_hp(:,3,4),'-sm','LineWidth',2)
% plot(var1,power_hp(:,3,5),'color',cc1,'LineWidth',2)
ylabel('Total Power (hp.) ')
xlabel('Forward Speed (knots)')
ylim([0 150])

figure(9)
hold on; grid on; box on;
set(gca,'FontName','Times','FontSize',fs);
plot(var1,rpm(:,1),'-sr','LineWidth',2)
plot(var1,rpm(:,2),'-sb','LineWidth',2)
% plot(var1,rpm(:,3),'-sk','LineWidth',2)
% plot(var1,rpm(:,4),'-sm','LineWidth',2)
% plot(var1,rpm(:,5),'color',cc1,'LineWidth',2)
ylabel('RPM ')
xlabel('Forward Speed (knots)')

figure(10)
hold on; grid on; box on;
set(gca,'FontName','Times','FontSize',fs);
plot(var1,Mtip(:,1),'-sr','LineWidth',2)
plot(var1,Mtip(:,2),'-sb','LineWidth',2)
% plot(var1,Mtip(:,3),'-sk','LineWidth',2)
% plot(var1,Mtip(:,4),'-sm','LineWidth',2)
% plot(var1,Mtip(:,5),'color',cc1,'LineWidth',2)
ylabel('Tip Mach Number ')
xlabel('Forward Speed (knots)')

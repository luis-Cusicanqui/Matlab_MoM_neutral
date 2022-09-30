%Only plot after computing the error
clear all
close all

testname =  "data/chapter5_norm/convergence/disc/test_weak_disc_u_convergence_rcx_eq";
error_name = [testname+"_1e4",testname+"_1e3",testname+"_1e2",testname+"_1e1",testname+"_1e0"];
error_name = error_name + "_1000";
%%
aspect_ratio = [1 1 1];
position = [0,0,0.5,0.5];
namename = "_rcx_eq_1e4";
linewidth = 1.5;
fontsize_label=13;
set(0,'DefaultLegendFontSize',13);
set(0,'DefaultLineLineWidth',linewidth);
%legend(["10^0","10^1","10^2","10^3","10^4"],"Interpreter",'latex');
collisionrates = [10^4,10^3,10^2,10^1,10^0];
error_filename= error_name+ "_errors";
legende = ["HME","QBME","LIN"];

figurepath = '../../../latex/official2/figures/ch-experiments/models/disc/';
savefig = 0;

%%
load(error_filename(1));
hme_errors(:,1) = [hme_rho_error1;hme_mom_error1;hme_energy_error1];
qbme_errors(:,1) = [qbme_rho_error1;qbme_mom_error1;qbme_energy_error1];
lin_errors(:,1) = [lin_rho_error1;lin_mom_error1;lin_energy_error1];

load(error_filename(2));
hme_errors(:,2) = [hme_rho_error1;hme_mom_error1;hme_energy_error1];
qbme_errors(:,2) = [qbme_rho_error1;qbme_mom_error1;qbme_energy_error1];
lin_errors(:,2) = [lin_rho_error1;lin_mom_error1;lin_energy_error1];

load(error_filename(3));
hme_errors(:,3) = [hme_rho_error1;hme_mom_error1;hme_energy_error1];
qbme_errors(:,3) = [qbme_rho_error1;qbme_mom_error1;qbme_energy_error1];
lin_errors(:,3) = [lin_rho_error1;lin_mom_error1;lin_energy_error1];

load(error_filename(4));
hme_errors(:,4) = [hme_rho_error1;hme_mom_error1;hme_energy_error1];
qbme_errors(:,4) = [qbme_rho_error1;qbme_mom_error1;qbme_energy_error1];
lin_errors(:,4) = [lin_rho_error1;lin_mom_error1;lin_energy_error1];

load(error_filename(5));
hme_errors(:,5) = [hme_rho_error1;hme_mom_error1;hme_energy_error1];
qbme_errors(:,5) = [qbme_rho_error1;qbme_mom_error1;qbme_energy_error1];
lin_errors(:,5) = [lin_rho_error1;lin_mom_error1;lin_energy_error1];

%%
fig = figure(1);
i=1;
plot(collisionrates,hme_errors(i,:),'+-')
hold on
plot(collisionrates,qbme_errors(i,:),'s-')
plot(collisionrates,lin_errors(i,:),'o-')
set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
xlabel("$\mathrm{R}_{cx}$","Interpreter","latex");
ylabel("Relative error");
legend(legende,'Location','northwest')
fig.Units = 'normalized';
fig.Position = position;
ax1=gca;
ax1.PlotBoxAspectRatio = [1 1 1];
xlabel("$\mathrm{R}_{cx}$",'FontSize', fontsize_label,"Interpreter","latex");
ylabel("Relative error","Interpreter","latex",'FontSize', fontsize_label);
grid on
box on
if savefig
   exportgraphics(gca,fullfile(figurepath,'density_error.png'),'Resolution',300); 
end


fig = figure(2);
i=2;
plot(collisionrates,hme_errors(i,:),'+-')
hold on
plot(collisionrates,qbme_errors(i,:),'s-')
plot(collisionrates,lin_errors(i,:),'o-')
set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
legend(legende)
ylabel("Relative error");

xlabel("$\mathrm{R}_{cx}$","Interpreter","latex");

fig.Units = 'normalized';
fig.Position = position;
ax2=gca;
ax2.PlotBoxAspectRatio = [1 1 1];
xlabel("$\mathrm{R}_{cx}$",'FontSize', fontsize_label,"Interpreter","latex");
ylabel("Relative error","Interpreter","latex",'FontSize', fontsize_label);
grid on
box on
if savefig
   exportgraphics(gca,fullfile(figurepath,'momentum_error.png'),'Resolution',300); 
end


fig = figure(3);
i=3;
plot(collisionrates,hme_errors(i,:),'+-')
hold on
plot(collisionrates,qbme_errors(i,:),'s-')
plot(collisionrates,lin_errors(i,:),'o-')
set(gca, 'XScale', 'log')

% set(gca, 'YScale', 'log')
legend(legende,'Location','southwest')

fig.Units = 'normalized';
fig.Position = position;
ax3=gca;
ax3.PlotBoxAspectRatio = [1 1 1];
xlabel("$\mathrm{R}_{cx}$",'FontSize', fontsize_label,"Interpreter","latex");
ylabel("Relative error","Interpreter","latex",'FontSize', fontsize_label);
grid on
box on

if savefig
   exportgraphics(gca,fullfile(figurepath,'energy_error.png'),'Resolution',300); 
end

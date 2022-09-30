%Only plot after computing the error
clear all
% close all

testname =  "data/models/disc/test_weak_disc_rho_1000_models_rcx_eq";
error_name = [testname+"_1e3",testname+"_1e2",testname+"_1e1",testname+"_1e0"];
%%
aspect_ratio = [1 1 1];
position = [0,0,0.5,0.5];
namename = "presentation_";
linewidth = 1.5;
fontsize_label=13;
set(0,'DefaultLegendFontSize',13);
set(0,'DefaultLineLineWidth',linewidth);
%legend(["10^0","10^1","10^2","10^3","10^4"],"Interpreter",'latex');
collisionrates = [10^3,10^2,10^1,10^0];
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
hme_errors(:,2) = [hme_rho_error2;hme_mom_error2;hme_energy_error2];
qbme_errors(:,2) = [qbme_rho_error2;qbme_mom_error2;qbme_energy_error2];
lin_errors(:,2) = [lin_rho_error2;lin_mom_error2;lin_energy_error2];

load(error_filename(3));
hme_errors(:,3) = [hme_rho_error3;hme_mom_error3;hme_energy_error3];
qbme_errors(:,3) = [qbme_rho_error3;qbme_mom_error3;qbme_energy_error3];
lin_errors(:,3) = [lin_rho_error3;lin_mom_error3;lin_energy_error3];

load(error_filename(4));
hme_errors(:,4) = [hme_rho_error4;hme_mom_error4;hme_energy_error4];
qbme_errors(:,4) = [qbme_rho_error4;qbme_mom_error4;qbme_energy_error4];
lin_errors(:,4) = [lin_rho_error4;lin_mom_error4;lin_energy_error4];


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
legend(legende,'Location','northeast')
fig.Units = 'normalized';
fig.Position = position;
ax1=gca;
ax1.PlotBoxAspectRatio = [1 1 1];
xlabel("$\mathrm{R}_{cx}$",'FontSize', fontsize_label,"Interpreter","latex");
ylabel("Relative error","Interpreter","latex",'FontSize', fontsize_label);
grid on
box on
if savefig
   exportgraphics(gca,fullfile(figurepath,namename+'density_error.png'),'Resolution',500); 
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
   exportgraphics(gca,fullfile(figurepath,namename+'momentum_error.png'),'Resolution',500); 
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
   exportgraphics(gca,fullfile(figurepath,namename+'energy_error.png'),'Resolution',500); 
end

%Only plot after computing the error
clear all
close all
rcx = 0;
nameplot = "rho";
testname =  "data/moments/disc/test_weak_disc_"+nameplot+"_1000_mom_rcx_1e"+rcx+"_mom_eq";
error_name1 = [testname+"_4",testname+"_5",testname+"_6",testname+"_7",testname+"_8"];
error_name2 = [testname+"_9",testname+"_10",testname+"_11",testname+"_12",testname+"_13"];
% error_name2 = [testname+"_14",testname+"_15",testname+"_16",testname+"_17",testname+"_18"];

%% Figure properties
aspect_ratio = [1 1.7 1];
position = [0,0,0.5,0.5];
namename = "disc_"+nameplot+"_rcx_eq_1e"+rcx+"_";
linewidth = 1.5;
fontsize_label=13;
set(0,'DefaultLegendFontSize',13);
set(0,'DefaultLineLineWidth',linewidth);
%legend(["10^0","10^1","10^2","10^3","10^4"],"Interpreter",'latex');
collisionrates = [4,5,6,7,8]%,9,10,11,12,13];

error_filename= [error_name1+ "_errors",error_name2+ "_errors"];
legende = ["HME","QBME","LIN"];

figurepath = '../../../latex/official2/figures/ch-experiments/moments/disc/';
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
load(error_filename(5));
hme_errors(:,5) = [hme_rho_error5;hme_mom_error5;hme_energy_error5];
qbme_errors(:,5) = [qbme_rho_error5;qbme_mom_error5;qbme_energy_error5];
lin_errors(:,5) = [lin_rho_error5;lin_mom_error5;lin_energy_error5];
%%
% load(error_filename(6));
% hme_errors(:,6) = [hme_rho_error1;hme_mom_error1;hme_energy_error1];
% qbme_errors(:,6) = [qbme_rho_error1;qbme_mom_error1;qbme_energy_error1];
% lin_errors(:,6) = [lin_rho_error1;lin_mom_error1;lin_energy_error1];
% 
% i=7;
% load(error_filename(i));
% hme_errors(:,i) = [hme_rho_error2;hme_mom_error2;hme_energy_error2];
% qbme_errors(:,i) = [qbme_rho_error2;qbme_mom_error2;qbme_energy_error2];
% lin_errors(:,i) = [lin_rho_error2;lin_mom_error2;lin_energy_error2];
% 
% load(error_filename(8));
% hme_errors(:,8) = [hme_rho_error3;hme_mom_error3;hme_energy_error3];
% qbme_errors(:,8) = [qbme_rho_error3;qbme_mom_error3;qbme_energy_error3];
% lin_errors(:,8) = [lin_rho_error3;lin_mom_error3;lin_energy_error3];
% 
% load(error_filename(9));
% hme_errors(:,9) = [hme_rho_error4;hme_mom_error4;hme_energy_error4];
% qbme_errors(:,9) = [qbme_rho_error4;qbme_mom_error4;qbme_energy_error4];
% lin_errors(:,9) = [lin_rho_error4;lin_mom_error4;lin_energy_error4];
% load(error_filename(10));
% hme_errors(:,10) = [hme_rho_error5;hme_mom_error5;hme_energy_error5];
% qbme_errors(:,10) = [qbme_rho_error5;qbme_mom_error5;qbme_energy_error5];
% lin_errors(:,10) = [lin_rho_error5;lin_mom_error5;lin_energy_error5];

%%


fig = figure(1);
i=1;
plot(collisionrates,hme_errors(i,:),'+-')
hold on
plot(collisionrates,qbme_errors(i,:),'s-')
plot(collisionrates,lin_errors(i,:),'o-')
% set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel("$\mathrm{R}_{cx}$","Interpreter","latex");
ylabel("Relative error");
legend(legende,'Location','northwest')
fig.Units = 'normalized';
fig.Position = position;
ax1=gca;
ax1.PlotBoxAspectRatio = aspect_ratio;
xlabel("$\mathrm{M}$",'FontSize', fontsize_label,"Interpreter","latex");
ylabel("Relative error","Interpreter","latex",'FontSize', fontsize_label);
xticks(collisionrates);
xlim([min(collisionrates),max(collisionrates)]);

grid on
box on
if savefig
   exportgraphics(gca,fullfile(figurepath,namename+'density_error.png'),'Resolution',300); 
end


fig = figure(2);
i=2;
plot(collisionrates,hme_errors(i,:),'+-')
hold on
plot(collisionrates,qbme_errors(i,:),'s-')
plot(collisionrates,lin_errors(i,:),'o-')
% set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
legend(legende,'Location','northwest')
ylabel("Relative error");

xlabel("$\mathrm{R}_{cx}$","Interpreter","latex");

fig.Units = 'normalized';
fig.Position = position;
ax2=gca;
ax2.PlotBoxAspectRatio = aspect_ratio;
xlabel("$\mathrm{M}$",'FontSize', fontsize_label,"Interpreter","latex");
ylabel("Relative error","Interpreter","latex",'FontSize', fontsize_label);
xticks(collisionrates);
xlim([min(collisionrates),max(collisionrates)]);

grid on
box on
if savefig
   exportgraphics(gca,fullfile(figurepath,namename+'momentum_error.png'),'Resolution',300); 
end


fig = figure(3);
i=3;
plot(collisionrates,hme_errors(i,:),'+-')
hold on
plot(collisionrates,qbme_errors(i,:),'s-')
plot(collisionrates,lin_errors(i,:),'o-')
% set(gca, 'XScale', 'log')

set(gca, 'YScale', 'log')
legend(legende,'Location','northwest')

fig.Units = 'normalized';
fig.Position = position;
ax3=gca;
ax3.PlotBoxAspectRatio = aspect_ratio;
xlabel("$\mathrm{M}$",'FontSize', fontsize_label,"Interpreter","latex");
ylabel("Relative error","Interpreter","latex",'FontSize', fontsize_label);
xticks(collisionrates);
xlim([min(collisionrates),max(collisionrates)]);
grid on
box on

if savefig
   exportgraphics(gca,fullfile(figurepath,namename+'energy_error.png'),'Resolution',300); 
end

close all
clear all
rcx = 4;
testname =  "data/moments/disc/test_weak_disc_theta_1000_mom_rcx_1e"+rcx+"_mom_eq";

plotname = [testname+"_4";testname+"_5";testname+"_6";testname+"_7";testname+"_8"];
plotname = [testname+"_9";testname+"_10";testname+"_11";testname+"_12";testname+"_13"];
% plotname = [testname+"_14";testname+"_15";testname+"_16";testname+"_17";testname+"_18"];

%% Get reference solution
plot_vince = 0;
filename_vince = "data/reference_solutions/sol_analogMC_FTS_m_testLuis_final_disrho"+rcx+"_L0_Smooth_Num_aNone_rno_BCper_N1e7_s1.h5";
x_vince = h5read(filename_vince,"/xval");
sol_vince = h5read(filename_vince,"/estimator");
if plot_vince
    figure(1) 
    plot(x_vince,sol_vince(:,1));
    figure(2) 
    plot(x_vince,sol_vince(:,2));
    figure(3) 
    plot(x_vince,sol_vince(:,3));
end
%% Figures

savefig = 0;
figurepath = '../../../latex/official2/figures/ch-experiments/moments/disc/';

plot_1 = 1;
plot_2 = 1;
plot_3 = 1;
plot_4 = 1;
plot_5 = 1;

plot_hme = 1;
plot_qbme= 0;
plot_lin = 0;

plot_density    =   0;
plot_momentum   =   0;
plot_energy     =   0;
plot_u          =   1;
plot_theta      =   1;

name4 = "un_hme.png";
name5 = "theta_hme.png";
legende = ["4","5","6","7","8"];%,"9","10","11","12","13",'Interpreter','latex'];
% legende = ["REF","9","10","11","12","13",'Interpreter','latex'];
% legende = ["14","15","16","17",'Interpreter','latex'];
% legende = ["12","13","14","15","16"]'
% legende = ["REF","HME","QBME","LIN"];
%% PLOT PROPERTIES
aspect_ratio = [1 1 1];
position = [0,0,0.5,0.5];
namename = "_rcx_eq_1e3";
linewidth = 1.5;
fontsize_label = 15;
set(0,'DefaultLineLineWidth',linewidth)
%%
if plot_hme

    if(plot_1)
        filename_hme = plotname(1)+"_hme.mat";
        load(filename_hme);
        [u_hme1,theta_hme1] = compute_macros(rho_hme_ave1,mom_hme_ave1,energy_hme_ave1);
        
        if plot_density
        figure(1)
        hold on
        plot(x,rho_hme_ave1)
        %title("Neutral mass density")
        end
        
        if plot_momentum
        figure(2)
        hold on
        plot(x,mom_hme_ave1)
        %title("Neutral particle momentum")
        end
        
        if plot_energy
        figure(3)
        hold on
        plot(x,energy_hme_ave1)
        %title("Neutral particle energy")
        end
        
        if plot_u
        figure(4)
        hold on
        plot(x,u_hme1);
        end
        
        if plot_theta
        figure(5)
        hold on
        plot(x,theta_hme1) 
        end
    end

    if(plot_2)
        filename_hme = plotname(2)+"_hme.mat";
        load(filename_hme);
        [u_hme2,theta_hme2] = compute_macros(rho_hme_ave2,mom_hme_ave2,energy_hme_ave2);
        
        if plot_density
        figure(1)
        hold on
        plot(x,rho_hme_ave2)
        %title("Neutral mass density")
        end
        
        if plot_momentum
        figure(2)
        hold on
        plot(x,mom_hme_ave2)
        %title("Neutral particle momentum")
        end
        
        if plot_energy
        figure(3)
        hold on
        plot(x,energy_hme_ave2)
        %title("Neutral particle energy")
        end
        
        if plot_u
        figure(4)
        hold on
        plot(x,u_hme2);
        end
        
        if plot_theta
        figure(5)
        hold on
        plot(x,theta_hme2) 
        end
    end
    if(plot_3)
        filename_hme = plotname(3)+"_hme.mat";
        load(filename_hme);
        [u_hme3,theta_hme3] = compute_macros(rho_hme_ave3,mom_hme_ave3,energy_hme_ave3);
        
        if plot_density
        figure(1)
        hold on
        plot(x,rho_hme_ave3)
        %title("Neutral mass density")
        end
        
        if plot_momentum
        figure(2)
        hold on
        plot(x,mom_hme_ave3)
        %title("Neutral particle momentum")
        end
        
        if plot_energy
        figure(3)
        hold on
        plot(x,energy_hme_ave3)
        %title("Neutral particle energy")
        end
        
        if plot_u
        figure(4)
        hold on
        plot(x,u_hme3);
        end
        
        if plot_theta
        figure(5)
        hold on
        plot(x,theta_hme3) 
        end
    end
    if(plot_4)
        filename_hme = plotname(4)+"_hme.mat";
        load(filename_hme);
        [u_hme4,theta_hme4] = compute_macros(rho_hme_ave4,mom_hme_ave4,energy_hme_ave4);
        
        if plot_density
        figure(1)
        hold on
        plot(x,rho_hme_ave4)
        %title("Neutral mass density")
        end
        
        if plot_momentum
        figure(2)
        hold on
        plot(x,mom_hme_ave4)
        %title("Neutral particle momentum")
        end
        
        if plot_energy
        figure(3)
        hold on
        plot(x,energy_hme_ave4)
        %title("Neutral particle energy")
        end
        
        if plot_u
        figure(4)
        hold on
        plot(x,u_hme4);
        end
        
        if plot_theta
        figure(5)
        hold on
        plot(x,theta_hme4) 
        end
    end
    if(plot_5)
        filename_hme = plotname(5)+"_hme.mat";
        load(filename_hme);
        [u_hme5,theta_hme5] = compute_macros(rho_hme_ave5,mom_hme_ave5,energy_hme_ave5);
        
        if plot_density
        figure(1)
        hold on
        plot(x,rho_hme_ave5)
        %title("Neutral mass density")
        end
        
        if plot_momentum
        figure(2)
        hold on
        plot(x,mom_hme_ave5)
        %title("Neutral particle momentum")
        end
        
        if plot_energy
        figure(3)
        hold on
        plot(x,energy_hme_ave5)
        %title("Neutral particle energy")
        end
        
        if plot_u
        figure(4)
        hold on
        plot(x,u_hme5);
        end
        
        if plot_theta
        figure(5)
        hold on
        plot(x,theta_hme5) 
        end
    end
end

%%
if plot_qbme

    if(plot_1)
        filename_qbme = plotname(1)+"_qbme.mat";
        load(filename_qbme);
        [u_qbme1,theta_qbme1] = compute_macros(rho_qbme_ave1,mom_qbme_ave1,energy_qbme_ave1);
        
        if plot_density
        figure(1)
        hold on
        plot(x,rho_qbme_ave1)
        %title("Neutral mass density")
        end
        
        if plot_momentum
        figure(2)
        hold on
        plot(x,mom_qbme_ave1)
        %title("Neutral particle momentum")
        end
        
        if plot_energy
        figure(3)
        hold on
        plot(x,energy_qbme_ave1)
        %title("Neutral particle energy")
        end
        
        if plot_u
        figure(4)
        hold on
        plot(x,u_qbme1);
        end
        
        if plot_theta
        figure(5)
        hold on
        plot(x,theta_qbme1) 
        end
    end

    if(plot_2)
        filename_qbme = plotname(2)+"_qbme.mat";
        load(filename_qbme);
        [u_qbme2,theta_qbme2] = compute_macros(rho_qbme_ave2,mom_qbme_ave2,energy_qbme_ave2);
        
        if plot_density
        figure(1)
        hold on
        plot(x,rho_qbme_ave2)
        %title("Neutral mass density")
        end
        
        if plot_momentum
        figure(2)
        hold on
        plot(x,mom_qbme_ave2)
        %title("Neutral particle momentum")
        end
        
        if plot_energy
        figure(3)
        hold on
        plot(x,energy_qbme_ave2)
        %title("Neutral particle energy")
        end
        
        if plot_u
        figure(4)
        hold on
        plot(x,u_qbme2);
        end
        
        if plot_theta
        figure(5)
        hold on
        plot(x,theta_qbme2) 
        end
    end
    if(plot_3)
        filename_qbme = plotname(3)+"_qbme.mat";
        load(filename_qbme);
        [u_qbme3,theta_qbme3] = compute_macros(rho_qbme_ave3,mom_qbme_ave3,energy_qbme_ave3);
        
        if plot_density
        figure(1)
        hold on
        plot(x,rho_qbme_ave3)
        %title("Neutral mass density")
        end
        
        if plot_momentum
        figure(2)
        hold on
        plot(x,mom_qbme_ave3)
        %title("Neutral particle momentum")
        end
        
        if plot_energy
        figure(3)
        hold on
        plot(x,energy_qbme_ave3)
        %title("Neutral particle energy")
        end
        
        if plot_u
        figure(4)
        hold on
        plot(x,u_qbme3);
        end
        
        if plot_theta
        figure(5)
        hold on
        plot(x,theta_qbme3) 
        end
    end
    if(plot_4)
        filename_qbme = plotname(4)+"_qbme.mat";
        load(filename_qbme);
        [u_qbme4,theta_qbme4] = compute_macros(rho_qbme_ave4,mom_qbme_ave4,energy_qbme_ave4);
        
        if plot_density
        figure(1)
        hold on
        plot(x,rho_qbme_ave4)
        %title("Neutral mass density")
        end
        
        if plot_momentum
        figure(2)
        hold on
        plot(x,mom_qbme_ave4)
        %title("Neutral particle momentum")
        end
        
        if plot_energy
        figure(3)
        hold on
        plot(x,energy_qbme_ave4)
        %title("Neutral particle energy")
        end
        
        if plot_u
        figure(4)
        hold on
        plot(x,u_qbme4);
        end
        
        if plot_theta
        figure(5)
        hold on
        plot(x,theta_qbme4) 
        end
    end
    if(plot_5)
        filename_qbme = plotname(5)+"_qbme.mat";
        load(filename_qbme);
        [u_qbme5,theta_qbme5] = compute_macros(rho_qbme_ave5,mom_qbme_ave5,energy_qbme_ave5);
        
        if plot_density
        figure(1)
        hold on
        plot(x,rho_qbme_ave5)
        %title("Neutral mass density")
        end
        
        if plot_momentum
        figure(2)
        hold on
        plot(x,mom_qbme_ave5)
        %title("Neutral particle momentum")
        end
        
        if plot_energy
        figure(3)
        hold on
        plot(x,energy_qbme_ave5)
        %title("Neutral particle energy")
        end
        
        if plot_u
        figure(4)
        hold on
        plot(x,u_qbme5);
        end
        
        if plot_theta
        figure(5)
        hold on
        plot(x,theta_qbme5) 
        end
    end
end

%%
if plot_lin

    if(plot_1)
        filename_lin = plotname(1)+"_lin.mat";
        load(filename_lin);
        [u_lin1,theta_lin1] = compute_macros(rho_lin_ave1,mom_lin_ave1,energy_lin_ave1);
        
        if plot_density
        figure(1)
        hold on
        plot(x,rho_lin_ave1)
        %title("Neutral mass density")
        end
        
        if plot_momentum
        figure(2)
        hold on
        plot(x,mom_lin_ave1)
        %title("Neutral particle momentum")
        end
        
        if plot_energy
        figure(3)
        hold on
        plot(x,energy_lin_ave1)
        %title("Neutral particle energy")
        end
        
        if plot_u
        figure(4)
        hold on
        plot(x,u_lin1);
        end
        
        if plot_theta
        figure(5)
        hold on
        plot(x,theta_lin1) 
        end
    end

    if(plot_2)
        filename_lin = plotname(2)+"_lin.mat";
        load(filename_lin);
        [u_lin2,theta_lin2] = compute_macros(rho_lin_ave2,mom_lin_ave2,energy_lin_ave2);
        
        if plot_density
        figure(1)
        hold on
        plot(x,rho_lin_ave2)
        %title("Neutral mass density")
        end
        
        if plot_momentum
        figure(2)
        hold on
        plot(x,mom_lin_ave2)
        %title("Neutral particle momentum")
        end
        
        if plot_energy
        figure(3)
        hold on
        plot(x,energy_lin_ave2)
        %title("Neutral particle energy")
        end
        
        if plot_u
        figure(4)
        hold on
        plot(x,u_lin2);
        end
        
        if plot_theta
        figure(5)
        hold on
        plot(x,theta_lin2) 
        end
    end
    if(plot_3)
        filename_lin = plotname(3)+"_lin.mat";
        load(filename_lin);
        [u_lin3,theta_lin3] = compute_macros(rho_lin_ave3,mom_lin_ave3,energy_lin_ave3);
        
        if plot_density
        figure(1)
        hold on
        plot(x,rho_lin_ave3)
        %title("Neutral mass density")
        end
        
        if plot_momentum
        figure(2)
        hold on
        plot(x,mom_lin_ave3)
        %title("Neutral particle momentum")
        end
        
        if plot_energy
        figure(3)
        hold on
        plot(x,energy_lin_ave3)
        %title("Neutral particle energy")
        end
        
        if plot_u
        figure(4)
        hold on
        plot(x,u_lin3);
        end
        
        if plot_theta
        figure(5)
        hold on
        plot(x,theta_lin3) 
        end
    end
    if(plot_4)
        filename_lin = plotname(4)+"_lin.mat";
        load(filename_lin);
        [u_lin4,theta_lin4] = compute_macros(rho_lin_ave4,mom_lin_ave4,energy_lin_ave4);
        
        if plot_density
        figure(1)
        hold on
        plot(x,rho_lin_ave4)
        %title("Neutral mass density")
        end
        
        if plot_momentum
        figure(2)
        hold on
        plot(x,mom_lin_ave4)
        %title("Neutral particle momentum")
        end
        
        if plot_energy
        figure(3)
        hold on
        plot(x,energy_lin_ave4)
        %title("Neutral particle energy")
        end
        
        if plot_u
        figure(4)
        hold on
        plot(x,u_lin4);
        end
        
        if plot_theta
        figure(5)
        hold on
        plot(x,theta_lin4) 
        end
    end
    if(plot_5)
        filename_lin = plotname(5)+"_lin.mat";
        load(filename_lin);
        [u_lin5,theta_lin5] = compute_macros(rho_lin_ave5,mom_lin_ave5,energy_lin_ave5);
        
        if plot_density
        figure(1)
        hold on
        plot(x,rho_lin_ave5)
        %title("Neutral mass density")
        end
        
        if plot_momentum
        figure(2)
        hold on
        plot(x,mom_lin_ave5)
        %title("Neutral particle momentum")
        end
        
        if plot_energy
        figure(3)
        hold on
        plot(x,energy_lin_ave5)
        %title("Neutral particle energy")
        end
        
        if plot_u
        figure(4)
        hold on
        plot(x,u_lin5);
        end
        
        if plot_theta
        figure(5)
        hold on
        plot(x,theta_lin5) 
        end
    end
end

%%
if plot_density
fig = figure(1);
fig.Units = 'normalized';

fig.Position = position;
ax1=gca;
ax1.PlotBoxAspectRatio = aspect_ratio;
grid on
box on
hline = findobj(gca, 'type', 'line');
% set(hline(1),'LineStyle','--');
% set(hline(2),'LineStyle',':');
% set(hline(3),'LineStyle','-');
xlabel("$x$",'FontSize', fontsize_label,"Interpreter","latex");
ylabel("$\rho_n$","Interpreter","latex",'FontSize', fontsize_label);
legend(legende,'Interpreter','latex');
if savefig
   exportgraphics(gca,fullfile(figurepath,'density'+namename+'.png'),'Resolution',300); 
end
end

if plot_momentum
fig = figure(2);
fig.Units = 'normalized';

fig.Position = position;
ax2=gca;
ax2.PlotBoxAspectRatio = aspect_ratio;
grid on
box on
hline = findobj(gca, 'type', 'line');
% set(hline(1),'LineStyle','--');
% set(hline(2),'LineStyle',':');
% set(hline(3),'LineStyle','-');
xlabel("$x$",'FontSize', fontsize_label,"Interpreter","latex");
ylabel("$m_n$","Interpreter","latex",'FontSize', fontsize_label);
legend(legende,'Interpreter','latex');
if savefig
   exportgraphics(gca,fullfile(figurepath,'momentum'+namename+'.png'),'Resolution',300); 
end
end

if plot_energy
fig = figure(3);
fig.Units = 'normalized';

fig.Position = position;ax3=gca;
ax3.PlotBoxAspectRatio = aspect_ratio;
grid on
box on
hline = findobj(gca, 'type', 'line');
% set(hline(1),'LineStyle','--');
% set(hline(2),'LineStyle',':');
% set(hline(3),'LineStyle','-');
xlabel("$x$",'FontSize', fontsize_label,"Interpreter","latex");
ylabel("$E_n$","Interpreter","latex",'FontSize', fontsize_label);
legend(legende,'Interpreter','latex');
if savefig
   exportgraphics(gca,fullfile(figurepath,'energy'+namename+'.png'),'Resolution',300); 
end
end

if plot_u
fig = figure(4);
fig.Units = 'normalized';

fig.Position = position;
ax2=gca;
ax2.PlotBoxAspectRatio = [1 1 1];
xlabel("$x$",'FontSize', fontsize_label,"Interpreter","latex");
ylabel("$u_n$","Interpreter","latex",'FontSize', fontsize_label);
grid on
box on
legend(legende,'Interpreter','latex');
if savefig
   exportgraphics(gca,fullfile(figurepath,name4),'Resolution',300); 
end
end

if plot_theta
fig = figure(5);
fig.Units = 'normalized';
fig.Position = position;

ax3=gca;
ax3.PlotBoxAspectRatio = [1 1 1];
grid on
box on
xlabel("$x$",'FontSize', fontsize_label,"Interpreter","latex");
ylabel("$\theta_n$","Interpreter","latex",'FontSize', fontsize_label);
legend(legende,'Location','southeast','Interpreter','latex');

if savefig
   exportgraphics(gca,fullfile(figurepath,name5),'Resolution',300); 
end
end


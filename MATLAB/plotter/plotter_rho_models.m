% close all
clear all


%Use to plot the results of the experiment with multiple collision rates.
%You can either show the results for multiple collision rates on top of
%each other or compare it to the reference solution.


testname =  "data/models/disc/test_weak_disc_rho_1000_models_rcx_eq";
% testname =  "data/models/smooth/test_smooth_rho_1000_models_rcx_eq";

plotname = [testname+"_1e3";testname+"_1e2";testname+"_1e1";testname+"_1e0"];
% 

%% Reference solution
plot_ref = 0; %only use if comparison with solution is necessary
collision = 4;
filename_ref = "data/reference_solutions/sol_analogMC_FTS_m_testLuis_final_disrho"+collision+"_L0_Smooth_Num_aNone_rno_BCper_N1e7_s1.h5";
x_ref = h5read(filename_ref,"/xval");
sol_ref = h5read(filename_ref,"/estimator");
if plot_ref
figure(1) 
plot(x_ref,sol_ref(:,1));
figure(2) 
plot(x_ref,sol_ref(:,2));
figure(3) 
plot(x_ref,sol_ref(:,3));
% legende = ["0.3","0.4","0.5","0.7","0.8"];
end
%% Figures
savefig =0;
figurepath = '../../../latex/official2/figures/ch-experiments/models/disc/';
plot_1 = 0;
plot_2 = 0;
plot_3 = 0;
plot_4 = 1;
plot_5 = 0;
plot_hme = 1;
plot_qbme= 1;
plot_lin = 1;

plot_density    =   1;
plot_momentum   =   1;
plot_energy     =   1;
plot_u          =   0;
plot_theta      =   0;

name4 = "un_hme.png";
name5 = "theta_hme.png";
% legende = ["$10^4$","$10^3$","$10^2$","$10^1$","$10^0$","Eq"];
% legende = ["REF","HME","QBME","LIN"];
% legende = ["HME","QBME","LIN"];
%% PLOT PROPERTIES
aspect_ratio = [1 1 1];
position = [0,0,0.5,0.5];
namename = "_rcx_eq_1e0";
linewidth = 1.5;
fontsize_label=13;
set(0,'DefaultLegendFontSize',13);
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
set(hline(1),'LineStyle','--');
set(hline(2),'LineStyle',':');
set(hline(3),'LineStyle','-');
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
set(hline(1),'LineStyle','--');
set(hline(2),'LineStyle',':');
set(hline(3),'LineStyle','-');
xlabel("$x$",'FontSize', fontsize_label,"Interpreter","latex");
ylabel("$\mathrm{m}_n$","Interpreter","latex",'FontSize', fontsize_label);
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
set(hline(1),'LineStyle','--');
set(hline(2),'LineStyle',':');
set(hline(3),'LineStyle','-');
xlabel("$x$",'FontSize', fontsize_label,"Interpreter","latex");
ylabel("$\mathrm{E}_n$","Interpreter","latex",'FontSize', fontsize_label);
legend(legende,'Interpreter','latex');
if savefig
   exportgraphics(gca,fullfile(figurepath,'energy'+namename+'.png'),'Resolution',300); 
end
end

if plot_u
fig = figure(4);
fig.Units = 'normalized';
% plot(x,u_p,'--');
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
% plot(x,theta_p,'--');

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

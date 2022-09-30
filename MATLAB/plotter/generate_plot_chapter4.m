pressure_hme4 = hme4(:,2).*hme4(:,4);
pressure_hme10 = hme10(:,2).*hme10(:,4);
pressure_qbme4 = qbme4(:,2).*qbme4(:,4);
pressure_qbme10 = qbme10(:,2).*qbme10(:,4);

close all
figurepath = '../../../latex/official2/figures/ch-softwareImpl/';



figure(1)
linewidth = 1.5;
set(0,'DefaultLineLineWidth',linewidth)
hold on
plot(x,U_sol1_hme(1,:,end),'-','Color','[0.21,0.61,0.21]')
plot(hme4(:,1),hme4(:,2),'--','Color','r');
ylim([0 8])
xlim([-1 1.5])
xticks(-1:0.5:1.5)
yticks(1:2:7)


ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';

grid on
% axis square
box on
legend(["$\textrm{Hme}$","$\textrm{Paper}$"],'Interpreter','latex')
saveas(gca,fullfile(figurepath,'density_hme4'),'png');
figure(2)
hold on
plot(x,pressure1_hme,'-','Color','[0.21,0.61,0.21]');
plot(x,pressure_hme4,'--','Color','r');
ylim([0 8])
xlim([-1 1.5])

xticks(-1:0.5:1.5)
yticks(1:2:7)

ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';

grid on
% axis square
box on
legend(["$\textrm{Hme}$","$\textrm{Paper}$"],'Interpreter','latex')
saveas(gca,fullfile(figurepath,'pressure_hme4'),'png');

figure(3)
hold on
plot(x,U_sol1_hme(2,:,end),'-','Color','[0.21,0.61,0.21]')
plot(hme4(:,1),hme4(:,3),'--','Color','r')
xticks(-1:0.5:1.5)
yticks(0:0.25:0.75);
ylim([-0.1250 0.8750])
xlim([-1 1.5])

ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';

grid on
% axis square
box on
% 
legend(["$\textrm{Hme}$","$\textrm{Paper}$"],'Interpreter','latex')

saveas(gca,fullfile(figurepath,'velocity_hme4'),'png');


%% 
figure(4)
linewidth = 1.5;
set(0,'DefaultLineLineWidth',linewidth)
hold on
plot(x,U_sol1_qbme(1,:,end),'-','Color','[0.21,0.61,0.21]')
plot(qbme4(:,1),qbme4(:,2),'--','Color','r');
ylim([0 8])
xlim([-1 1.5])
xticks(-1:0.5:1.5)
yticks(1:2:7)


ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';

grid on
% axis square
box on
legend(["$\textrm{Qbme}$","$\textrm{Paper}$"],'Interpreter','latex')
saveas(gca,fullfile(figurepath,'density_qbme4'),'png');

figure(5)
hold on
plot(x,pressure1_qbme,'-','Color','[0.21,0.61,0.21]');
plot(x,pressure_qbme4,'--','Color','r');
ylim([0 8])
xlim([-1 1.5])

xticks(-1:0.5:1.5)
yticks(1:2:7)

ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';

grid on
% axis square
box on
legend(["$\textrm{Qbme}$","$\textrm{Paper}$"],'Interpreter','latex')
saveas(gca,fullfile(figurepath,'pressure_qbme4'),'png');

figure(6)
hold on
plot(x,U_sol1_qbme(2,:,end),'-','Color','[0.21,0.61,0.21]')
plot(qbme4(:,1),qbme4(:,3),'--','Color','r')
xticks(-1:0.5:1.5)
yticks(0:0.25:0.75);
ylim([-0.1250 0.8750])
xlim([-1 1.5])

ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';

grid on
box on
legend(["$\textrm{Qbme}$","$\textrm{Paper}$"],'Interpreter','latex')

saveas(gca,fullfile(figurepath,'velocity_qbme4'),'png');

% 
legend(["$\textrm{Qbme}$","$\textrm{Paper}$"],'Interpreter','latex')

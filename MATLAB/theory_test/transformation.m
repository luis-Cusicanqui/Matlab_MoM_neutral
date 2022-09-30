syms x theta rho u
assume(rho,'real');
assume(theta,'real');
assume(theta,'positive');
assume(u,'real');
her = hermite(x,10);
f_max = 1/sqrt(2*sym(pi)*theta)*exp(-(x-u)^2/(2*theta));

u_trans = [2 -2 3]
theta_trans = [5 10 0.5]
f_max1 = subs(f_max,[u,theta],[u_trans(1),theta_trans(1)]);
f_max2 = subs(f_max,[u,theta],[u_trans(2),theta_trans(2)]);
f_max3 = subs(f_max,[u,theta],[u_trans(3),theta_trans(3)]);

f_1 = subs(f_max1,x,x*sqrt(theta_trans(1))+u_trans(1));
f_2 = subs(f_max2,x,x*sqrt(theta_trans(2))+u_trans(2));
f_3 = subs(f_max3,x,x*sqrt(theta_trans(3))+u_trans(3));

legende = ["$\textbf{M}_1$","$\textbf{M}_2$","$\textbf{M}_3$"];
figurepath = '../../../Latex/official2/figures/ch-math/transformed_bgk/';

figure
fplot(f_max1,[-10,10]);
hold on
fplot(f_max2,[-10,10]);
fplot(f_max3,[-10,10]);
xlabel('$v$','Interpreter','latex');
ylabel('$f(t,x,v)$','Interpreter','latex');
legend(legende,'Interpreter','latex');
grid on
exportgraphics(gca,fullfile(figurepath,'v_space.png'),'Resolution',300);
figure 
fplot(f_1)
hold on
fplot(f_2)
fplot(f_3)
xlabel('$\xi$','Interpreter','latex')
ylabel('$f(t,x,\xi)$','Interpreter','latex');
grid on
legend(legende,'Interpreter','latex');
exportgraphics(gca,fullfile(figurepath,'xi_space.png'),'Resolution',300);

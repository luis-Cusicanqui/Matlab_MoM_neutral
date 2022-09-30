% close all
clear all

%This file only computes the error of the test that end with models
collision_rate = 0; %can be 0 1 2 or 3
%Choose the location of the data of the test as defined on the test file
testname =  "data/moments/disc/test_weak_disc_rho_1000_mom_rcx_1e"+collision_rate+"_mom_eq"; 

%Define the collision rates tested. This should also match the
%corresponding collision rates as in the test file.
moments = [4 5 6 7 8];
% moments = [9 10 11 12 13];

plotname = [testname+"_"+moments(1); ...
    testname+"_"+moments(2); ...
    testname+"_"+moments(3); ...
    testname+"_"+moments(4); ...
    testname+"_"+moments(5)];

%CHOOSE ONE OF THE FOLLOWING 3
testname_ref = "data/reference_solutions/sol_analogMC_FTS_m_testLuis_final_disrho";
% testname_ref = "data/reference_solutions/sol_analogMC_FTS_m_testLuis_final_disu";
% testname_ref = "data/reference_solutions/sol_analogMC_FTS_m_testLuis_final_dise";

filename_ref = testname_ref+collision_rate+"_L0_Smooth_Num_aNone_rno_BCper_N1e7_s1.h5";
x_ref = h5read(filename_ref,"/xval");
sol_ref = h5read(filename_ref,"/estimator");
%% Variables to decide which error to compute
compute_error_1 = 1; % M= 4 or 9
compute_error_2 = 1; % M= 5 or 10
compute_error_3 = 1; % M= 6 or 11
compute_error_4 = 1; % M= 7 or 12
compute_error_5 = 1; % M= 8 or 13



compute_hme_error = 1;
compute_lin_error = 1;
compute_qbme_error = 1;

save_errors =1;

norm_number = 2;

%% Get the computed solution


if compute_hme_error
    if(compute_error_1)
        filename_hme = plotname(1)+"_hme.mat";
        load(filename_hme);
        [u_hme1,theta_hme1] = compute_macros(rho_hme_ave1,mom_hme_ave1,energy_hme_ave1);
    end

    if(compute_error_2)
        filename_hme = plotname(2)+"_hme.mat";
        load(filename_hme);
        [u_hme2,theta_hme2] = compute_macros(rho_hme_ave2,mom_hme_ave2,energy_hme_ave2);
 
    end
    if(compute_error_3)
        filename_hme = plotname(3)+"_hme.mat";
        load(filename_hme);
        [u_hme3,theta_hme3] = compute_macros(rho_hme_ave3,mom_hme_ave3,energy_hme_ave3);
    end
    if(compute_error_4)
        filename_hme = plotname(4)+"_hme.mat";
        load(filename_hme);
        [u_hme4,theta_hme4] = compute_macros(rho_hme_ave4,mom_hme_ave4,energy_hme_ave4);
 
    end
    if(compute_error_5)
        filename_hme = plotname(5)+"_hme.mat";
        load(filename_hme);
        [u_hme5,theta_hme5] = compute_macros(rho_hme_ave5,mom_hme_ave5,energy_hme_ave5);
 
    end
end

if compute_qbme_error

    if(compute_error_1)
        filename_qbme = plotname(1)+"_qbme.mat";
        load(filename_qbme);
        [u_qbme1,theta_qbme1] = compute_macros(rho_qbme_ave1,mom_qbme_ave1,energy_qbme_ave1);
        
    end

    if(compute_error_2)
        filename_qbme = plotname(2)+"_qbme.mat";
        load(filename_qbme);
        [u_qbme2,theta_qbme2] = compute_macros(rho_qbme_ave2,mom_qbme_ave2,energy_qbme_ave2);
    end
    if(compute_error_3)
        filename_qbme = plotname(3)+"_qbme.mat";
        load(filename_qbme);
        [u_qbme3,theta_qbme3] = compute_macros(rho_qbme_ave3,mom_qbme_ave3,energy_qbme_ave3);
        

    end
    if(compute_error_4)
        filename_qbme = plotname(4)+"_qbme.mat";
        load(filename_qbme);
        [u_qbme4,theta_qbme4] = compute_macros(rho_qbme_ave4,mom_qbme_ave4,energy_qbme_ave4);
    end
    if(compute_error_5)
        filename_qbme = plotname(5)+"_qbme.mat";
        load(filename_qbme);
        [u_qbme5,theta_qbme5] = compute_macros(rho_qbme_ave5,mom_qbme_ave5,energy_qbme_ave5);
    end
end

if compute_lin_error

    if(compute_error_1)
        filename_lin = plotname(1)+"_lin.mat";
        load(filename_lin);
        [u_lin1,theta_lin1] = compute_macros(rho_lin_ave1,mom_lin_ave1,energy_lin_ave1);
    end

    if(compute_error_2)
        filename_lin = plotname(2)+"_lin.mat";
        load(filename_lin);
        [u_lin2,theta_lin2] = compute_macros(rho_lin_ave2,mom_lin_ave2,energy_lin_ave2);
    end
    if(compute_error_3)
        filename_lin = plotname(3)+"_lin.mat";
        load(filename_lin);
        [u_lin3,theta_lin3] = compute_macros(rho_lin_ave3,mom_lin_ave3,energy_lin_ave3);
    end
    if(compute_error_4)
        filename_lin = plotname(4)+"_lin.mat";
        load(filename_lin);
        [u_lin4,theta_lin4] = compute_macros(rho_lin_ave4,mom_lin_ave4,energy_lin_ave4);
    end
    if(compute_error_5)
        filename_lin = plotname(5)+"_lin.mat";
        load(filename_lin);
        [u_lin5,theta_lin5] = compute_macros(rho_lin_ave5,mom_lin_ave5,energy_lin_ave5);
    end
end

%% Compute error
if compute_hme_error
    if compute_error_1
        hme_rho_diff1 = (rho_hme_ave1-sol_ref(:,1));
        hme_mom_diff1 = (mom_hme_ave1-sol_ref(:,2));
        hme_energy_diff1 = (energy_hme_ave1-sol_ref(:,3)); 
        hme_rho_error1 = norm(hme_rho_diff1,norm_number)/norm(sol_ref(:,1),norm_number);
        hme_mom_error1 = norm(hme_mom_diff1,norm_number)/norm(sol_ref(:,2),norm_number);
        hme_energy_error1 = norm(hme_energy_diff1,norm_number)/norm(sol_ref(:,3),norm_number);
    end
    if compute_error_2
        hme_rho_diff2 = (rho_hme_ave2-sol_ref(:,1));
        hme_mom_diff2 = (mom_hme_ave2-sol_ref(:,2));
        hme_energy_diff2 = (energy_hme_ave2-sol_ref(:,3)); 
        hme_rho_error2 = norm(hme_rho_diff2,norm_number)/norm(sol_ref(:,1),norm_number);
        hme_mom_error2 = norm(hme_mom_diff2,norm_number)/norm(sol_ref(:,2),norm_number);
        hme_energy_error2 = norm(hme_energy_diff2,norm_number)/norm(sol_ref(:,3),norm_number);
    end
    if compute_error_3
        hme_rho_diff3 = (rho_hme_ave3-sol_ref(:,1));
        hme_mom_diff3 = (mom_hme_ave3-sol_ref(:,2));
        hme_energy_diff3 = (energy_hme_ave3-sol_ref(:,3)); 
        hme_rho_error3 = norm(hme_rho_diff3,norm_number)/norm(sol_ref(:,1),norm_number);
        hme_mom_error3 = norm(hme_mom_diff3,norm_number)/norm(sol_ref(:,2),norm_number);
        hme_energy_error3 = norm(hme_energy_diff3,norm_number)/norm(sol_ref(:,3),norm_number);
    end
    if compute_error_4
        hme_rho_diff4 = (rho_hme_ave4-sol_ref(:,1));
        hme_mom_diff4 = (mom_hme_ave4-sol_ref(:,2));
        hme_energy_diff4 = (energy_hme_ave4-sol_ref(:,3)); 
        hme_rho_error4 = norm(hme_rho_diff4,norm_number)/norm(sol_ref(:,1),norm_number);
        hme_mom_error4 = norm(hme_mom_diff4,norm_number)/norm(sol_ref(:,2),norm_number);
        hme_energy_error4 = norm(hme_energy_diff4,norm_number)/norm(sol_ref(:,3),norm_number);
    end
    if compute_error_5
        hme_rho_diff5 = (rho_hme_ave5-sol_ref(:,1));
        hme_mom_diff5 = (mom_hme_ave5-sol_ref(:,2));
        hme_energy_diff5 = (energy_hme_ave5-sol_ref(:,3)); 
        hme_rho_error5 = norm(hme_rho_diff5,norm_number)/norm(sol_ref(:,1),norm_number);
        hme_mom_error5 = norm(hme_mom_diff5,norm_number)/norm(sol_ref(:,2),norm_number);
        hme_energy_error5 = norm(hme_energy_diff5,norm_number)/norm(sol_ref(:,3),norm_number);
    end
end
if compute_qbme_error
    if compute_error_1
        qbme_rho_diff1 = (rho_qbme_ave1-sol_ref(:,1));
        qbme_mom_diff1 = (mom_qbme_ave1-sol_ref(:,2));
        qbme_energy_diff1 = (energy_qbme_ave1-sol_ref(:,3)); 
        qbme_rho_error1 = norm(qbme_rho_diff1,norm_number)/norm(sol_ref(:,1),norm_number);
        qbme_mom_error1 = norm(qbme_mom_diff1,norm_number)/norm(sol_ref(:,2),norm_number);
        qbme_energy_error1 = norm(qbme_energy_diff1,norm_number)/norm(sol_ref(:,3),norm_number);
    end
    if compute_error_2
        qbme_rho_diff2 = (rho_qbme_ave2-sol_ref(:,1));
        qbme_mom_diff2 = (mom_qbme_ave2-sol_ref(:,2));
        qbme_energy_diff2 = (energy_qbme_ave2-sol_ref(:,3)); 
        qbme_rho_error2 = norm(qbme_rho_diff2,norm_number)/norm(sol_ref(:,1),norm_number);
        qbme_mom_error2 = norm(qbme_mom_diff2,norm_number)/norm(sol_ref(:,2),norm_number);
        qbme_energy_error2 = norm(qbme_energy_diff2,norm_number)/norm(sol_ref(:,3),norm_number);
    end
    if compute_error_3
        qbme_rho_diff3 = (rho_qbme_ave3-sol_ref(:,1));
        qbme_mom_diff3 = (mom_qbme_ave3-sol_ref(:,2));
        qbme_energy_diff3 = (energy_qbme_ave3-sol_ref(:,3)); 
        qbme_rho_error3 = norm(qbme_rho_diff3,norm_number)/norm(sol_ref(:,1),norm_number);
        qbme_mom_error3 = norm(qbme_mom_diff3,norm_number)/norm(sol_ref(:,2),norm_number);
        qbme_energy_error3 = norm(qbme_energy_diff3,norm_number)/norm(sol_ref(:,3),norm_number);
    end
    if compute_error_4
        qbme_rho_diff4 = (rho_qbme_ave4-sol_ref(:,1));
        qbme_mom_diff4 = (mom_qbme_ave4-sol_ref(:,2));
        qbme_energy_diff4 = (energy_qbme_ave4-sol_ref(:,3)); 
        qbme_rho_error4 = norm(qbme_rho_diff4,norm_number)/norm(sol_ref(:,1),norm_number);
        qbme_mom_error4 = norm(qbme_mom_diff4,norm_number)/norm(sol_ref(:,2),norm_number);
        qbme_energy_error4 = norm(qbme_energy_diff4,norm_number)/norm(sol_ref(:,3),norm_number);
    end
    if compute_error_5
        qbme_rho_diff5 = (rho_qbme_ave5-sol_ref(:,1));
        qbme_mom_diff5 = (mom_qbme_ave5-sol_ref(:,2));
        qbme_energy_diff5 = (energy_qbme_ave5-sol_ref(:,3)); 
        qbme_rho_error5 = norm(qbme_rho_diff5,norm_number)/norm(sol_ref(:,1),norm_number);
        qbme_mom_error5 = norm(qbme_mom_diff5,norm_number)/norm(sol_ref(:,2),norm_number);
        qbme_energy_error5 = norm(qbme_energy_diff5,norm_number)/norm(sol_ref(:,3),norm_number);
    end
end
if compute_lin_error
    if compute_error_1
        lin_rho_diff1 = (rho_lin_ave1-sol_ref(:,1));
        lin_mom_diff1 = (mom_lin_ave1-sol_ref(:,2));
        lin_energy_diff1 = (energy_lin_ave1-sol_ref(:,3)); 
        lin_rho_error1 = norm(lin_rho_diff1,norm_number)/norm(sol_ref(:,1),norm_number);
        lin_mom_error1 = norm(lin_mom_diff1,norm_number)/norm(sol_ref(:,2),norm_number);
        lin_energy_error1 = norm(lin_energy_diff1,norm_number)/norm(sol_ref(:,3),norm_number);
    end
    if compute_error_2
        lin_rho_diff2 = (rho_lin_ave2-sol_ref(:,1));
        lin_mom_diff2 = (mom_lin_ave2-sol_ref(:,2));
        lin_energy_diff2 = (energy_lin_ave2-sol_ref(:,3)); 
        lin_rho_error2 = norm(lin_rho_diff2,norm_number)/norm(sol_ref(:,1),norm_number);
        lin_mom_error2 = norm(lin_mom_diff2,norm_number)/norm(sol_ref(:,2),norm_number);
        lin_energy_error2 = norm(lin_energy_diff2,norm_number)/norm(sol_ref(:,3),norm_number);
    end
    if compute_error_3
        lin_rho_diff3 = (rho_lin_ave3-sol_ref(:,1));
        lin_mom_diff3 = (mom_lin_ave3-sol_ref(:,2));
        lin_energy_diff3 = (energy_lin_ave3-sol_ref(:,3)); 
        lin_rho_error3 = norm(lin_rho_diff3,norm_number)/norm(sol_ref(:,1),norm_number);
        lin_mom_error3 = norm(lin_mom_diff3,norm_number)/norm(sol_ref(:,2),norm_number);
        lin_energy_error3 = norm(lin_energy_diff3,norm_number)/norm(sol_ref(:,3),norm_number);
    end
    if compute_error_4
        lin_rho_diff4 = (rho_lin_ave4-sol_ref(:,1));
        lin_mom_diff4 = (mom_lin_ave4-sol_ref(:,2));
        lin_energy_diff4 = (energy_lin_ave4-sol_ref(:,3)); 
        lin_rho_error4 = norm(lin_rho_diff4,norm_number)/norm(sol_ref(:,1),norm_number);
        lin_mom_error4 = norm(lin_mom_diff4,norm_number)/norm(sol_ref(:,2),norm_number);
        lin_energy_error4 = norm(lin_energy_diff4,norm_number)/norm(sol_ref(:,3),norm_number);
    end
    if compute_error_5
        lin_rho_diff5 = (rho_lin_ave5-sol_ref(:,1));
        lin_mom_diff5 = (mom_lin_ave5-sol_ref(:,2));
        lin_energy_diff5 = (energy_lin_ave5-sol_ref(:,3)); 
        lin_rho_error5 = norm(lin_rho_diff5,norm_number)/norm(sol_ref(:,1),norm_number);
        lin_mom_error5 = norm(lin_mom_diff5,norm_number)/norm(sol_ref(:,2),norm_number);
        lin_energy_error5 = norm(lin_energy_diff5,norm_number)/norm(sol_ref(:,3),norm_number);
    end
end

if save_errors
    if compute_error_1
       save(plotname(1)+"_errors", ...
            "hme_rho_error1","lin_rho_error1","qbme_rho_error1", ...
            "hme_mom_error1","lin_mom_error1","qbme_mom_error1", ...
            "hme_energy_error1","lin_energy_error1","qbme_energy_error1"); 
    end
    if compute_error_2
        save(plotname(2)+"_errors",...
            "hme_rho_error2","lin_rho_error2","qbme_rho_error2",...
            "hme_mom_error2","lin_mom_error2","qbme_mom_error2",...
            "hme_energy_error2","lin_energy_error2","qbme_energy_error2");
    end
    if compute_error_3
        save(plotname(3)+"_errors",...
            "hme_rho_error3","lin_rho_error3","qbme_rho_error3",...
            "hme_mom_error3","lin_mom_error3","qbme_mom_error3",...
            "hme_energy_error3","lin_energy_error3","qbme_energy_error3");
    end
    if compute_error_4
        save(plotname(4)+"_errors",...
            "hme_rho_error4","lin_rho_error4","qbme_rho_error4",...
            "hme_mom_error4","lin_mom_error4","qbme_mom_error4",...
            "hme_energy_error4","lin_energy_error4","qbme_energy_error4");
    end
    if compute_error_5
        save(plotname(5)+"_errors",...
            "hme_rho_error5","lin_rho_error5","qbme_rho_error5",...
            "hme_mom_error5","lin_mom_error5","qbme_mom_error5",...
            "hme_energy_error5","lin_energy_error5","qbme_energy_error5");
    end
end

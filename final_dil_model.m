
%% Figure Hiding
set(0,'DefaultFigureVisible','off');
%% figure off and variable clearance
clc
clear variables
h = [];
%% Initial Conditions
%sol5c = bvp5c(@bvpfcn, @bcfcn, solinit, opts);
td__c = 0.01;
theta_cc_mul=[0.01,0.1,1,2];
theta_cc_mul_nd = theta_cc_mul/td__c;
%% for positive & Negative AEM
col_dat = ["#0072BD","#D95319", 	"#7E2F8E","#77AC30"];
h(1) = figure;
omega__c = 1;
%Claculating BCs
C__posL2 = 1;
C__posR2 = 10;
diff_cons = 0.5;
i__Val = 0;
Pe__Val = 0;
bcs_0 = zeros(length(theta_cc_mul),1);
bcs_1 = zeros(length(theta_cc_mul),1);
sol4c_yd = [];
for i = 1:length(theta_cc_mul)
    bcs_0(i) = -omega__c*theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posL2^2)/2;
    bcs_1(i) = -omega__c*theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posR2^2)/2;
end
for i=1:length(theta_cc_mul)
    xi = linspace(0,1,200);
    xmesh = linspace(0, 1, 1000);
    opts = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",diff_cons,i__Val,Pe__Val),'RelTol',0.001,'AbsTol',0.001,'Stats','on');
    solinit = bvpinit(xmesh,[0.8; 1]);
    sol4c = bvp4c(bvp_call_mod_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",diff_cons,i__Val,Pe__Val), bcf_call(bcs_0(i),bcs_1(i)),solinit, opts);
    sol4c_yd(:,i) = sol4c.y(1,:);
    plot(sol4c.x,sol4c.y(1,:),'LineWidth',3,'Color',col_dat(i));
    hold on;
end
title("C_{pos}^{*}(\xi) vs \xi for AEM")
    legend("\theta^* = 1","\theta^* = 10","\theta^* = 100","\theta^* = 200",'location', 'best','Orientation','horizontal')
    xlabel('\xi')
    ylabel('C_{pos}^{*}(\xi)')
    hold off;
    saveas(h(1),"plots_sv/f3_cpos_aem.png")
    fprintf(2,"Plot: "+ num2str(1)+" Done\n");
% for negative AEM
h(2) = figure;
for i=1:length(theta_cc_mul)
    plot(sol4c.x,sol4c_yd(:,i)+ omega__c*theta_cc_mul_nd(i),'LineWidth',3,'Color',col_dat(i));
    hold on;
end
title("C_{neg}^{*}(\xi) vs \xi for AEM")
    legend("\theta^* = 1","\theta^* = 10","\theta^* = 100","\theta^* = 200",'location', 'best','Orientation','horizontal')
    xlabel('\xi')
    ylabel('C_{neg}^{*}(\xi)')
    hold off;
 saveas(h(2),"plots_sv/f3_cneg_aem.png")
 fprintf(2,"Plot: "+ num2str(2)+" Done\n");
%% for positive & negative CEM
%for positive CEM
h(3) = figure;
omega__c = -1;
diff_cons = 0.5;
%Claculating BCs
C__posL2 = 1;
C__posR2 = 10;
bcs_0 = [];
bcs_1 = [];
sol4c_yd = [];
for i = 1:length(theta_cc_mul)
    bcs_0(i) = -omega__c*theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posL2^2)/2;
    bcs_1(i) = -omega__c*theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posR2^2)/2;
end
for i=1:length(theta_cc_mul)
    xmesh = linspace(0, 1, 1000);
    opts = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",diff_cons,0,0),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
    solinit = bvpinit(xmesh,[0; -1]);
    sol4c = bvp4c(bvp_call_mod_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",diff_cons,0,0), bcf_call(bcs_0(i),bcs_1(i)),solinit, opts);
    sol4c_yd(:,i) = sol4c.y(1,:); 
    plot(sol4c.x,sol4c.y(1,:),'LineWidth',3,'Color',col_dat(i));
    hold on;
end
    title("C_{pos}^{*}(\xi) vs \xi for CEM")
    legend("\theta^* = 1","\theta^* = 10","\theta^* = 100","\theta^* = 200",'location', 'best','Orientation','horizontal')
    xlabel('\xi')
    ylabel('C_{pos}^{*}(\xi)')
    hold off;
    saveas(h(3),"plots_sv/f3_cpos_cem.png")
    fprintf(2,"Plot: "+ num2str(3)+" Done\n");  
% For negative CEM
h(4)=figure;
for i=1:length(theta_cc_mul)
    plot(sol4c.x,sol4c_yd(:,i)+omega__c*theta_cc_mul_nd(i),'LineWidth',3,'Color',col_dat(i));
    hold on;
end
    title("C_{neg}^{*}(\xi) vs \xi for CEM")
    legend("\theta^* = 1","\theta^* = 10","\theta^* = 100","\theta^* = 200",'location', 'best','Orientation','horizontal')
    xlabel('\xi')
    ylabel('C_{neg}^{*}(\xi)')
    hold off;
    saveas(h(4),"plots_sv/f3_cneg_cem.png")
    fprintf(2,"Plot: "+ num2str(4)+" Done\n");
%% For  Calculating positive Flux in AEM
% For positive Flux in AEM
shp = '-o';
diff_cons = 0.5;
Pe__val = 0;
i__val = 0;
theta_cc_mul=linspace(0.001,4,80);
theta_cc_mul_nd = theta_cc_mul/td__c;
h(5) = figure;
omega__c = 1;
%Claculating BCs
C__posL2 = 1;
C__posR2 = 10;
N_poss = [];
bcs_0 = [];
bcs_1 = [];
for i = 1:length(theta_cc_mul)
    bcs_0(i) = -omega__c*theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posL2^2)/2;
    bcs_1(i) = -omega__c*theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posR2^2)/2;
end
for i=1:length(theta_cc_mul)
    xmesh = linspace(0, 1, 100);
    opts = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",diff_cons,i__val,Pe__val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
    solinit = bvpinit(xmesh,[100; -0.01]);
    sol4c = bvp4c(bvp_call_mod_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",diff_cons,i__val,Pe__val), bcf_call(bcs_0(i),bcs_1(i)),solinit, opts);
    [sol4c.y(1,1),sol4c.y(2,1)]
    N_poss_0 = zeros(length(sol4c.x),1);
    d__pos = diff_cons;
    A = 0;
    B = theta_cc_mul_nd(i);
    for k = 1:length(sol4c.x)
        C__1 = sol4c.y(1,k);
        C__2 = sol4c.y(2,k);
        xi = sol4c.x(1,k);
      N_poss_0(k) = -(d__pos * C__2) + d__pos * (-omega__c * A + ((-1 + d__pos) * C__2) + omega__c * (A * (xi - 0.5e0) + B) * Pe__val + i__val) * C__1 / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) + C__1 * Pe__val;
    end
    N_poss(i) = mean(N_poss_0);
end
plot(theta_cc_mul_nd,N_poss,shp,'LineWidth',3)
title("N_{pos}^{*} vs \theta^* for AEM")
    xlabel('\theta^*')
    ylabel('N_{pos}^{*}(\theta^*)')
    hold off;
    saveas(h(5),"plots_sv/f3_Npos_aem.png")
    fprintf(2,"Plot: "+ num2str(5)+" Done\n");
%% For Calculating negative Flux in CEM
h(6) = figure;
omega__c = -1;
diff_cons = 0.5;
%Claculating BCs
C__posL2 = 1;
C__posR2 = 10;
Pe_val = 0;
bcs_0 = [];
bcs_1 = [];
sol4c_yd = [];
for i = 1:length(theta_cc_mul)
    bcs_0(i) = -omega__c*theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posL2^2)/2;
    bcs_1(i) = -omega__c*theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posR2^2)/2;
end
for i=1:length(theta_cc_mul)
    xmesh = linspace(0, 1, 100);
    opts = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",diff_cons,0,Pe_val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
    solinit = bvpinit(xmesh,[300; -0]);
    sol4c = bvp4c(bvp_call_mod_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",diff_cons,0,Pe_val), bcf_call(bcs_0(i),bcs_1(i)),solinit, opts);
    N_negs_0 = zeros(length(sol4c.x),1);
    A = 0;
    B = theta_cc_mul_nd(i);
    d_pos = diff_cons;
    for k = 1:length(sol4c.x)
        C__1 = sol4c.y(1,k);
        C__2 = sol4c.y(2,k);
        xi = sol4c.x(1,k);
        N_negs_0(k) = (-d__pos * ((A * (xi - 0.5e0) + B) * omega__c + (2 * C__1)) * C__2 - A * C__1 * omega__c * d__pos + (Pe__val * (0.1e1 + d__pos) * C__1 - i__val) * (C__1 + (A * (xi - 0.5e0) + B) * omega__c)) / ((A * (xi - 0.5e0) + B) * omega__c + C__1 * (0.1e1 + d__pos));
    end
    N_negs(i) = mean(N_negs_0);
end
plot(theta_cc_mul_nd,N_negs,'-o','LineWidth',3)
title("N_{neg}^{*} vs \theta^* for AEM")
    xlabel('\theta^*')
    ylabel('N_{neg}^{*}(\theta^*)')
    hold off;
    saveas(h(6),"plots_sv/f3_Nneg_cem.png")
    fprintf(2,"Plot: "+ num2str(6)+" Done\n");
%% Combined postive and negative flux
h(7) = figure;
plot(theta_cc_mul_nd,N_poss,shp,'LineWidth',3);hold on;
plot(theta_cc_mul_nd,N_negs,shp,'LineWidth',3)
title("N^{*} vs \theta^*")
    xlabel('\theta')
    ylabel('N^{*}(\theta^*)')
    legend("N_{pos} for AEM","N_{neg} for CEM",'location', 'best')
    hold off;
    saveas(h(7),"plots_sv/f3_Nneg_Npos_combined.png")
    fprintf(2,"Plot: "+ num2str(7)+" Done\n");
%% Diffusion coeff variation for pos for AEM @ theta = 2
h(8) = figure;
Pe__val = 0;
omega__c = 1;
theta_const = 2/td__c;
diff_cnst = [0.025 1 10 40];
i_const = 0;
%Claculating BCs
C__posL2 = 1;
C__posR2 = 10;
    bcs_0 = -omega__c*theta_const/2 + sqrt(theta_const^2 + 4*C__posL2^2)/2;
    bcs_1 = -omega__c*theta_const/2 + sqrt(theta_const^2 + 4*C__posR2^2)/2;
for i=1:length(diff_cnst)
    xmesh = linspace(0, 1, 100);
    opts = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_const,omega__c,"pos",diff_cnst(i),i_const,Pe__val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
    solinit = bvpinit(xmesh,[0;0.15]);
    sol4c = bvp4c(bvp_call_mod_theta_pe_var(0,theta_const,omega__c,"pos",diff_cnst(i),i_const,Pe__val), bcf_call(bcs_0,bcs_1),solinit, opts);
    plot(sol4c.x,sol4c.y(1,:),'LineWidth',3,'Color',col_dat(i));
    hold on;
end
    title("C_{pos}^{*}(\xi) vs \xi in AEM for different ratios of the diffusion coefficient")
    legend("D_{+}/D_{-}=0.025","D_{+}/D_{-}=1","D_{+}/D_{-}=10","D_{+}/D_{-}=40",'location', 'best')
    xlabel('\xi')
    ylabel('C_{pos}^{*}(\xi)')
    hold off;
    saveas(h(8),"plots_sv/f4_diff_cpos_aem.png")
    fprintf(2,"Plot: "+ num2str(8)+" Done\n");
%% Diffusion coeff variation for pos for CEM @ theta = 2
h(9) = figure;
Pe__val = 0;
omega__c = -1;
theta_const = 2/td__c;
i_const  = 0;
diff_cnst = [0.025 1 10 40];
%Claculating BCs
C__posL2 = 1;
C__posR2 = 10;
    bcs_0 = -omega__c*theta_const/2 + sqrt(theta_const^2 + 4*C__posL2^2)/2;
    bcs_1 = -omega__c*theta_const/2 + sqrt(theta_const^2 + 4*C__posR2^2)/2;
for i=1:length(diff_cnst)
    xmesh = linspace(0, 1, 100);
    opts = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_const,omega__c,"pos",diff_cnst(i),i_const,Pe__val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
    solinit = bvpinit(xmesh,[200.005; 0.4938]);
    sol4c = bvp4c(bvp_call_mod_theta_pe_var(0,theta_const,omega__c,"pos",diff_cnst(i),i_const,Pe__val), bcf_call(bcs_0,bcs_1),solinit, opts);
    plot(sol4c.x,sol4c.y(1,:),'LineWidth',3,'Color',col_dat(i));
    hold on;
end
    title("C_{pos}^{*}(\xi) vs \xi in CEM for different ratios of the diffusion coefficient")
    legend("D_{+}/D_{-}=0.025","D_{+}/D_{-}=1","D_{+}/D_{-}=10","D_{+}/D_{-}=40",'location', 'best')
    xlabel('\xi')
    ylabel('C_{pos}^{*}(\xi)')
    hold off;
    saveas(h(9),"plots_sv/f4_diff_cpos_cem.png")
    fprintf(2,"Plot: "+ num2str(9)+" Done\n");
%% For Calculating positive Flux in CEM & AEM Charge
col_dat_22 = ["#0072BD","#D95319", 	"#7E2F8E","#77AC30"];
%for i =500,-500,1000,-1000 for CEM
i_act = [500,1000]; % Actual Current Density
i_nond = 1;
i_dens = i_act/i_nond; % Dimensionless Current Density
% For  Calculating positive Flux in CEM
theta_cc_mul=linspace(0.01,4,80);
theta_cc_mul_nd = theta_cc_mul/td__c;
h(10) = figure;
diff_cons = 0.5;
Pe__val = 0;
omega__c_rng = [1,-1];
%Claculating BCs
C__posL2 = 1;
C__posR2 = 10;
for ii = 1:2
    if ii==1
        ii_c = [1 2];
    else
        ii_c = [3 4];
    end
    omega__c = omega__c_rng(ii);
for j = 1:length(i_act)
    N_poss_cem_1 = zeros(length(theta_cc_mul),1);
    bcs_0 = zeros(length(theta_cc_mul),1);
    bcs_1 = zeros(length(theta_cc_mul),1);
    for i = 1:length(theta_cc_mul)
        bcs_0(i) = -omega__c*theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posL2^2)/2;
        bcs_1(i) = -omega__c*theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posR2^2)/2;
    end
    for i=1:length(theta_cc_mul)
        xmesh = linspace(0, 1, 100);
        opts_0 = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",diff_cons,0,Pe__val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
        solinit_0 = bvpinit(xmesh,[10; 1]);
        sol4c_0 = bvp4c(bvp_call_mod_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",diff_cons,0,Pe__val), bcf_call(bcs_0(i),bcs_1(i)),solinit_0, opts_0);
        if omega__c == 1
        solinit = bvpinit(xmesh,[10;1]);
        else
        solinit = bvpinit(xmesh,[sol4c_0.y(1,1); sol4c_0.y(2,1)]);
        end
        opts = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",diff_cons,i_dens(j),Pe__val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
        sol4c = bvp4c(bvp_call_mod_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",diff_cons,i_dens(j),Pe__val), bcf_call(bcs_0(i),bcs_1(i)),solinit, opts);
        N_poss_cem_1_0 = zeros(length(sol4c.x),1);
        i__val = i_dens(j);
        A = 0;
        B = theta_cc_mul_nd(i);
        d__pos = diff_cons;
        for k = 1:length(sol4c.x)
            C__1 = sol4c.y(1,k);
            C__2 = sol4c.y(2,k);
            xi = sol4c.x(1,k);
            N_poss_cem_1_0 = -(d__pos * C__2) + d__pos * (-omega__c * A + ((-1 + d__pos) * C__2) + omega__c * (A * (xi - 0.5e0) + B) * Pe__val + i__val) * C__1 / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) + C__1 * Pe__val;
        end
        N_poss_cem_1(i) = mean(N_poss_cem_1_0);
    end
plot(theta_cc_mul_nd,N_poss_cem_1,'-o','LineWidth',3,'Color',col_dat_22(ii_c(j)));hold on;
end
end
title("Cation Fluxes for various current densities in ")
legend("N_{pos}^* for CEM(500)","N_{pos}^*  for CEM(1000)","N_{pos}^*  for AEM(500)","N_{pos}^*  for AEM(1000)",'location', 'best')
    xlabel('\theta')
    ylabel('N_{pos}^{*}(\theta^*)')
saveas(h(10),"plots_sv/f7_pos_i_flux_cem_charge.png")
fprintf(2,"Plot: "+ num2str(10)+" Done\n");
hold off;
%% For Calculating positive Flux in CEM & AEM disCharge
col_dat_22 = ["#0072BD","#D95319", 	"#7E2F8E","#77AC30"];
%for i =500,-500,1000,-1000 for CEM
i_act = [-500,-1000]; % Actual Current Density
i_nond = 1;
i_dens = i_act/i_nond; % Dimensionless Current Density
% For  Calculating positive Flux in CEM
theta_cc_mul=linspace(0.01,4,80);
theta_cc_mul_nd = theta_cc_mul/td__c;
h(11) = figure;
diff_cons = 0.5;
Pe__val = 0;
omega__c_rng = [1,-1];
%Claculating BCs
C__posL2 = 1;
C__posR2 = 10;
for ii = 1:2
    if ii==1
        ii_c = [1 2];
    else
        ii_c = [3 4];
    end
    omega__c = omega__c_rng(ii);
for j = 1:length(i_act)
    N_poss_cem_1 = zeros(length(theta_cc_mul),1);
    bcs_0 = zeros(length(theta_cc_mul),1);
    bcs_1 = zeros(length(theta_cc_mul),1);
    for i = 1:length(theta_cc_mul)
        bcs_0(i) = -omega__c*theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posL2^2)/2;
        bcs_1(i) = -omega__c*theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posR2^2)/2;
    end
    for i=1:length(theta_cc_mul)
        xmesh = linspace(0, 1, 100);
        opts_0 = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",diff_cons,0,Pe__val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
        solinit_0 = bvpinit(xmesh,[10; 1]);
        sol4c_0 = bvp5c(bvp_call_mod_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",diff_cons,0,Pe__val), bcf_call(bcs_0(i),bcs_1(i)),solinit_0, opts_0);
        if omega__c == -1
        solinit = bvpinit(xmesh,[100;-0.01]);
        else
        solinit = bvpinit(xmesh,[sol4c_0.y(1,1); sol4c_0.y(2,1)]);
        end
        opts = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",diff_cons,i_dens(j),Pe__val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
        sol4c = bvp4c(bvp_call_mod_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",diff_cons,i_dens(j),Pe__val), bcf_call(bcs_0(i),bcs_1(i)),solinit, opts);
        N_poss_cem_1_0 = zeros(length(sol4c.x),1);
        i__val = i_dens(j);
        A = 0;
        B = theta_cc_mul_nd(i);
        d__pos = diff_cons;
        for k = 1:length(sol4c.x)
            C__1 = sol4c.y(1,k);
            C__2 = sol4c.y(2,k);
            xi = sol4c.x(1,k);
            N_poss_cem_1_0 = -(d__pos * C__2) + d__pos * (-omega__c * A + ((-1 + d__pos) * C__2) + omega__c * (A * (xi - 0.5e0) + B) * Pe__val + i__val) * C__1 / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) + C__1 * Pe__val;
        end
        N_poss_cem_1(i) = mean(N_poss_cem_1_0);
    end
plot(theta_cc_mul_nd,N_poss_cem_1,'-o','LineWidth',3,'Color',col_dat_22(ii_c(j)));hold on;
end
end
title("Cation Fluxes for various current densities in CEM")
legend("N_{pos}^* for CEM(500)","N_{pos}^*  for CEM(1000)","N_{pos}^*  for AEM(500)","N_{pos}^*  for AEM(1000)",'location', 'best')
    xlabel('\theta')
    ylabel('N_{pos}^{*}(\theta^*)')
saveas(h(11),"plots_sv/f7_pos_i_flux_cem_discharge.png")
fprintf(2,"Plot: "+ num2str(10)+" Done\n");
hold off;
%% For co-ion and counter ion Fluxes in CEM as function of current density
i_act = linspace(-1000,1000,10); % Actual Current Density
i_dens = i_act; % Dimensionless Current Density
xmesh = linspace(0, 1, 100);
Pe__val = 0;
diff_cons = 0.5;
% For  Calculating positive Flux in CEM
theta_cc_mul=2;
theta_cc_mul_nd = theta_cc_mul/td__c;
h(12) = figure;
omega__c = -1;
%Claculating BCs
C__posL2 = 1;
C__posR2 = 10;
N_poss_cem_1 = zeros(length(i_act),1);
N_negs_cem_1 = zeros(length(i_act),1);
for j = 1:length(i_act)
        bcs_0 = theta_cc_mul_nd/2 + sqrt(theta_cc_mul_nd^2 + 4*C__posL2^2)/2;
        bcs_1 = theta_cc_mul_nd/2 + sqrt(theta_cc_mul_nd^2 + 4*C__posR2^2)/2;
    opts_0 = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_cc_mul_nd,omega__c,"pos",diff_cons,0,Pe__Val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
    solinit_0 = bvpinit(xmesh,[0; 1]);
    sol4c_0 = bvp4c(bvp_call_mod_theta_pe_var(0,theta_cc_mul_nd,omega__c,"pos",diff_cons,0,Pe__Val), bcf_call(bcs_0,bcs_1),solinit_0, opts_0);
    if i_dens(j)>=0
    solinit = bvpinit(xmesh,[sol4c_0.y(1,1); sol4c_0.y(2,1)]);
    else
    solinit = bvpinit(xmesh,[10; 1]);
    end
    opts = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_cc_mul_nd,omega__c,"pos",diff_cons,i_dens(j),Pe__Val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
    sol4c = bvp4c(bvp_call_mod_theta_pe_var(0,theta_cc_mul_nd,omega__c,"pos",diff_cons,i_dens(j),Pe__Val), bcf_call(bcs_0,bcs_1),solinit, opts);
    N_poss_cem_1_0 = zeros(length(sol4c.x),1);
    N_negs_cem_1_0 = zeros(length(sol4c.x),1);
    i__val = i_dens(j);
    d__pos = diff_cons;
    A = 0;
    B = theta_cc_mul_nd;
    for k = 1:length(sol4c.x)
        C__1 = sol4c.y(1,k);
        C__2 = sol4c.y(2,k);
        xi = sol4c.x(1,k);
        N_poss_cem_1_0(k) = -(d__pos * C__2) + d__pos * (-omega__c * A + ((-1 + d__pos) * C__2) + omega__c * (A * (xi - 0.5e0) + B) * Pe__val + i__val) * C__1 / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) + C__1 * Pe__val;
        N_negs_cem_1_0(k) = (-d__pos * ((A * (xi - 0.5e0) + B) * omega__c + (2 * C__1)) * C__2 - A * C__1 * omega__c * d__pos + (Pe__val * (0.1e1 + d__pos) * C__1 - i__val) * (C__1 + (A * (xi - 0.5e0) + B) * omega__c)) / ((A * (xi - 0.5e0) + B) * omega__c + C__1 * (0.1e1 + d__pos));
    end
N_poss_cem_1(j) = mean(N_poss_cem_1_0);
N_negs_cem_1(j) = mean(N_negs_cem_1_0);
end
plot(i_act,N_poss_cem_1,'-o','LineWidth',3);hold on;
plot(i_act,N_negs_cem_1,'-o','LineWidth',3);
 title("Co-ion and Counter ion flux as func of i_e*")
legend("N_{pos}^* for CEM","N_{neg}^* for CEM",'location', 'best')
    xlabel('i_e^*(current density)')
    ylabel('N_{i}^{*}')
hold off;
saveas(h(12),"plots_sv/f7_pos_i_var_flux_cem.png")
fprintf(2,"Plot: "+ num2str(14)+" Done\n");
%% For co-ion and counter ion Fluxes in AEM as function of current density
i_act = linspace(-1000,1000,10); % Actual Current Density
i_dens = i_act; % Dimensionless Current Density
xmesh = linspace(0, 1, 1000);
Pe__val = 0;
diff_cons = 0.5;
% For  Calculating positive Flux in CEM
theta_cc_mul=2;
theta_cc_mul_nd = theta_cc_mul/td__c;
h(13) = figure;
omega__c = 1;
%Claculating BCs
C__posL2 = 1;
C__posR2 = 10;
N_poss_cem_1 = zeros(length(i_act),1);
N_negs_cem_1 = zeros(length(i_act),1);
for j = 1:length(i_act)
        bcs_0 = -theta_cc_mul_nd/2 + sqrt(theta_cc_mul_nd^2 + 4*C__posL2^2)/2;
        bcs_1 = -theta_cc_mul_nd/2 + sqrt(theta_cc_mul_nd^2 + 4*C__posR2^2)/2;
    opts_0 = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_cc_mul_nd,omega__c,"pos",diff_cons,0,Pe__Val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
    solinit_0 = bvpinit(xmesh,[0; 1]);
    sol4c_0 = bvp4c(bvp_call_mod_theta_pe_var(0,theta_cc_mul_nd,omega__c,"pos",diff_cons,0,Pe__Val), bcf_call(bcs_0,bcs_1),solinit_0, opts_0);
    if i_dens(j)>=0
    solinit = bvpinit(xmesh,[sol4c_0.y(1,1); sol4c_0.y(2,1)]);
    else
    solinit = bvpinit(xmesh,[10; -0.01]);
    end
    opts = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_cc_mul_nd,omega__c,"pos",diff_cons,i_dens(j),Pe__Val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
    sol4c = bvp4c(bvp_call_mod_theta_pe_var(0,theta_cc_mul_nd,omega__c,"pos",diff_cons,i_dens(j),Pe__Val), bcf_call(bcs_0,bcs_1),solinit, opts);
    N_poss_cem_1_0 = zeros(length(sol4c.x),1);
    N_negs_cem_1_0 = zeros(length(sol4c.x),1);
    i__val = i_dens(j);
    d__pos = diff_cons;
    A = 0;
    B = theta_cc_mul_nd;
    for k = 1:length(sol4c.x)
        C__1 = sol4c.y(1,k);
        C__2 = sol4c.y(2,k);
        xi = sol4c.x(1,k);
        N_poss_cem_1_0(k) = -(d__pos * C__2) + d__pos * (-omega__c * A + ((-1 + d__pos) * C__2) + omega__c * (A * (xi - 0.5e0) + B) * Pe__val + i__val) * C__1 / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) + C__1 * Pe__val;
        N_negs_cem_1_0(k) = (-d__pos * ((A * (xi - 0.5e0) + B) * omega__c + (2 * C__1)) * C__2 - A * C__1 * omega__c * d__pos + (Pe__val * (0.1e1 + d__pos) * C__1 - i__val) * (C__1 + (A * (xi - 0.5e0) + B) * omega__c)) / ((A * (xi - 0.5e0) + B) * omega__c + C__1 * (0.1e1 + d__pos));
    end
N_poss_cem_1(j) = mean(N_poss_cem_1_0);
N_negs_cem_1(j) = mean(N_negs_cem_1_0);
end
plot(i_act,N_poss_cem_1,'-o','LineWidth',3,'Color',"blue");hold on;
plot(i_act,N_negs_cem_1,'-o','LineWidth',3,'Color',"cyan");
 title("Co-ion and Counter ion flux as func of i_e*")
legend("N_{pos}^* for CEM","N_{neg}^* for AEM",'location', 'best')
    xlabel('i_e^*(current density)')
    ylabel('N_{i}^{*}')
hold off;
saveas(h(13),"plots_sv/f7_pos_i_var_flux_aem.png")
fprintf(2,"Plot: "+ num2str(14)+" Done\n");
%% Diffusion coeff variation for pos for AEM @ theta = 2 for i charge
col_dat_22 = ["#0072BD","#D95319", 	"#7E2F8E","#77AC30"];
omega_rng = [-1,1];
i_rng = [-500,500];
for w = 1:length(omega_rng)
for ii = 1:length(i_rng)
h(14,w-1+ii) = figure;
Pe__val = 0;
omega__c = omega_rng(w);
theta_const = 2/td__c;
diff_cnst = [0.025 1 10 40];
i_orig = i_rng(ii);
i_nond = 1;
i_const  = i_orig/i_nond;
%Claculating BCs
C__posL2 = 1;
C__posR2 = 10;
    bcs_0 = -omega__c*theta_const/2 + sqrt(theta_const^2 + 4*C__posL2^2)/2;
    bcs_1 = -omega__c*theta_const/2 + sqrt(theta_const^2 + 4*C__posR2^2)/2;
for i=1:length(diff_cnst)
    xmesh = linspace(0, 1, 1000);
    opts = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_const,omega__c,"pos",diff_cnst(i),i_const,Pe__val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
    solinit = bvpinit(xmesh,[0.995;0.5348]);
    sol4c = bvp4c(bvp_call_mod_theta_pe_var(0,theta_const,omega__c,"pos",diff_cnst(i),i_const,Pe__val), bcf_call(bcs_0,bcs_1),solinit, opts);
    plot(sol4c.x,sol4c.y(1,:),'LineWidth',3,'Color',col_dat_22(i));
    hold on;
end
if omega__c == 1
title("C_{pos}^{*}(\xi) vs \xi in AEM @ i= "+num2str(i_orig)+" for various D_{+}/D_{-}")
else
title("C_{pos}^{*}(\xi) vs \xi in CEM @ i= "+num2str(i_orig)+" for various D_{+}/D_{-}")
end
    legend("D_{+}/D_{-}=0.025","D_{+}/D_{-}=1","D_{+}/D_{-}=10","D_{+}/D_{-}=40",'location', 'best')
    xlabel('\xi')
    ylabel('C_{pos}^{*}(\xi)')
    hold off;
    saveas(h(14,w-1+ii),"plots_sv/f8_diff_cpos_i_ch_w="+num2str(omega__c)+num2str(i_const)+".png")
    fprintf(2,"Plot: "+ num2str(12.1)+" Done\n");
end
end
%% For  Calculating positive Flux in AEM for inhomogeneous one side 
        xmesh = linspace(0, 1,100);
        i__val = 0;
        xi_m = xmesh;
        diff_cnst = 0.5;
        Pe__val = 0;
        h(15) = figure;
        omega_rng = [1 -1];
        A_rng = linspace(0.000002,400,40);
        B_rng = linspace(0.000001,200,40);
        for w =1:length(omega_rng)
            N_poss_0_t = zeros(length(A_rng),1);
            N_negs_0_t = zeros(length(A_rng),1);
            omega__c = omega_rng(w);
            for t = 1:length(A_rng)
                A = A_rng(t);
                B = B_rng(t);
                theta_cc_mul = A*(xi_m-0.5)+B;
                theta_cc_mul_nd = theta_cc_mul;
                %Claculating BCs
                C__posL2 = 1;
                C__posR2 = 10;
                    bcs_0 = -omega__c*theta_cc_mul_nd(1)/2 + sqrt(theta_cc_mul_nd(1)^2 + 4*C__posL2^2)/2;
                    bcs_1 = -omega__c*theta_cc_mul_nd(end)/2 + sqrt(theta_cc_mul_nd(end)^2 + 4*C__posR2^2)/2;
                    opts = bvpset('FJacobian',jac_call_theta_pe_var(A,B,omega__c,"pos",diff_cnst,i__val,Pe__val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
                     if omega__c == 1
                        solinit = bvpinit(xmesh,[10; -0.01]);
                     else
                        solinit = bvpinit(xmesh,[300; -0.01]);
                     end               
                    sol4c = bvp4c(bvp_call_mod_theta_pe_var(A,B,omega__c,"pos",diff_cnst,i__val,Pe__val), bcf_call(bcs_0,bcs_1),solinit, opts);
                    N_poss_0 = zeros(length(sol4c.x),1);
                    N_negs_0 = zeros(length(sol4c.x),1);
                    d__pos = diff_cnst;
                    for k = 1:length(sol4c.x)
                        C__1 = sol4c.y(1,k);
                        C__2 = sol4c.y(2,k);
                        xi_vl = sol4c.x(1,k);
                        xi = xi_vl;
                        if omega__c == 1
                        N_poss_0(k) = -(d__pos * C__2) + d__pos * (-omega__c * A + ((-1 + d__pos) * C__2) + omega__c * (A * (xi - 0.5e0) + B) * Pe__val + i__val) * C__1 / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) + C__1 * Pe__val;
                        else
                        N_negs_0(k) = (-d__pos * ((A * (xi - 0.5e0) + B) * omega__c + (2 * C__1)) * C__2 - A * C__1 * omega__c * d__pos + (Pe__val * (0.1e1 + d__pos) * C__1 - i__val) * (C__1 + (A * (xi - 0.5e0) + B) * omega__c)) / ((A * (xi - 0.5e0) + B) * omega__c + C__1 * (0.1e1 + d__pos));   
                        end
                    N_poss_0_t(t) = mean(N_poss_0);
                    N_negs_0_t(t) = mean(N_negs_0);
                    end
            end
                    if omega__c == 1
                    plot(A_rng,N_poss_0_t,'-o','LineWidth',3);hold on;
                    else
                    plot(A_rng,N_negs_0_t,'-o','LineWidth',3);hold on;
                    end
        end
                title("N^{*} vs \theta^* for AEM")
                legend("N_{pos}^* for AEM","N_{neg}^* for CEM",'location', 'best')
                xlabel('\theta^*')
                ylabel('N^{*}(\theta^*)')
                hold off;
                saveas(h(15),"plots_sv/f10_Npos_aem_cem_theta_var_inhm.png")
 %% For  Calculating positive Flux in AEM for inhomogeneous one side charge and discharge
        xmesh = linspace(0, 1,1000);
        xi_m = xmesh;
        diff_cnst = 0.5;
        i_rng = [-500,500];
        Pe__val = 0;
        omega_rng = [1 -1];
        A_rng = linspace(0.002,400,40);
        B_rng = linspace(0.001,200,40);
        for ii = 1:length(i_rng)
            h(15+ii) = figure;
            i__val = i_rng(ii);
            for w =1:length(omega_rng)
                N_poss_0_t = zeros(length(A_rng),1);
                omega__c = omega_rng(w);
                for t = 1:length(A_rng)
                    A = A_rng(t);
                    B = B_rng(t);
                    theta_cc_mul = A*(xi_m-0.5)+B;
                    theta_cc_mul_nd = theta_cc_mul;
                    %Claculating BCs
                    C__posL2 = 1;
                    C__posR2 = 10;
                        bcs_0 = -omega__c*theta_cc_mul_nd(1)/2 + sqrt(theta_cc_mul_nd(1)^2 + 4*C__posL2^2)/2;
                        bcs_1 = -omega__c*theta_cc_mul_nd(end)/2 + sqrt(theta_cc_mul_nd(end)^2 + 4*C__posR2^2)/2;
                        opts = bvpset('FJacobian',jac_call_theta_pe_var(A,B,omega__c,"pos",diff_cnst,i__val,Pe__val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
                        if i__val <=0
                            if omega__c == -1
                             solinit = bvpinit(xmesh,[300; 0]);
                            else
                             solinit = bvpinit(xmesh,[400; -0.01]);
                            end
                        else
                            if omega__c == -1
                             solinit = bvpinit(xmesh,[500; -0.01]);
                            else
                             solinit = bvpinit(xmesh,[100; -0.01]);
                            end                            
                        end
                        sol4c = bvp4c(bvp_call_mod_theta_pe_var(A,B,omega__c,"pos",diff_cnst,i__val,Pe__val), bcf_call(bcs_0,bcs_1),solinit, opts);
                        N_poss_0 = zeros(length(sol4c.x),1);
                        d__pos = diff_cnst;
                        for k = 1:length(sol4c.x)
                            C__1 = sol4c.y(1,k);
                            C__2 = sol4c.y(2,k);
                            xi_vl = sol4c.x(1,k);
                            xi = xi_vl;
                        N_poss_0(k) = -(d__pos * C__2) + d__pos * (-omega__c * A + ((-1 + d__pos) * C__2) + omega__c * (A * (xi - 0.5e0) + B) * Pe__val + i__val) * C__1 / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) + C__1 * Pe__val;
                        end
                 N_poss_0_t(t) = mean(N_poss_0);
                end
                        plot(A_rng,N_poss_0_t,'-o','LineWidth',3);hold on;
            end
                    title("N^{*} vs \theta^* for i = "+num2str(i__val))
                    legend("N_{pos}^* for AEM","N_{pos}^* for CEM",'location', 'best')
                    xlabel('\theta^*')
                    ylabel('N^{*}(\theta^*)')
                    hold off;
                    saveas(h(15+ii),"plots_sv/f10_Npos_aem_cem_theta_var_inhm_ivar"+num2str(i__val)+".png")
        end
%  %% For  Calculating positive Flux in AEM for inhomogeneous two side variation
%         xmesh = linspace(0, 1,100);
%         i__val = 0;
%         xi_m = xmesh;
%         diff_cnst = 0.5;
%         Pe__val = 0;
%         h(18) = figure;
%         omega_rng = [-1 1];
%         A_rng = linspace(0,200,40);
%         for w =1:length(omega_rng)
%             N_poss_0_t = zeros(length(A_rng),1);
%             N_negs_0_t = zeros(length(A_rng),1);
%             omega__c = omega_rng(w);
%             for t = 1:length(A_rng)
%                 A = A_rng(t);
%                 B = 0;
%                 A_1 = A_rng(t);
%                 theta_cc_mul = zeros(length(xmesh),1);
%             for i = 1:length(xmesh)
%                 if xi_m(i)<=0.5
%                 theta_cc_mul(i)= -A_1*(xi_m(i)-0.5);
%                 else
%                 theta_cc_mul(i) = A_1*(xi_m(i)-0.5);
%                 end
%             end
%             theta_cc_mul_nd = theta_cc_mul;
%                 %Claculating BCs
%                 C__posL2 = 1;
%                 C__posR2 = 10;
%                     bcs_0 = -omega__c*theta_cc_mul_nd(1)/2 + sqrt(theta_cc_mul_nd(1)^2 + 4*C__posL2^2)/2;
%                     bcs_1 = -omega__c*theta_cc_mul_nd(end)/2 + sqrt(theta_cc_mul_nd(end)^2 + 4*C__posR2^2)/2;
%                     opts = bvpset('FJacobian',jac_call_theta_pe_var(A,B,omega__c,"pos_inhm",diff_cnst,i__val,Pe__val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
%                      if omega__c == 1
%                         solinit = bvpinit(xmesh,[100; -0.01]);
%                      else
%                         solinit = bvpinit(xmesh,[300; -0.01]);
%                      end               
%                     sol4c = bvp4c(bvp_call_mod_theta_pe_var(A,B,omega__c,"pos_inhm",diff_cnst,i__val,Pe__val), bcf_call(bcs_0,bcs_1),solinit, opts);
%                     [sol4c.y(1,1),sol4c.y(2,1)]
%                     N_poss_0 = zeros(length(sol4c.x),1);
%                     N_negs_0 = zeros(length(sol4c.x),1);
%                     d__pos = diff_cnst;
%                     for k = 1:length(sol4c.x)
%                         C__1 = sol4c.y(1,k);
%                         C__2 = sol4c.y(2,k);
%                         xi_vl = sol4c.x(1,k);
%                         xi = xi_vl;
%                         if omega__c == 1
%                             if xi_vl<=0.5
%                                 A = A_1;
%                                 B = 0;
%                                 N_poss_0(k) = -(d__pos * C__2) + d__pos * (-omega__c * A + ((-1 + d__pos) * C__2) + omega__c * (A * (xi - 0.5e0) + B) * Pe__val + i__val) * C__1 / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) + C__1 * Pe__val;
%                             else
%                                 A = -A_1;
%                                 B = 0;
%                                 N_poss_0(k) = -(d__pos * C__2) + d__pos * (-omega__c * A + ((-1 + d__pos) * C__2) + omega__c * (A * (xi - 0.5e0) + B) * Pe__val + i__val) * C__1 / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) + C__1 * Pe__val;
%                             end
%                         else
%                             if xi_vl<=0.5
%                                 A = A_1;
%                                 B = 0;
%                                 N_negs_0(k) = (-d__pos * ((A * (xi - 0.5e0) + B) * omega__c + (2 * C__1)) * C__2 - A * C__1 * omega__c * d__pos + (Pe__val * (0.1e1 + d__pos) * C__1 - i__val) * (C__1 + (A * (xi - 0.5e0) + B) * omega__c)) / ((A * (xi - 0.5e0) + B) * omega__c + C__1 * (0.1e1 + d__pos));   
%                             else
%                                 A = -A_1;
%                                 B = 0;
%                                 N_negs_0(k) = (-d__pos * ((A * (xi - 0.5e0) + B) * omega__c + (2 * C__1)) * C__2 - A * C__1 * omega__c * d__pos + (Pe__val * (0.1e1 + d__pos) * C__1 - i__val) * (C__1 + (A * (xi - 0.5e0) + B) * omega__c)) / ((A * (xi - 0.5e0) + B) * omega__c + C__1 * (0.1e1 + d__pos));   
%                             end
%                         end
%                     N_poss_0_t(t) = mean(N_poss_0);
%                     N_negs_0_t(t) = mean(N_negs_0);
%                     end
%             end
%                     if omega__c == 1
%                     plot(A_rng,N_poss_0_t,'-o','LineWidth',3);hold on;
%                     else
%                     plot(A_rng,-N_negs_0_t,'-o','LineWidth',3);hold on;
%                     end
%         end
%                 title("N^{*} vs \theta^* for AEM")
%                 legend("N_{pos}^* for AEM","N_{neg}^* for CEM")
%                 xlabel('\theta^*')
%                 ylabel('N^{*}(\theta^*)')
%                 hold off;
%                 saveas(h(18),"plots_sv/f10_Npos_aem_cem_theta_var_inhm_2sid.png")
    %% For  Calculating positive Flux in AEM for varying Pe for Diff diff ratios      
        td__c = 0.01;
        dpos_rng = [0.1 0.5];
        Pe_rng = [-0.1,-1,-10,-100];
        i__val = 0;
        omega_rng = [1 -1];
       theta_cc_mul=linspace(0.01,4,80);
       theta_cc_mul_nd = theta_cc_mul/td__c;
        for ii = 1:length(Pe_rng)
            Pe__val = Pe_rng(ii);
            if Pe__val == -100
              xmesh = linspace(0,1,10000);
            else
                xmesh = linspace(0,1,1000);
            end
            h(18+ii) = figure;
             xi_m = xmesh;
            for dd = 1:length(dpos_rng)
                diff_cnst = dpos_rng(dd);
            for w =1:length(omega_rng)
                N_poss_0_t = zeros(length(theta_cc_mul_nd),1);
                N_negs_0_t = zeros(length(theta_cc_mul_nd),1);
                omega__c = omega_rng(w);
                for t = 1:length(theta_cc_mul_nd)
                    A = 0;
                    B = theta_cc_mul_nd(t);
                    %Claculating BCs
                    C__posL2 = 1;
                    C__posR2 = 10;
                        bcs_0 = -omega__c*theta_cc_mul_nd(t)/2 + sqrt(theta_cc_mul_nd(t)^2 + 4*C__posL2^2)/2;
                        bcs_1 = -omega__c*theta_cc_mul_nd(t)/2 + sqrt(theta_cc_mul_nd(t)^2 + 4*C__posR2^2)/2;
                        opts = bvpset('FJacobian',jac_call_theta_pe_var(A,B,omega__c,"pos",diff_cnst,i__val,Pe__val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
                        if i__val <=0
                            if omega__c == -1
                             solinit = bvpinit(xmesh,[300; 0]);
                            else
                             solinit = bvpinit(xmesh,[400; -0.01]);
                            end
                        else
                            if omega__c == -1
                             solinit = bvpinit(xmesh,[500; -0.01]);
                            else
                             solinit = bvpinit(xmesh,[100; -0.01]);
                            end                            
                        end
                        sol4c = bvp4c(bvp_call_mod_theta_pe_var(A,B,omega__c,"pos",diff_cnst,i__val,Pe__val), bcf_call(bcs_0,bcs_1),solinit, opts);
                        N_poss_0 = zeros(length(sol4c.x),1);
                        N_negs_0 = zeros(length(sol4c.x),1);
                        d__pos = diff_cnst;
                        for k = 1:length(sol4c.x)
                            C__1 = sol4c.y(1,k);
                            C__2 = sol4c.y(2,k);
                            xi_vl = sol4c.x(1,k);
                            xi = xi_vl;
                        if omega__c == 1
                        N_poss_0(k) = -(d__pos * C__2) + d__pos * (-omega__c * A + ((-1 + d__pos) * C__2) + omega__c * (A * (xi - 0.5e0) + B) * Pe__val + i__val) * C__1 / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) + C__1 * Pe__val;
                        else
                        N_negs_0(k) = (-d__pos * ((A * (xi - 0.5e0) + B) * omega__c + (2 * C__1)) * C__2 - A * C__1 * omega__c * d__pos + (Pe__val * (0.1e1 + d__pos) * C__1 - i__val) * (C__1 + (A * (xi - 0.5e0) + B) * omega__c)) / ((A * (xi - 0.5e0) + B) * omega__c + C__1 * (0.1e1 + d__pos));   
                        end
                        end
                    N_poss_0_t(t) = mean(N_poss_0);
                    N_negs_0_t(t) = mean(N_negs_0);
                end
                    if omega__c == 1
                    plot(theta_cc_mul_nd,N_poss_0_t,'-o','LineWidth',3);hold on;
                    else
                    plot(theta_cc_mul_nd,N_negs_0_t,'-o','LineWidth',3);hold on;
                    end
            end
            end
                    title("N^{*} vs \theta^* for Pe = "+num2str(-Pe__val))
                    legend("N_{pos}^* AEM (D_{+}=0.1)","N_{neg}^* CEM (D_{+}=0.1)","N_{pos}^* AEM (D_{+}=0.5)","N_{neg}^* CEM (D_{+}=0.5)",'location', 'best')
                    xlabel('\theta^*')
                    ylabel('N^{*}(\theta^*)')
                    hold off;
                    saveas(h(18+ii),"plots_sv/S3_N_aem_cem_Pevar"+num2str(Pe__val)+".png")
        end
  %% For  Calculating positive Flux in AEM for varying Pe for i =1000 charge
          td__c = 0.01;
        dpos_rng = 1;
        Pe_rng = [-10,-100];
        i__val = 1000;
        omega_rng = [1 -1];
       theta_cc_mul=linspace(0.01,4,80);
       theta_cc_mul_nd = theta_cc_mul/td__c;
        for ii = 1:length(Pe_rng)
            Pe__val = Pe_rng(ii);
            if Pe__val == -100
              xmesh = linspace(0,1,1000);
            else
                xmesh = linspace(0,1,1000);
            end
            h(21+ii) = figure;
             xi_m = xmesh;
            for dd = 1:length(dpos_rng)
                diff_cnst = dpos_rng(dd);
            for w =1:length(omega_rng)
                N_poss_0_t = zeros(length(theta_cc_mul_nd),1);
                N_negs_0_t = zeros(length(theta_cc_mul_nd),1);
                omega__c = omega_rng(w);
                if omega__c == 1
                    Pe__val = -Pe_rng(ii);
                else
                    Pe__val = Pe_rng(ii);
                end
                for t = 1:length(theta_cc_mul_nd)
                    A = 0;
                    B = theta_cc_mul_nd(t);
                    %Claculating BCs
                    C__posL2 = 1;
                    C__posR2 = 10;
                        bcs_0 = -omega__c*theta_cc_mul_nd(t)/2 + sqrt(theta_cc_mul_nd(t)^2 + 4*C__posL2^2)/2;
                        bcs_1 = -omega__c*theta_cc_mul_nd(t)/2 + sqrt(theta_cc_mul_nd(t)^2 + 4*C__posR2^2)/2;
                        opts = bvpset('FJacobian',jac_call_theta_pe_var(A,B,omega__c,"pos",diff_cnst,i__val,Pe__val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
                        if i__val <=0
                            if omega__c == -1
                             solinit = bvpinit(xmesh,[300; 0]);
                            else
                             solinit = bvpinit(xmesh,[400; -0.01]);
                            end
                        else
                            if omega__c == -1
                             solinit = bvpinit(xmesh,[500; -0.01]);
                            else
                             solinit = bvpinit(xmesh,[100; -0.01]);
                            end                            
                        end
                        sol4c = bvp4c(bvp_call_mod_theta_pe_var(A,B,omega__c,"pos",diff_cnst,i__val,Pe__val), bcf_call(bcs_0,bcs_1),solinit, opts);
                        N_poss_0 = zeros(length(sol4c.x),1);
                        N_negs_0 = zeros(length(sol4c.x),1);
                        d__pos = diff_cnst;
                        for k = 1:length(sol4c.x)
                            C__1 = sol4c.y(1,k);
                            C__2 = sol4c.y(2,k);
                            xi_vl = sol4c.x(1,k);
                            xi = xi_vl;
                        if omega__c == 1
                        N_poss_0(k) = -(d__pos * C__2) + d__pos * (-omega__c * A + ((-1 + d__pos) * C__2) + omega__c * (A * (xi - 0.5e0) + B) * Pe__val + i__val) * C__1 / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) + C__1 * Pe__val;
                        else
                        N_negs_0(k) = (-d__pos * ((A * (xi - 0.5e0) + B) * omega__c + (2 * C__1)) * C__2 - A * C__1 * omega__c * d__pos + (Pe__val * (0.1e1 + d__pos) * C__1 - i__val) * (C__1 + (A * (xi - 0.5e0) + B) * omega__c)) / ((A * (xi - 0.5e0) + B) * omega__c + C__1 * (0.1e1 + d__pos));   
                        end
                        end
                    N_poss_0_t(t) = mean(N_poss_0);
                    N_negs_0_t(t) = mean(N_negs_0);
                end
                    if omega__c == 1
                    plot(theta_cc_mul_nd,N_poss_0_t,'-o','LineWidth',3);hold on;
                    else
                    plot(theta_cc_mul_nd,N_negs_0_t,'-o','LineWidth',3);hold on;
                    end
            end
            end
                    title("N^{*} vs \theta^* for Pe = "+num2str(-Pe__val)+" for i ="+num2str(i__val))
                    legend("N_{pos}^* AEM ","N_{neg}^* CEM",'location', 'best')
                    xlabel('\theta^*')
                    ylabel('N^{*}(\theta^*)')
                    hold off;
                    saveas(h(21+ii),"plots_sv/S6_N_aem_cem_Pevar_ivar"+num2str(Pe__val)+".png")
        end
    %% For  Calculating positive Flux in AEM for varying Pe for i parametric sweep    
        td__c = 0.01;
        dpos_rng = 1;
        Pe_rng = [-0.0001,-1,-10,-100];
        i__rng = linspace(0,1000,20);
        omega_rng = [1 -1];
       theta_cc_mul=2;
       theta_cc_mul_nd = theta_cc_mul/td__c;
        for ii = 1:length(Pe_rng)
            Pe__val = Pe_rng(ii);
            if Pe__val == -100
              xmesh = linspace(0,1,10000);
            else
                xmesh = linspace(0,1,1000);
            end
            h(22+ii) = figure;
             xi_m = xmesh;
            for dd = 1:length(dpos_rng)
                diff_cnst = dpos_rng(dd);
            for w =1:length(omega_rng)
                N_poss_0_t = zeros(length(i__rng),1);
                N_negs_0_t = zeros(length(i__rng),1);
                omega__c = omega_rng(w);
                if omega__c == 1
                    Pe__val = -Pe_rng(ii);
                else
                    Pe__val = Pe_rng(ii);
                end
                for t = 1:length(i__rng)
                    i__val = i__rng(t);
                    A = 0;
                    B = theta_cc_mul_nd;
                    %Claculating BCs
                    C__posL2 = 1;
                    C__posR2 = 10;
                        bcs_0 = -omega__c*theta_cc_mul_nd/2 + sqrt(theta_cc_mul_nd^2 + 4*C__posL2^2)/2;
                        bcs_1 = -omega__c*theta_cc_mul_nd/2 + sqrt(theta_cc_mul_nd^2 + 4*C__posR2^2)/2;
                        opts = bvpset('FJacobian',jac_call_theta_pe_var(A,B,omega__c,"pos",diff_cnst,i__val,Pe__val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
                        if i__val <=0
                            if omega__c == -1
                             solinit = bvpinit(xmesh,[300; -0.01]);
                            else
                             solinit = bvpinit(xmesh,[400; -0.01]);
                            end
                        else
                            if omega__c == -1
                             solinit = bvpinit(xmesh,[500; -0.01]);
                            else
                             solinit = bvpinit(xmesh,[100; -0.01]);
                            end                            
                        end
                        sol4c = bvp4c(bvp_call_mod_theta_pe_var(A,B,omega__c,"pos",diff_cnst,i__val,Pe__val), bcf_call(bcs_0,bcs_1),solinit, opts);
                        N_poss_0 = zeros(length(sol4c.x),1);
                        N_negs_0 = zeros(length(sol4c.x),1);
                        d__pos = diff_cnst;
                        for k = 1:length(sol4c.x)
                            C__1 = sol4c.y(1,k);
                            C__2 = sol4c.y(2,k);
                            xi_vl = sol4c.x(1,k);
                            xi = xi_vl;
                        if omega__c == 1
                        N_poss_0(k) = -(d__pos * C__2) + d__pos * (-omega__c * A + ((-1 + d__pos) * C__2) + omega__c * (A * (xi - 0.5e0) + B) * Pe__val + i__val) * C__1 / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) + C__1 * Pe__val;
                        else
                        N_negs_0(k) = (-d__pos * ((A * (xi - 0.5e0) + B) * omega__c + (2 * C__1)) * C__2 - A * C__1 * omega__c * d__pos + (Pe__val * (0.1e1 + d__pos) * C__1 - i__val) * (C__1 + (A * (xi - 0.5e0) + B) * omega__c)) / ((A * (xi - 0.5e0) + B) * omega__c + C__1 * (0.1e1 + d__pos));   
                        end
                        end
                    N_poss_0_t(t) = mean(N_poss_0);
                    N_negs_0_t(t) = mean(N_negs_0);
                end
                    if omega__c == 1
                    plot(i__rng,N_poss_0_t,'-o','LineWidth',3);hold on;
                    else
                    plot(i__rng,N_negs_0_t,'-o','LineWidth',3);hold on;
                    end
            end
            end
                    title("N^{*} vs i^* for Pe = "+num2str(-Pe__val))
                    legend("N_{pos}^* AEM","N_{neg}^* CEM",'location', 'best')
                    xlabel('i^*')
                    ylabel('N^{*}')
                    hold off;
                    saveas(h(22+ii),"plots_sv/S7_N_aem_cem_Pevar_i_sweep"+num2str(Pe__val)+".png")
        end
  %% For  Calculating positive Flux in AEM for varying Pe for i =1000 charge & inhomogeneous modif on one side
        td__c = 0.01;
        xmesh = linspace(0,1,1000);
        dpos_rng = [0.5,1];
        Pe_rng = [-1,-10];
        omega_rng = [1 -1];
        A_rng = linspace(0.000002,400,40);
        B_rng = linspace(0.000001,200,40);
        for ii = 1:length(Pe_rng)
            Pe__val = Pe_rng(ii);
            if Pe__val == -10
              i__val = 100;
              diff_cnst = dpos_rng(ii);
            else
                i__val = 0;
                diff_cnst = dpos_rng(ii);
            end
            h(23+ii) = figure;
             xi_m = xmesh;
            for w =1:length(omega_rng)
                N_poss_0_t = zeros(length(A_rng),1);
                N_negs_0_t = zeros(length(A_rng),1);
                omega__c = omega_rng(w);
                if omega__c == 1 && Pe_rng(ii) == -10
                    Pe__val = -Pe_rng(ii);
                elseif omega__c == -1 && Pe_rng(ii) == -10
                    Pe__val = Pe_rng(ii);
                end
                for t = 1:length(A_rng)
                    A = A_rng(t);
                    B = B_rng(t);
                    theta_cc_mul = A*(xi_m-0.5)+B;
                    theta_cc_mul_nd = theta_cc_mul;
                    %Claculating BCs
                    C__posL2 = 1;
                    C__posR2 = 10;
                        bcs_0 = -omega__c*theta_cc_mul_nd(1)/2 + sqrt(theta_cc_mul_nd(1)^2 + 4*C__posL2^2)/2;
                        bcs_1 = -omega__c*theta_cc_mul_nd(end)/2 + sqrt(theta_cc_mul_nd(end)^2 + 4*C__posR2^2)/2;
                        opts = bvpset('FJacobian',jac_call_theta_pe_var(A,B,omega__c,"pos",diff_cnst,i__val,Pe__val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
                        if i__val <=0
                            if omega__c == 1
                             solinit = bvpinit(xmesh,[10; -0.01]);
                            else
                             solinit = bvpinit(xmesh,[300; -0.01]);
                            end
                        else
                            if omega__c == -1
                             solinit = bvpinit(xmesh,[500; -0.01]);
                            else
                             solinit = bvpinit(xmesh,[100; -0.01]);
                            end                            
                        end
                        sol4c = bvp4c(bvp_call_mod_theta_pe_var(A,B,omega__c,"pos",diff_cnst,i__val,Pe__val), bcf_call(bcs_0,bcs_1),solinit, opts);
                        N_poss_0 = zeros(length(sol4c.x),1);
                        N_negs_0 = zeros(length(sol4c.x),1);
                        d__pos = diff_cnst;
                        for k = 1:length(sol4c.x)
                            C__1 = sol4c.y(1,k);
                            C__2 = sol4c.y(2,k);
                            xi_vl = sol4c.x(1,k);
                            xi = xi_vl;
                        if omega__c == 1
                        N_poss_0(k) = -(d__pos * C__2) + d__pos * (-omega__c * A + ((-1 + d__pos) * C__2) + omega__c * (A * (xi - 0.5e0) + B) * Pe__val + i__val) * C__1 / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) + C__1 * Pe__val;
                        else
                        N_negs_0(k) = (-d__pos * ((A * (xi - 0.5e0) + B) * omega__c + (2 * C__1)) * C__2 - A * C__1 * omega__c * d__pos + (Pe__val * (0.1e1 + d__pos) * C__1 - i__val) * (C__1 + (A * (xi - 0.5e0) + B) * omega__c)) / ((A * (xi - 0.5e0) + B) * omega__c + C__1 * (0.1e1 + d__pos));   
                        end
                        end
                    N_poss_0_t(t) = mean(N_poss_0);
                    N_negs_0_t(t) = mean(N_negs_0);
                end
                    if omega__c == 1
                    plot(A_rng,N_poss_0_t,'-o','LineWidth',3);hold on;
                    else
                    plot(A_rng,N_negs_0_t,'-o','LineWidth',3);hold on;
                    end
            end
                    title("N^{*} vs \theta^* for Pe = "+num2str(-Pe__val)+" for i ="+num2str(i__val)+" inhomo 1 side")
                    legend("N_{pos}^* AEM ","N_{neg}^* CEM",'location', 'best')
                    xlabel('\theta^*')
                    ylabel('N^{*}(\theta^*)')
                    hold off;
                    saveas(h(23+ii),"plots_sv/S8_N_aem_cem_Pevar_ivar_inhom_1_side"+num2str(Pe__val)+".png")
        end
  %% For  Calculating positive Flux in AEM for varying Pe for i =1000 charge & inhomogeneous modif on both side
        td__c = 0.01;
        xmesh = linspace(0,1,1000);
        dpos_rng = [0.5];
        Pe_rng = [-10];
        omega_rng = [-1,1];
        A_rng = linspace(0.000002,400,40);
        B_rng = 0;
        for ii = 1:length(Pe_rng)
            Pe__val = Pe_rng(ii);
            if Pe__val == -10
              i__val = 100;
            else
                i__val = 0;
            end
            h(24+ii) = figure;
             xi_m = xmesh;
            for dd = 1:length(dpos_rng)
                diff_cnst = dpos_rng(dd);
            for w =1:length(omega_rng)
                N_poss_0_t = zeros(length(A_rng),1);
                N_negs_0_t = zeros(length(A_rng),1);
                omega__c = omega_rng(w);
                if omega__c == 1 && Pe_rng(ii) == -10
                    Pe__val = -Pe_rng(ii);
                elseif omega__c == -1 && Pe_rng(ii) == -10
                    Pe__val = Pe_rng(ii);
                end
                for t = 1:length(A_rng)
                    A = A_rng(t);
                    B = 0;
                    for i = 1:length(xmesh)
                        if xi_m(i)<=0.5
                        theta_cc_mul(i)= -A_1*(xi_m(i)-0.5);
                        else
                        theta_cc_mul(i) = A_1*(xi_m(i)-0.5);
                        end
                    end
                theta_cc_mul_nd = theta_cc_mul;
                    %Claculating BCs
                    C__posL2 = 1;
                    C__posR2 = 10;
                        bcs_0 = -omega__c*theta_cc_mul_nd(1)/2 + sqrt(theta_cc_mul_nd(1)^2 + 4*C__posL2^2)/2;
                        bcs_1 = -omega__c*theta_cc_mul_nd(end)/2 + sqrt(theta_cc_mul_nd(end)^2 + 4*C__posR2^2)/2;
                        opts = bvpset('FJacobian',jac_call_theta_pe_var(A,B,omega__c,"pos_inhm",diff_cnst,i__val,Pe__val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
                        if i__val <=0
                            if omega__c == 1
                             solinit = bvpinit(xmesh,[100; -0.01]);
                            else
                             solinit = bvpinit(xmesh,[300; -0.01]);
                            end
                        else
                            if omega__c == -1
                             solinit = bvpinit(xmesh,[500; -0.01]);
                            else
                             solinit = bvpinit(xmesh,[100; -0.01]);
                            end                            
                        end
                        sol4c = bvp4c(bvp_call_mod_theta_pe_var(A,B,omega__c,"pos_inhm",diff_cnst,i__val,Pe__val), bcf_call(bcs_0,bcs_1),solinit, opts);
                        N_poss_0 = zeros(length(sol4c.x),1);
                        N_negs_0 = zeros(length(sol4c.x),1);
                        d__pos = diff_cnst;
                        for k = 1:length(sol4c.x)
                            C__1 = sol4c.y(1,k);
                            C__2 = sol4c.y(2,k);
                            xi_vl = sol4c.x(1,k);
                            xi = xi_vl;
                        if omega__c == 1
                            if xi_vl<=0.5
                                A = A_1;
                                B = 0;
                                N_poss_0(k) = -(d__pos * C__2) + d__pos * (-omega__c * A + ((-1 + d__pos) * C__2) + omega__c * (A * (xi - 0.5e0) + B) * Pe__val + i__val) * C__1 / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) + C__1 * Pe__val;
                            else
                                A = -A_1;
                                B = 0;
                                N_poss_0(k) = -(d__pos * C__2) + d__pos * (-omega__c * A + ((-1 + d__pos) * C__2) + omega__c * (A * (xi - 0.5e0) + B) * Pe__val + i__val) * C__1 / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) + C__1 * Pe__val;
                            end
                        else
                            if xi_vl<=0.5
                                A = A_1;
                                B = 0;
                                N_negs_0(k) = (-d__pos * ((A * (xi - 0.5e0) + B) * omega__c + (2 * C__1)) * C__2 - A * C__1 * omega__c * d__pos + (Pe__val * (0.1e1 + d__pos) * C__1 - i__val) * (C__1 + (A * (xi - 0.5e0) + B) * omega__c)) / ((A * (xi - 0.5e0) + B) * omega__c + C__1 * (0.1e1 + d__pos));   
                            else
                                A = -A_1;
                                B = 0;
                                N_negs_0(k) = (-d__pos * ((A * (xi - 0.5e0) + B) * omega__c + (2 * C__1)) * C__2 - A * C__1 * omega__c * d__pos + (Pe__val * (0.1e1 + d__pos) * C__1 - i__val) * (C__1 + (A * (xi - 0.5e0) + B) * omega__c)) / ((A * (xi - 0.5e0) + B) * omega__c + C__1 * (0.1e1 + d__pos));   
                            end
                        end
                        end
                    N_poss_0_t(t) = mean(N_poss_0);
                    N_negs_0_t(t) = mean(N_negs_0);
                 end
            end
                    if omega__c == 1
                    plot(A_rng,N_poss_0_t,'-o','LineWidth',3);hold on;
                    else
                    plot(A_rng,-N_negs_0_t,'-o','LineWidth',3);hold on;
                    end
         end
                    title("N^{*} vs \theta^* for Pe = "+num2str(-Pe__val)+" for i ="+num2str(i__val)+" inhomo 1 side")
                    legend("N_{pos}^* AEM ","N_{neg}^* CEM",'location', 'best')
                    xlabel('\theta^*')
                    ylabel('N^{*}(\theta^*)')
                    hold off;
                    saveas(h(24+ii),"plots_sv/S8_N_aem_cem_Pevar_ivar_inhom_2_side"+num2str(Pe__val)+".png")
        end        

%% Figure Hiding
set(0,'DefaultFigureVisible','on');
%% figure off and variable clearance
clc
clear variables
h = [];
%% Initial Conditions
%sol5c = bvp5c(@bvpfcn, @bcfcn, solinit, opts);
td__c = 0.01;
theta_cc_mul=[0.01,0.1,1,2];
theta_cc_mul_nd = theta_cc_mul/td__c;
col_dat = ["red", "blue", "green", "magenta"];
%% for positive & Negative AEM
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
    legend("\theta^* = 1","\theta^* = 10","\theta^* = 100","\theta^* = 200")
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
    legend("\theta^* = 1","\theta^* = 10","\theta^* = 100","\theta^* = 200")
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
    bcs_0(i) = theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posL2^2)/2;
    bcs_1(i) = theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posR2^2)/2;
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
    legend("\theta^* = 1","\theta^* = 10","\theta^* = 100","\theta^* = 200")
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
    legend("\theta^* = 1","\theta^* = 10","\theta^* = 100","\theta^* = 200")
    xlabel('\xi')
    ylabel('C_{neg}^{*}(\xi)')
    hold off;
    saveas(h(4),"plots_sv/f3_cneg_cem.png")
    fprintf(2,"Plot: "+ num2str(4)+" Done\n");
%% For  Calculating positive Flux in AEM & negative flux in CEM
% For positive Flux in AEM
shp = '-o';
diff_cons = 0.5;
Pe__val = 0.1;
i__val = 0;
theta_cc_mul=linspace(0.01,4,80);
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
    bcs_0(i) = -theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posL2^2)/2;
    bcs_1(i) = -theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posR2^2)/2;
end
for i=1:length(theta_cc_mul)
    xmesh = linspace(0, 1, 100);
    opts = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",diff_cons,i__val,Pe__val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
    solinit = bvpinit(xmesh,[0; 1]);
    sol4c = bvp4c(bvp_call_mod_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",diff_cons,i__val,Pe__val), bcf_call(bcs_0(i),bcs_1(i)),solinit, opts);
    N_poss_0 = zeros(length(sol4c.x),1);
    d__pos = diff_cons;
    A = 0;
    B = theta_cc_mul_nd(i);
    for k = 1:length(sol4c.x)
        C__1 = sol4c.y(1,k);
        C__2 = sol4c.y(2,k);
        xi = sol4c.x(1,k);
      N_poss_0 = -(d__pos * C__2) + d__pos * (-omega__c * A + ((-1 + d__pos) * C__2) + omega__c * (A * (xi - 0.5e0) + B) * Pe__val + i__val) * C__1 / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) + C__1 * Pe__val;
    end
    N_poss(i) = mean(N_poss_0);
end
plot(theta_cc_mul_nd,N_poss,shp,'LineWidth',3)
title("N_{pos}^{*} vs \theta^* for AEM")
    xlabel('\theta')
    ylabel('N_{pos}^{*}(\theta)')
    hold off;
    saveas(h(5),"plots_sv/f3_Npos_aem.png")
    fprintf(2,"Plot: "+ num2str(5)+" Done\n");
%% For Calculating negative Flux in CEM
h(6) = figure;
theta_cc_mul=linspace(0.01,4,800);
theta_cc_mul_nd = theta_cc_mul/td__c;
omega__c = -1;
diff_cons = 0.5;
%Claculating BCs
C__posL2 = 1;
C__posR2 = 10;
Pe_val = -1;
bcs_0 = [];
bcs_1 = [];
sol4c_yd = [];
for i = 1:length(theta_cc_mul)
    bcs_0(i) = theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posL2^2)/2;
    bcs_1(i) = theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posR2^2)/2;
end
for i=1:length(theta_cc_mul)
    xmesh = linspace(0, 1, 1000);
    opts = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",diff_cons,0,Pe_val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
    solinit = bvpinit(xmesh,[300; -0]);
    sol4c = bvp4c(bvp_call_mod_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",diff_cons,0,Pe_val), bcf_call(bcs_0(i),bcs_1(i)),solinit, opts);
    [sol4c.y(1,1),sol4c.y(2,1)]
    N_negs_0 = zeros(length(sol4c.x),1);
    A = 0;
    B = theta_cc_mul_nd(i);
    d_pos = diff_cons;
    for k = 1:length(sol4c.x)
        C__1 = sol4c.y(1,k);
        C__2 = sol4c.y(2,k);
        xi = sol4c.x(1,k);
        N_negs_0 = (-d__pos * ((A * (xi - 0.5e0) + B) * omega__c + (2 * C__1)) * C__2 - A * C__1 * omega__c * d__pos + (Pe__val * (0.1e1 + d__pos) * C__1 - i__val) * (C__1 + (A * (xi - 0.5e0) + B) * omega__c)) / ((A * (xi - 0.5e0) + B) * omega__c + C__1 * (0.1e1 + d__pos));
    end
    N_negs(i) = mean(N_negs_0);
end
plot(theta_cc_mul_nd,N_negs,'-o','LineWidth',3)
title("N_{neg}^{*} vs \theta^* for AEM")
    xlabel('\theta')
    ylabel('N_{neg}^{*}(\theta)')
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
    legend("N_{pos} for AEM","N_{neg} for CEM")
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
    xmesh = linspace(0, 1, 10000);
    opts = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_const,omega__c,"pos",diff_cnst(i),i_const,Pe__val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
    solinit = bvpinit(xmesh,[0;0.15]);
    sol4c = bvp4c(bvp_call_mod_theta_pe_var(0,theta_const,omega__c,"pos",diff_cnst(i),i_const,Pe__val), bcf_call(bcs_0,bcs_1),solinit, opts);
    plot(sol4c.x,sol4c.y(1,:),'LineWidth',3,'Color',col_dat(i));
    hold on;
end
    title("C_{pos}^{*}(\xi) vs \xi in AEM for different ratios of the diffusion coefficient")
    legend("D_{+}/D_{-}=0.025","D_{+}/D_{-}=1","D_{+}/D_{-}=10","D_{+}/D_{-}=40")%,"D_{+}/D_{-}=1","D_{+}/D_{-}=8")
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
    legend("D_{+}/D_{-}=0.025","D_{+}/D_{-}=1","D_{+}/D_{-}=10","D_{+}/D_{-}=40")
    xlabel('\xi')
    ylabel('C_{pos}^{*}(\xi)')
    hold off;
    saveas(h(9),"plots_sv/f4_diff_cpos_cem.png")
    fprintf(2,"Plot: "+ num2str(9)+" Done\n");
%% For Calculating positive Flux in CEM & AEM Charge
%for i =500,-500,1000,-1000 for CEM
i_act = [500,800]; % Actual Current Density
i_nond = 1;
i_dens = i_act/i_nond; % Dimensionless Current Density
% For  Calculating positive Flux in CEM
theta_cc_mul=linspace(0.01,4,20);
theta_cc_mul_nd = theta_cc_mul/td__c;
h(10) = figure;
diff_cons = 0.5;
Pe__val = 0;
omega__c_rng = [-1,1];
%Claculating BCs
C__posL2 = 1;
C__posR2 = 10;
for ii = 1:2
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
        solinit_0 = bvpinit(xmesh,[0; 1]);
        sol4c_0 = bvp4c(bvp_call_mod_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",diff_cons,0,Pe__val), bcf_call(bcs_0(i),bcs_1(i)),solinit_0, opts_0);
        if i_dens(j)>=0
        solinit = bvpinit(xmesh,[sol4c_0.y(1,1); sol4c_0.y(2,1)]);
        else
        solinit = bvpinit(xmesh,[-sol4c_0.y(1,1); -sol4c_0.y(2,1)]);
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
            N_poss_cem_1_0 = -(d__pos * C__2) + d__pos * (-omega__c * A + ((-1 + d__pos) * C__2) + omega__c * (A * (xi - 0.5e0) + B) * Pe__val + i__val) * C__1 / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) + C__1 * Pe__val;
        end
        N_poss_cem_1(i) = mean(N_poss_cem_1_0);
    end
plot(theta_cc_mul_nd,N_poss_cem_1,'-o','LineWidth',3);hold on;
end
end
title("Cation Fluxes for various current densities in CEM")
legend("N_{pos}^* for CEM(500)","N_{pos}^*  for CEM(1000)","N_{pos}^*  for AEM(500)","N_{pos}^*  for AEM(1000)")
    xlabel('\theta')
    ylabel('N_{pos}^{*}(\theta^*)')
saveas(h(10),"plots_sv/f7_pos_i_flux_cem.png")
fprintf(2,"Plot: "+ num2str(10)+" Done\n");
hold off;
%% For Calculating positive Flux in AEM
%for i =500,-500,1000,-1000 for AEM
i_act = [-1000,-500,500,1000]; % Actual Current Density
i_nond = 100;
i_dens = i_act/i_nond; % Dimensionless Current Density
% For  Calculating positive Flux in AEM
theta_cc_mul=linspace(0.01,4,20);
theta_cc_mul_nd = theta_cc_mul/td__c;
h(11) = figure;
diff_cons = 0.5;
Pe__Val = 0;
omega__c = 1;
%Claculating BCs
C__posL2 = 1;
C__posR2 = 10;
for j = 1:length(i_act)
    N_poss_aem_1 = zeros(length(theta_cc_mul),1);
    bcs_0 = zeros(length(theta_cc_mul),1);
    bcs_1 = zeros(length(theta_cc_mul),1);
    for i = 1:length(theta_cc_mul)
        bcs_0(i) = -theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posL2^2)/2;
        bcs_1(i) = -theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posR2^2)/2;
    end
    for i=1:length(theta_cc_mul)
        xmesh = linspace(0, 1, 10000);
        opts_0 = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",diff_cons,0,Pe__Val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
        solinit_0 = bvpinit(xmesh,[0; 1]);
        sol4c_0 = bvp4c(bvp_call_mod_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",diff_cons,0,Pe__Val), bcf_call(bcs_0(i),bcs_1(i)),solinit_0, opts_0);
        if i_dens(j)<=0
        solinit = bvpinit(xmesh,[sol4c_0.y(1,1); sol4c_0.y(2,1)]);
        else
        solinit = bvpinit(xmesh,[sol4c_0.y(1,1); sol4c_0.y(2,1)]);
        end
        opts = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",diff_cons,i_dens(j),Pe__Val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
        sol4c = bvp4c(bvp_call_mod_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",diff_cons,i_dens(j),Pe__Val), bcf_call(bcs_0(i),bcs_1(i)),solinit, opts);
        N_poss_aem_1_0 =[];
        i__val = i_dens(j);
        theta__c = theta_cc_mul_nd(i);
        A__1 = 0;
        B__1 = theta_cc_mul_nd(i);
        i__Val = i_dens(j);
        d__pos = diff_cons;
        for k = 1:length(sol4c.x)
            C__1 = sol4c.y(1,k);
            C__2 = sol4c.y(2,k);
            N_poss_aem_1_0(k) = -C__2*d__pos - d__pos*C__1*((1 - d__pos)*C__2 - omega__c*(A__1*(xi - 0.5) + B__1)*Pe__Val + omega__c*A__1 - i__Val)/((A__1*(xi - 0.5) + B__1)*omega__c + C__1*(1 + d__pos)) + C__1*Pe__Val;
        end
        N_poss_aem_1(i) = mean(N_poss_aem_1_0);
    end
plot(theta_cc_mul,N_poss_aem_1,'-o','LineWidth',3);hold on;
end
title("Cation Fluxes for various current densities in AEM")
legend("N_{pos} for AEM_(-1000)","N_{pos} for AEM_(-500)","N_{pos} for AEM_(500)","N_{pos} for AEM_(1000)")
    xlabel('\theta')
    ylabel('N_{pos}^{*}(\theta)')
hold off;
saveas(h(11),"plots_sv/f7_pos_i_flux_aem.png")
fprintf(2,"Plot: "+ num2str(11)+" Done\n");
%% Diffusion coeff variation for pos for AEM @ theta = 2 for i charge
h(12,1) = figure;
Pe__Val = 0;
omega__c = 1;
theta_const = 2/td__c;
diff_cnst = [0.025 1 10 40];
i_orig = 500;
i_nond = 1;
i_const  = i_orig/i_nond;
%Claculating BCs
C__posL2 = 1;
C__posR2 = 10;
    bcs_0 = -omega__c*theta_const/2 + sqrt(theta_const^2 + 4*C__posL2^2)/2;
    bcs_1 = -omega__c*theta_const/2 + sqrt(theta_const^2 + 4*C__posR2^2)/2;
for i=1:length(diff_cnst)
    xmesh = linspace(0, 1, 10000);
    opts = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_const,omega__c,"pos",diff_cnst(i),i_const,Pe__Val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
    solinit = bvpinit(xmesh,[0.995;0.5348]);
    sol4c = bvp4c(bvp_call_mod_theta_pe_var(0,theta_const,omega__c,"pos",diff_cnst(i),i_const,Pe__Val), bcf_call(bcs_0,bcs_1),solinit, opts);
    plot(sol4c.x,sol4c.y(1,:),'LineWidth',3,'Color',col_dat(i));
    hold on;
end
title("C_{pos}^{*}(\xi) vs \xi in AEM @ i= "+num2str(i_orig)+" for various D_{+}/D_{-}")
    legend("D_{+}/D_{-}=0.025","D_{+}/D_{-}=1","D_{+}/D_{-}=10","D_{+}/D_{-}=40")%,"D_{+}/D_{-}=1","D_{+}/D_{-}=8")
    xlabel('\xi')
    ylabel('C_{pos}^{*}(\xi)')
    hold off;
    saveas(h(12,1),"plots_sv/f8_diff_cpos_aem_i_ch.png")
    fprintf(2,"Plot: "+ num2str(12.1)+" Done\n");
%% Diffusion coeff variation for pos for AEM @ theta = 2 for i discharge
h(12,2) = figure;
Pe__Val = 0;
omega__c = 1;
theta_const = 2/td__c;
diff_cnst = [0.025 1 10 40];
i_orig = -500;
i_nond = 1;
i_const  = i_orig/i_nond;
%Claculating BCs
C__posL2 = 1;
C__posR2 = 10;
    bcs_0 = -omega__c*theta_const/2 + sqrt(theta_const^2 + 4*C__posL2^2)/2;
    bcs_1 = -omega__c*theta_const/2 + sqrt(theta_const^2 + 4*C__posR2^2)/2;
for i=1:length(diff_cnst)
    xmesh = linspace(0, 1, 10000);
    opts = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_const,omega__c,"pos",diff_cnst(i),i_const,Pe__Val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
    solinit = bvpinit(xmesh,[0.995;0.5348]);
    sol4c = bvp4c(bvp_call_mod_theta_pe_var(0,theta_const,omega__c,"pos",diff_cnst(i),i_const,Pe__Val), bcf_call(bcs_0,bcs_1),solinit, opts);
    plot(sol4c.x,sol4c.y(1,:),'LineWidth',3,'Color',col_dat(i));
    hold on;
end
title("C_{pos}^{*}(\xi) vs \xi in AEM @ i= "+num2str(i_orig)+" for various D_{+}/D_{-}")
    legend("D_{+}/D_{-}=0.025","D_{+}/D_{-}=1","D_{+}/D_{-}=10","D_{+}/D_{-}=40")%,"D_{+}/D_{-}=1","D_{+}/D_{-}=8")
    xlabel('\xi')
    ylabel('C_{pos}^{*}(\xi)')
    hold off;
    saveas(h(12,2),"plots_sv/f8_diff_cpos_aem_i_dch.png")
    fprintf(2,"Plot: "+ num2str(12.2)+" Done\n");
%% Diffusion coeff variation for pos for CEM @ theta = 2 for i charge
h(13,1) = figure;
Pe__Val = 0;
omega__c = -1;
theta_const = 2/td__c;
i_orig = 500;
i_nond = 1;
i_const  = i_orig/i_nond;
diff_cnst = [0.025 1 10 40];
%Claculating BCs
C__posL2 = 1;
C__posR2 = 10;
    bcs_0 = -omega__c*theta_const/2 + sqrt(theta_const^2 + 4*C__posL2^2)/2;
    bcs_1 = -omega__c*theta_const/2 + sqrt(theta_const^2 + 4*C__posR2^2)/2;
for i=1:length(diff_cnst)
    xmesh = linspace(0, 1, 1000);
    opts = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_const,omega__c,"pos",diff_cnst(i),i_const,Pe__Val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
    solinit = bvpinit(xmesh,[200.005; 0.4938]);
    sol4c = bvp4c(bvp_call_mod_theta_pe_var(0,theta_const,omega__c,"pos",diff_cnst(i),i_const,Pe__Val), bcf_call(bcs_0,bcs_1),solinit, opts);
    plot(sol4c.x,sol4c.y(1,:),'LineWidth',3,'Color',col_dat(i));
    hold on;
end
title("C_{pos}^{*}(\xi) vs \xi in CEM @ i= "+num2str(i_orig)+" for various D_{+}/D_{-}")
    legend("D_{+}/D_{-}=0.025","D_{+}/D_{-}=1","D_{+}/D_{-}=10","D_{+}/D_{-}=40")
    xlabel('\xi')
    ylabel('C_{pos}^{*}(\xi)')
    hold off;
    saveas(h(13,1),"plots_sv/f8_diff_cpos_cem_i_ch.png")
    fprintf(2,"Plot: "+ num2str(13.1)+" Done\n");
%% Diffusion coeff variation for pos for CEM @ theta = 2 for i discharge
h(13,2) = figure;
Pe__Val = 0;
omega__c = -1;
theta_const = 2/td__c;
i_orig = -500;
i_nond = 1;
i_const  = i_orig/i_nond;
diff_cnst = [0.025 1 10 40]*100;
%Claculating BCs
C__posL2 = 1;
C__posR2 = 10;
    bcs_0 = -omega__c*theta_const/2 + sqrt(theta_const^2 + 4*C__posL2^2)/2;
    bcs_1 = -omega__c*theta_const/2 + sqrt(theta_const^2 + 4*C__posR2^2)/2;
for i=1:length(diff_cnst)
    xmesh = linspace(0, 1, 1000);
    opts = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_const,omega__c,"pos",diff_cnst(i),i_const,Pe__Val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
    %opts = bvpset('FJacobian',jac_call(theta_const,omega__c,"pos",diff_cnst(i),i_const),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
    solinit = bvpinit(xmesh,[200; 0]);
    sol4c = bvp4c(bvp_call_mod_theta_pe_var(0,theta_const,omega__c,"pos",diff_cnst(i),i_const,Pe__Val), bcf_call(bcs_0,bcs_1),solinit, opts);
    %sol4c = bvp4c(bvp_call_mod(theta_const,omega__c,"pos",diff_cnst(i),i_const), bcf_call(bcs_0(i),bcs_1(i)),solinit, opts);
    plot(sol4c.x,sol4c.y(1,:),'LineWidth',3,'Color',col_dat(i));
    hold on;
end
title("C_{pos}^{*}(\xi) vs \xi in CEM @ i= "+num2str(i_orig)+" for various D_{+}/D_{-}")
    legend("D_{+}/D_{-}=0.025","D_{+}/D_{-}=1","D_{+}/D_{-}=10","D_{+}/D_{-}=40")
    xlabel('\xi')
    ylabel('C_{pos}^{*}(\xi)')
    hold off;
    saveas(h(13,2),"plots_sv/f8_diff_cpos_cem_i_dch.png")
    fprintf(2,"Plot: "+ num2str(13.2)+" Done\n");
%% For co-ion and counter ion Fluxes in CEM as function of current density
i_act = linspace(-1000,1000,10); % Actual Current Density
i_dens = i_act; % Dimensionless Current Density
xmesh = linspace(0, 1, 10000);
Pe__val = 0;
diff_cons = 0.5;
% For  Calculating positive Flux in CEM
theta_cc_mul=2;
theta_cc_mul_nd = theta_cc_mul/td__c;
h(14) = figure;
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
    solinit = bvpinit(xmesh,[-sol4c_0.y(1,1); -sol4c_0.y(2,1)]);
    end
    opts = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_cc_mul_nd,omega__c,"pos",diff_cons,i_dens(j),Pe__Val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
    sol4c = bvp4c(bvp_call_mod_theta_pe_var(0,theta_cc_mul_nd,omega__c,"pos",diff_cons,i_dens(j),Pe__Val), bcf_call(bcs_0,bcs_1),solinit, opts);
    N_poss_cem_1_0 = zeros(length(sol4c.x),1);
    N_negs_cem_1_0 = zeros(length(sol4c.x),1);
    i__val = i_dens(j);
    d__pos = diff_cons;
    A__1 = 0;
    B__1 = theta_cc_mul_nd;
    for k = 1:length(sol4c.x)
        C__1 = sol4c.y(1,k);
        C__2 = sol4c.y(2,k);
        N_poss_cem_1_0(k) = -(d__pos * C__2) + d__pos * (-omega__c * A + ((-1 + d__pos) * C__2) + omega__c * (A * (xi - 0.5e0) + B) * Pe__val + i__val) * C__1 / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) + C__1 * Pe__val;
        N_negs_cem_1_0(k) = (-d__pos * ((A * (xi - 0.5e0) + B) * omega__c + (2 * C__1)) * C__2 - A * C__1 * omega__c * d__pos + (Pe__val * (0.1e1 + d__pos) * C__1 - i__val) * (C__1 + (A * (xi - 0.5e0) + B) * omega__c)) / ((A * (xi - 0.5e0) + B) * omega__c + C__1 * (0.1e1 + d__pos));
    end
N_poss_cem_1(j) = mean(N_poss_cem_1_0);
N_negs_cem_1(j) = mean(N_negs_cem_1_0);
end
plot(i_act,N_poss_cem_1,'-o','LineWidth',3,'Color',"red");hold on;
plot(i_act,N_negs_cem_1,'-o','LineWidth',3,'Color',"yellow");
 title("Co-ion and Counter ion flux as func of i")
legend("N_{pos} for CEM","N_{neg} for CEM")
    xlabel('i(cuurent density)')
    ylabel('N_{i}^{*}')
hold off;
saveas(h(14),"plots_sv/f7_pos_i_var_flux_cem.png")
fprintf(2,"Plot: "+ num2str(14)+" Done\n");
%% For co-ion and counter ion Fluxes in AEM as function of current density
h(15) = figure;
i_act = linspace(-1000,1000,10); % Actual Current Density
i_dens = i_act; % Dimensionless Current Density
xmesh = linspace(0, 1, 100);
Pe__val = 0;
diff_cons = 0.5;
% For  Calculating positive Flux in CEM
theta_cc_mul=2;
theta_cc_mul_nd = theta_cc_mul/td__c;
h(14) = figure;
omega__c = 1;
%Claculating BCs
C__posL2 = 1;
C__posR2 = 10;
N_poss_cem_1 = zeros(length(i_act),1);
N_negs_cem_1 = zeros(length(i_act),1);
for j = 1:length(i_act)
        bcs_0 = theta_cc_mul_nd/2 + sqrt(theta_cc_mul_nd^2 + 4*C__posL2^2)/2;
        bcs_1 = theta_cc_mul_nd/2 + sqrt(theta_cc_mul_nd^2 + 4*C__posR2^2)/2;
    opts_0 = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_cc_mul_nd,omega__c,"pos",diff_cons,0,Pe__val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
    solinit_0 = bvpinit(xmesh,[0; 1]);
    sol4c_0 = bvp4c(bvp_call_mod_theta_pe_var(0,theta_cc_mul_nd,omega__c,"pos",diff_cons,0,Pe__val), bcf_call(bcs_0,bcs_1),solinit_0, opts_0);
    if i_dens(j)>=0
    solinit = bvpinit(xmesh,[sol4c_0.y(1,1); sol4c_0.y(2,1)]);
    else
    solinit = bvpinit(xmesh,[-sol4c_0.y(1,1); -sol4c_0.y(2,1)]);
    end
    opts = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_cc_mul_nd,omega__c,"pos",diff_cons,i_dens(j),Pe__Val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
    sol4c = bvp4c(bvp_call_mod_theta_pe_var(0,theta_cc_mul_nd,omega__c,"pos",diff_cons,i_dens(j),Pe__Val), bcf_call(bcs_0,bcs_1),solinit, opts);
    N_poss_cem_1_0 = zeros(length(sol4c.x),1);
    N_negs_cem_1_0 = zeros(length(sol4c.x),1);
    i__val = i_dens(j);
    d__pos = diff_cons;
    A__1 = 0;
    B__1 = theta_cc_mul_nd;
    for k = 1:length(sol4c.x)
        C__1 = sol4c.y(1,k);
        C__2 = sol4c.y(2,k);
        N_poss_cem_1_0(k) = -(d__pos * C__2) + d__pos * (-omega__c * A + ((-1 + d__pos) * C__2) + omega__c * (A * (xi - 0.5e0) + B) * Pe__val + i__val) * C__1 / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) + C__1 * Pe__val;
        N_negs_cem_1_0(k) = (-d__pos * ((A * (xi - 0.5e0) + B) * omega__c + (2 * C__1)) * C__2 - A * C__1 * omega__c * d__pos + (Pe__val * (0.1e1 + d__pos) * C__1 - i__val) * (C__1 + (A * (xi - 0.5e0) + B) * omega__c)) / ((A * (xi - 0.5e0) + B) * omega__c + C__1 * (0.1e1 + d__pos));
    end
N_poss_cem_1(j) = mean(N_poss_cem_1_0);
N_negs_cem_1(j) = mean(N_negs_cem_1_0);
end
plot(i_act,N_poss_cem_1,'-o','LineWidth',3,'Color',"red");hold on;
plot(i_act,N_negs_cem_1,'-o','LineWidth',3,'Color',"yellow");
     title("Co-ion and Counter ion flux as func of i")
    legend("N_{pos} for CEM","N_{neg} for CEM")
        xlabel('i(current density)')
        ylabel('N_{i}^{*}')
    hold off;
    saveas(h(15),"plots_sv/f7_pos_i_var_flux_aem.png")
    fprintf(2,"Plot: "+ num2str(15)+" Done\n");
    %% For  Calculating positive Flux in AEM for varying Pe for Diff diff ratios
    pec_vals = [-0.1,-1];
for pec  = 1:length(pec_vals)
    theta_cc_mul=linspace(0.01,4,80);
    theta_cc_mul_nd = theta_cc_mul/td__c;
    xmesh = linspace(0, 1, 100);
    Pe__val = pec_vals(pec);
    i__val = 0;
    d__pos = 0.5;
    h(16+pec-1)= figure;
    omega__c = 1;
    %Claculating BCs
    C__posL2 = 1;
    C__posR2 = 10;
    N_poss_D1 = [];
    bcs_0 = [];
    bcs_1 = [];
    for i = 1:length(theta_cc_mul)
        bcs_0(i) = -theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posL2^2)/2;
        bcs_1(i) = -theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posR2^2)/2;
    end
    for i=1:length(theta_cc_mul)
        opts = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",d__pos,0,Pe__val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
        solinit = bvpinit(xmesh,[0; 1]);
        sol4c = bvp4c(bvp_call_mod_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",d__pos,0,Pe__val), bcf_call(bcs_0(i),bcs_1(i)),solinit, opts);
        N_poss_0 = zeros(length(sol4c.x),1);
        for k = 1:length(sol4c.x)
            A = 0;
            B = theta_cc_mul_nd(i);
            C__1 = sol4c.y(1,k);
            C__2 = sol4c.y(2,k);
            xi = sol4c.x(1,k);
            N_poss_0(k) = -(d__pos * C__2) + d__pos * (-omega__c * A + ((-1 + d__pos) * C__2) + omega__c * (A * (xi - 0.5e0) + B) * Pe__val + i__val) * C__1 / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) + C__1 * Pe__val;
        end
        N_poss_D1(i) = mean(N_poss_0);
    end
    omega__c = -1;
    %Claculating BCs
    C__posL2 = 1;
    C__posR2 = 10;
    N_negs_D1 = [];
    bcs_0 = [];
    bcs_1 = [];
    for i = 1:length(theta_cc_mul)
        bcs_0(i) = theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posL2^2)/2;
        bcs_1(i) = theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posR2^2)/2;
    end
    for i=1:length(theta_cc_mul)
        opts = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",d__pos,0,Pe__val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
        solinit = bvpinit(xmesh,[300; 0.5]);
        sol4c = bvp4c(bvp_call_mod_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",d__pos,0,Pe__val), bcf_call(bcs_0(i),bcs_1(i)),solinit, opts);
        N_negs_0 = zeros(length(sol4c.x),1);
        for k = 1:length(sol4c.x)
            A = 0;
            B = theta_cc_mul_nd(i);
            C__1 = sol4c.y(1,k);
            C__2 = sol4c.y(2,k);
            xi = sol4c.x(1,k);
            N_negs_0(k) = (-d__pos * ((A * (xi - 0.5e0) + B) * omega__c + (2 * C__1)) * C__2 - A * C__1 * omega__c * d__pos + (Pe__val * (0.1e1 + d__pos) * C__1 - i__val) * (C__1 + (A * (xi - 0.5e0) + B) * omega__c)) / ((A * (xi - 0.5e0) + B) * omega__c + C__1 * (0.1e1 + d__pos));
        end
        N_negs_D1(i) = mean(N_negs_0);
    end
    % For Diff ratio 0.1
    d__pos = 0.1;
    omega__c = 1;
    %Claculating BCs
    C__posL2 = 1;
    C__posR2 = 10;
    N_poss_D2 = [];
    bcs_0 = [];
    bcs_1 = [];
    for i = 1:length(theta_cc_mul)
        bcs_0(i) = -theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posL2^2)/2;
        bcs_1(i) = -theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posR2^2)/2;
    end
    for i=1:length(theta_cc_mul)
        opts = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",d__pos,0,Pe__val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
        solinit = bvpinit(xmesh,[0; 1]);
        sol4c = bvp4c(bvp_call_mod_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",d__pos,0,Pe__val), bcf_call(bcs_0(i),bcs_1(i)),solinit, opts);
        N_poss_0 = zeros(length(sol4c.x),1);
        for k = 1:length(sol4c.x)
            A = 0;
            B = theta_cc_mul_nd(i);
            C__1 = sol4c.y(1,k);
            C__2 = sol4c.y(2,k);
            N_poss_0(k) = -(d__pos * C__2) + d__pos * (-omega__c * A + ((-1 + d__pos) * C__2) + omega__c * (A * (xi - 0.5e0) + B) * Pe__val + i__val) * C__1 / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) + C__1 * Pe__val;
        end
        N_poss_D2(i) = mean(N_poss_0);
    end
    omega__c = -1;
    %Claculating BCs
    C__posL2 = 1;
    C__posR2 = 10;
    N_negs_D2 = [];
    bcs_0 = [];
    bcs_1 = [];
    for i = 1:length(theta_cc_mul)
        bcs_0(i) = theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posL2^2)/2;
        bcs_1(i) = theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posR2^2)/2;
    end
    for i=1:length(theta_cc_mul)
        opts = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",d__pos,0,Pe__val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
        solinit = bvpinit(xmesh,[300; 0.5]);
        sol4c = bvp4c(bvp_call_mod_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",d__pos,0,Pe__val), bcf_call(bcs_0(i),bcs_1(i)),solinit, opts);
        N_negs_0 = zeros(length(sol4c.x),1);
        for k = 1:length(sol4c.x)
            A = 0;
            B = theta_cc_mul_nd(i);
            C__1 = sol4c.y(1,k);
            C__2 = sol4c.y(2,k);
            xi = sol4c.x(1,k);
            N_negs_0(k) = (-d__pos * ((A * (xi - 0.5e0) + B) * omega__c + (2 * C__1)) * C__2 - A * C__1 * omega__c * d__pos + (Pe__val * (0.1e1 + d__pos) * C__1 - i__val) * (C__1 + (A * (xi - 0.5e0) + B) * omega__c)) / ((A * (xi - 0.5e0) + B) * omega__c + C__1 * (0.1e1 + d__pos));
        end
        N_negs_D2(i) = mean(N_negs_0);
    end
    shp = "-o";
    plot(theta_cc_mul,N_poss_D1,shp,'LineWidth',3);hold on;
    plot(theta_cc_mul,N_negs_D1,shp,'LineWidth',3);hold on;
    plot(theta_cc_mul,N_poss_D2,shp,'LineWidth',3);hold on;
    plot(theta_cc_mul,N_negs_D2,shp,'LineWidth',3);hold on;
    % Combined postive and negative flux
    title("N^{*} vs \theta for Pe ="+num2str(Pe__val))
        xlabel('\theta')
        ylabel('N^{*}(\theta)')
        legend("N_{pos} for AEM D1","N_{neg} for CEM D1","N_{pos} for AEM D2","N_{neg} for CEM D2")
        hold off;
        saveas(h(16+pec-1),"plots_sv/s3_Nneg_Npos_combined_Pe_Var"+num2str(Pe__val)+".png")
        fprintf(2,"Plot: "+ num2str(16+pec-1)+" Done\n");
end
    %% For  Calculating positive Flux in AEM for varying Pe for i change ratios
    shp = '-o';
    pec_vals = [-1,-10];
for pec  = 1:length(pec_vals)
    theta_cc_mul=linspace(0.01,4,80);
    theta_cc_mul_nd = theta_cc_mul/td__c;
    xi = linspace(0,1,100);
    xmesh = linspace(0, 1, 100);
    Pe__Val = pec_vals(pec);
    i__Val = 1000;
    d__pos = 0.5;
    h(20+pec-1)=figure;
    omega__c = 1;
    %Claculating BCs
    C__posL2 = 1;
    C__posR2 = 10;
    N_poss_D1 = [];
    bcs_0 = [];
    bcs_1 = [];
    for i = 1:length(theta_cc_mul)
        bcs_0(i) = -theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posL2^2)/2;
        bcs_1(i) = -theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posR2^2)/2;
    end
    for i=1:length(theta_cc_mul)
        opts = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",d__pos,0,Pe__Val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
        solinit = bvpinit(xmesh,[0; 1]);
        sol4c = bvp4c(bvp_call_mod_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",d__pos,0,Pe__Val), bcf_call(bcs_0(i),bcs_1(i)),solinit, opts);
        N_poss_0 = zeros(length(sol4c.x),1);
        for k = 1:length(sol4c.x)
            A__1 = 0;
            B__1 = theta_cc_mul_nd(i);
            C__1 = sol4c.y(1,k);
            C__2 = sol4c.y(2,k);
            N_poss_0(k) = -C__2*d__pos - d__pos*C__1*((1 - d__pos)*C__2 - omega__c*(A__1*(xi - 0.5) + B__1)*Pe__Val + omega__c*A__1 - i__Val)/((A__1*(xi - 0.5) + B__1)*omega__c + C__1*(1 + d__pos)) + C__1*Pe__Val; 
        end
        N_poss_D1(i) = mean(N_poss_0);
    end
    % For Calculating negative Flux
    omega__c = -1;
    %Claculating BCs
    N_negs_D1 = [];
    bcs_0 = [];
    bcs_1 = [];
    for i = 1:length(theta_cc_mul)
        bcs_0(i) = theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posL2^2)/2;
        bcs_1(i) = theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posR2^2)/2;
    end
    for i=1:length(theta_cc_mul)
        opts = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",d__pos,0,Pe__Val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
        solinit = bvpinit(xmesh,[300; 0.5]);
        sol4c = bvp4c(bvp_call_mod_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",d__pos,0,Pe__Val), bcf_call(bcs_0(i),bcs_1(i)),solinit, opts);
        N_negs_0 = zeros(length(sol4c.x),1);
        for k = 1:length(sol4c.x)
            A__1 = 0;
            B__1 = theta_cc_mul_nd(i);
            C__1 = sol4c.y(1,k);
            C__2 = sol4c.y(2,k);
            N_negs_0(k) = (((((-0.1e1 + (xi - 0.5e0) * Pe__Val) * A__1 + Pe__Val * B__1) * C__1 + (-0.1e1 * xi + 0.5e0) * C__2 * A__1 - 0.1e1 * C__2 * B__1) * d__pos + ((xi - 0.5e0) * Pe__Val * A__1 + Pe__Val * B__1) * C__1) * omega__c + (0.2e1 * C__2 * C__1 - 0.1e1 * C__1 ^ 2 * Pe__Val) * d__pos + C__1 * i__Val - 0.1e1 * C__1 ^ 2 * Pe__Val) / ((A__1 * (xi - 0.5e0) + B__1) * omega__c * d__pos - 0.1e1 * C__1 * d__pos - 0.1e1 * C__1);
        end
        N_negs_D1(i) = mean(N_negs_0);
    end
    plot(theta_cc_mul,N_poss_D1,shp,'LineWidth',3);hold on;
    plot(theta_cc_mul,N_negs_D1,shp,'LineWidth',3);hold on;
    % Combined postive and negative flux
    title("N^{*} vs \theta for Pe ="+num2str(Pe__Val))
        xlabel('\theta')
        ylabel('N^{*}(\theta)')
        legend("N_{pos} for AEM D1","N_{neg} for CEM D1")
        hold off;
        saveas(h(20+pec-1),"plots_sv/s6_Nneg_Npos_combined_Pe_Var_dpos_ich"+num2str(Pe__Val)+".png")
         fprintf(2,"Plot: "+ num2str(20+pec-1)+" Done\n");
end
    %% For  Calculating positive Flux in AEM for varying Pe for i change ranging 0:1000
    pec_vals = [0,-1,-10,-100];
for pec  = 1:length(pec_vals)
    theta_cc_mul=2*ones(80,1);
    theta_cc_mul_nd = theta_cc_mul/td__c;
    xi = linspace(0,1,100);
    xmesh = linspace(0, 1, 100);
    Pe__Val = pec_vals(pec);
    i__Val_rng = linspace(0,1000,length(theta_cc_mul));
    d__pos = 0.5;
    h(22+pec-1)=figure;
    omega__c = 1;
    %Claculating BCs
    C__posL2 = 1;
    C__posR2 = 10;
    N_poss_D1 = [];
    bcs_0 = [];
    bcs_1 = [];
    for i = 1:length(theta_cc_mul)
        bcs_0(i) = -theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posL2^2)/2;
        bcs_1(i) = -theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posR2^2)/2;
    end
    for i=1:length(theta_cc_mul)
        i__Val = i__Val_rng(i);
        opts = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",d__pos,i__Val,Pe__Val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
        solinit = bvpinit(xmesh,[0; 1]);
        sol4c = bvp4c(bvp_call_mod_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",d__pos,i__Val,Pe__Val), bcf_call(bcs_0(i),bcs_1(i)),solinit, opts);
        N_poss_0 = zeros(length(sol4c.x),1);
        for k = 1:length(sol4c.x)
            A__1 = 0;
            B__1 = theta_cc_mul_nd(i);
            C__1 = sol4c.y(1,k);
            C__2 = sol4c.y(2,k);
            N_poss_0(k) = -C__2*d__pos - d__pos*C__1*((1 - d__pos)*C__2 - omega__c*(A__1*(xi - 0.5) + B__1)*Pe__Val + omega__c*A__1 - i__Val)/((A__1*(xi - 0.5) + B__1)*omega__c + C__1*(1 + d__pos)) + C__1*Pe__Val; 
        end
        N_poss_D1(i) = mean(N_poss_0);
    end
    % For Calculating negative Flux
    omega__c = -1;
    %Claculating BCs
    C__negL2 = 1;
    C__negR2 = 10;
    N_negs_D1 = [];
    bcs_0 = [];
    bcs_1 = [];
    for i = 1:length(theta_cc_mul)
        bcs_0(i) = theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posL2^2)/2;
        bcs_1(i) = theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posR2^2)/2;
    end
    for i=1:length(theta_cc_mul)
        opts = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",d__pos,0,Pe__Val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
        solinit = bvpinit(xmesh,[300; 0.5]);
        sol4c = bvp4c(bvp_call_mod_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",d__pos,0,Pe__Val), bcf_call(bcs_0(i),bcs_1(i)),solinit, opts);
        N_negs_0 = zeros(length(sol4c.x),1);
        for k = 1:length(sol4c.x)
            A__1 = 0;
            B__1 = theta_cc_mul_nd(i);
            C__1 = sol4c.y(1,k);
            C__2 = sol4c.y(2,k);
            N_negs_0(k) = (((((-0.1e1 + (xi - 0.5e0) * Pe__Val) * A__1 + Pe__Val * B__1) * C__1 + (-0.1e1 * xi + 0.5e0) * C__2 * A__1 - 0.1e1 * C__2 * B__1) * d__pos + ((xi - 0.5e0) * Pe__Val * A__1 + Pe__Val * B__1) * C__1) * omega__c + (0.2e1 * C__2 * C__1 - 0.1e1 * C__1 ^ 2 * Pe__Val) * d__pos + C__1 * i__Val - 0.1e1 * C__1 ^ 2 * Pe__Val) / ((A__1 * (xi - 0.5e0) + B__1) * omega__c * d__pos - 0.1e1 * C__1 * d__pos - 0.1e1 * C__1);
        end
        N_negs_D1(i) = mean(N_negs_0);
    end
    plot(i__Val_rng,N_poss_D1,shp,'LineWidth',3);hold on;
    plot(i__Val_rng,N_negs_D1,shp,'LineWidth',3);hold on;
    % Combined postive and negative flux
    title("N^{*} vs \theta for Pe ="+num2str(Pe__Val))
        xlabel('i')
        ylabel('N^{*}(\theta)')
        legend("N_{pos} for AEM D1","N_{neg} for CEM D1")
        hold off;
        saveas(h(22+pec-1),"plots_sv/s7_Nneg_Npos_combined_Pe_Var_irng"+num2str(Pe__Val)+".png")
        fprintf(2,"Plot: "+ num2str(22+pec-1)+" Done\n");

end
%% For  Calculating positive Flux in AEM for varying Pe for inhomo change on one side
pec_vals = [-1,-10];
for pec  = 1:length(pec_vals)
    xi = linspace(0,1,100);
    xmesh = linspace(0.01, 1, 100);
    A = 400;
    B = 200;
    theta_cc_mul = A*(xi-0.5)+B;
    theta_cc_mul_nd = theta_cc_mul;
    Pe__val = pec_vals(pec);
    i__val = 1000;
    d__pos = 0.5;
    h(26+pec-1)=figure;
    omega__c = 1;
    if pec_vals(pec) == -10
              i__val = -1000;
    end
    %Claculating BCs
    C__posL2 = 1;
    C__posR2 = 10;
    N_poss_D1 = [];
    bcs_0 = [];
    bcs_1 = [];
    for i = 1:length(theta_cc_mul)
        bcs_0(i) = -theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posL2^2)/2;
        bcs_1(i) = -theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posR2^2)/2;
    end
    for i=1:length(theta_cc_mul)
        opts = bvpset('FJacobian',jac_call_theta_pe_var(A,B,omega__c,"pos",0.5,i__val,Pe__val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
        solinit = bvpinit(xmesh,[0.005; -1]);
        sol4c = bvp4c(bvp_call_mod_theta_pe_var(A,B,omega__c,"pos",0.5,i__val,Pe__val), bcf_call(bcs_0(i),bcs_1(i)),solinit, opts);
        N_poss_0 = zeros(length(sol4c.x),1);
        for k = 1:length(sol4c.x)
            C__1 = sol4c.y(1,k);
            C__2 = sol4c.y(2,k);
            xi = sol4c.x(1,k);
            N_poss_0(k) = -(d__pos * C__2) + d__pos * (-omega__c * A + ((-1 + d__pos) * C__2) + omega__c * (A * (xi - 0.5e0) + B) * Pe__val + i__val) * C__1 / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) + C__1 * Pe__val;
        end
        N_poss_D1(i) = mean(N_poss_0);
    end
    % For Calculating negative Flux
    omega__c = -1;
    %Claculating BCs
    C__negL2 = 1;
    C__negR2 = 10;
    if pec_vals(pec) == -10
        i__Val = -1000;
    end
    N_negs_D1 = [];
    bcs_0 = [];
    bcs_1 = [];
    for i = 1:length(theta_cc_mul)
        bcs_0(i) = theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posL2^2)/2;
        bcs_1(i) = theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posR2^2)/2;
    end
    for i=1:length(theta_cc_mul)
        opts = bvpset('FJacobian',jac_call_theta_pe_var(A,B,omega__c,"pos",0.5,i__val,Pe__val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
        solinit = bvpinit(xmesh,[3000; -0.5 ...
            ]);
        sol4c = bvp4c(bvp_call_mod_theta_pe_var(A,B,omega__c,"pos",0.5,i__val,Pe__val), bcf_call(bcs_0(i),bcs_1(i)),solinit, opts);
        N_negs_0 = zeros(length(sol4c.x),1);
        for k = 1:length(sol4c.x)
            C__1 = sol4c.y(1,k);
            C__2 = sol4c.y(2,k);
            xi = sol4c.x(1,k);
            N_negs_0(k) = (-d__pos * ((A * (xi - 0.5e0) + B) * omega__c + (2 * C__1)) * C__2 - A * C__1 * omega__c * d__pos + (Pe__val * (0.1e1 + d__pos) * C__1 - i__val) * (C__1 + (A * (xi - 0.5e0) + B) * omega__c)) / ((A * (xi - 0.5e0) + B) * omega__c + C__1 * (0.1e1 + d__pos));
        end
        N_negs_D1(i) = mean(N_negs_0);
    end
    shp = '-o';
    plot(theta_cc_mul_nd,N_poss_D1,shp,'LineWidth',3);hold on;
    plot(theta_cc_mul_nd,N_negs_D1,shp,'LineWidth',3);hold on;
    % Combined postive and negative flux
    title("N^{*} vs \theta for Pe ="+num2str(Pe__val))
        xlabel('i')
        ylabel('N^{*}(\theta)')
        legend("N_{pos} for AEM D1","N_{neg} for CEM D1")
        hold off;
        saveas(h(26+pec-1),"plots_sv/s8_Nneg_Npos_combined_Pe_Var_inh_oneside"+num2str(Pe__val)+".png")
        fprintf(2,"Plot: "+ num2str(26+pec-1)+" Done\n");
end
%% For  Calculating positive Flux in AEM for varying Pe for inhomo change on both side
pec_vals = [-1,-10];
for pec  = 1:length(pec_vals)
    xmesh = linspace(0.0001, 1, 80);
    xi = xmesh;
    A = 400;
    B = 0;
    theta_cc_mul = [];
    for i = 1:length(xmesh)
        if xi(i)<=0.5
        theta_cc_mul(i)= -A*(xi(i)-0.5);
        else
        theta_cc_mul(i) = A*(xi(i)-0.5);
        end
    end
    theta_cc_mul_nd = theta_cc_mul;
    theta_cc_mul_dim = theta_cc_mul*td__c;
    Pe__val = pec_vals(pec);
    i__val = 0;
    d__pos = 0.5;
    h(28+pec-1)=figure;
    omega__c = 1;
    if pec_vals(pec) == -10
              i__val = 1000;
    end
    %Claculating BCs
    C__posL2 = 1;
    C__posR2 = 10;
    N_poss_D1 = [];
    bcs_0 = [];
    bcs_1 = [];
    for i = 1:length(theta_cc_mul)
        bcs_0(i) = -theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posL2^2)/2;
        bcs_1(i) = -theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posR2^2)/2;
    end
    for i=1:length(theta_cc_mul)
        opts = bvpset('FJacobian',jac_call_theta_pe_var(A,B,omega__c,"pos_inhm",0.5,i__Val,Pe__val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
        solinit = bvpinit(xmesh,[0.005; -1]);
        sol4c = bvp4c(bvp_call_mod_theta_pe_var(A,B,omega__c,"pos_inhm",0.5,i__Val,Pe__val), bcf_call(bcs_0(i),bcs_1(i)),solinit, opts);
        %[sol4c.y(1,1) sol4c.y(2,1)]
        N_poss_0 = zeros(length(sol4c.x),1);
        for k = 1:length(sol4c.x)
            C__1 = sol4c.y(1,k);
            C__2 = sol4c.y(2,k);
            xi = sol4c.x(1,k);
            if xi <=0.5
                A = A;
            N_poss_0(k) = -(d__pos * C__2) + d__pos * (-omega__c * A + ((-1 + d__pos) * C__2) + omega__c * (A * (xi - 0.5e0) + B) * Pe__val + i__val) * C__1 / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) + C__1 * Pe__val;
            else
                A = -A;
            N_poss_0(k) = -(d__pos * C__2) + d__pos * (-omega__c * A + ((-1 + d__pos) * C__2) + omega__c * (A * (xi - 0.5e0) + B) * Pe__val + i__val) * C__1 / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) + C__1 * Pe__val;                
            end
        end
        N_poss_D1(i) = mean(N_poss_0);
    end
    % For Calculating negative Flux
    omega__c = -1;
    %Claculating BCs
    C__negL2 = 1;
    C__negR2 = 10;
    if pec_vals(pec) == -10
        i__val = 1000;
    end
    N_negs_D1 = [];
    bcs_0 = [];
    bcs_1 = [];
    for i = 1:length(theta_cc_mul)
        bcs_0(i) = theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posL2^2)/2;
        bcs_1(i) = theta_cc_mul_nd(i)/2 + sqrt(theta_cc_mul_nd(i)^2 + 4*C__posR2^2)/2;
    end
    for i=1:length(theta_cc_mul)
        opts = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",d__pos,0,Pe__val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
        solinit = bvpinit(xmesh,[300; 0.5]);
        sol4c = bvp4c(bvp_call_mod_theta_pe_var(0,theta_cc_mul_nd(i),omega__c,"pos",d__pos,0,Pe__val), bcf_call(bcs_0(i),bcs_1(i)),solinit, opts);
        N_negs_0 = zeros(length(sol4c.x),1);
        for k = 1:length(sol4c.x)
            C__1 = sol4c.y(1,k);
            C__2 = sol4c.y(2,k);
            xi = sol4c.x(1,k);
             if xi <=0.5
                A = A;
            N_negs_0(k) = (-d__pos * ((A * (xi - 0.5e0) + B) * omega__c + (2 * C__1)) * C__2 - A * C__1 * omega__c * d__pos + (Pe__val * (0.1e1 + d__pos) * C__1 - i__val) * (C__1 + (A * (xi - 0.5e0) + B) * omega__c)) / ((A * (xi - 0.5e0) + B) * omega__c + C__1 * (0.1e1 + d__pos));
            else
                A = -A;
            N_negs_0(k) = (-d__pos * ((A * (xi - 0.5e0) + B) * omega__c + (2 * C__1)) * C__2 - A * C__1 * omega__c * d__pos + (Pe__val * (0.1e1 + d__pos) * C__1 - i__val) * (C__1 + (A * (xi - 0.5e0) + B) * omega__c)) / ((A * (xi - 0.5e0) + B) * omega__c + C__1 * (0.1e1 + d__pos));
             end
        end
        N_negs_D1(i) = mean(N_negs_0);
    end
    shp = '-o';
    plot(theta_cc_mul_nd,N_poss_D1,shp,'LineWidth',3);hold on;
    plot(theta_cc_mul_nd,N_negs_D1,shp,'LineWidth',3);hold on;
    % Combined postive and negative flux
    title("N^{*} vs \theta for Pe ="+num2str(Pe__val))
        xlabel('\theta')
        ylabel('N^{*}(\theta)')
        legend("N_{pos} for AEM D1","N_{neg} for CEM D1")
        hold off;
        saveas(h(28+pec-1),"plots_sv/s9_Nneg_Npos_combined_Pe_Var_inh_2side"+num2str(Pe__val)+".png")
        fprintf(2,"Plot: "+ num2str(28+pec-1)+" Done\n");
end
%% For  Calculating positive Flux in AEM for inhomogeneous both side
xmesh = linspace(0.01, 1, 200);
td__c = 0.01;
xi = xmesh;
h(39) = figure;
omega__c = 1;
A_1_rng = linspace(5,400,100);
N_poss_0_t = [];
i__val = 0;
Pe_val = 0;
d__pos = 0.5;
for t = 1:length(A_1_rng)
    A_1 = A_1_rng(t);
    B_1 = 0;
    theta_cc_mul = [];
        for i = 1:length(xmesh)
            if xi(i)<=0.5
            theta_cc_mul(i)= -A_1*(xi(i)-0.5);
            else
            theta_cc_mul(i) = A_1*(xi(i)-0.5);
            end
        end
    theta_cc_mul_nd = theta_cc_mul;
    %Claculating BCs
    C__posL2 = 1;
    C__posR2 = 10;
        bcs_0 = -theta_cc_mul_nd(1)/2 + sqrt(theta_cc_mul_nd(1)^2 + 4*C__posL2^2)/2;
        bcs_1 = -theta_cc_mul_nd(end)/2 + sqrt(theta_cc_mul_nd(end)^2 + 4*C__posR2^2)/2;
        opts = bvpset('FJacobian',jac_call_theta_pe_var(A_1,B_1,omega__c,"pos_inhm",0.5,i__val,0),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
        solinit = bvpinit(xmesh,[0; 1]);
        sol4c = bvp4c(bvp_call_mod_theta_pe_var(A_1,B_1,omega__c,"pos_inhm",0.5,i__val,0), bcf_call(bcs_0,bcs_1),solinit, opts);
        N_poss_0 = zeros(length(sol4c.x),1);
        theta_cc_mul_dim_plt = zeros(length(sol4c.x),1);
        for k = 1:length(sol4c.x)
            C__1 = sol4c.y(1,k);
            C__2 = sol4c.y(2,k);
            xi_vl = sol4c.x(1,k);
            if xi_vl<=0.5
                A = -A_1;
               N_poss_0(k) = -(d__pos * C__2) + d__pos * (-omega__c * A + ((-1 + d__pos) * C__2) + omega__c * (A * (xi - 0.5e0) + B) * Pe__val + i__val) * C__1 / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) + C__1 * Pe__val;
            else
                A = A_1;
              N_poss_0(k) = -(d__pos * C__2) + d__pos * (-omega__c * A + ((-1 + d__pos) * C__2) + omega__c * (A * (xi - 0.5e0) + B) * Pe__val + i__val) * C__1 / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) + C__1 * Pe__val;
            end
        end
        N_poss_0_t(t) = mean(N_poss_0);
end
plot(A_1_rng*td__c,N_poss_0_t,'-o','LineWidth',3)
title("N_{pos}^{*} vs \theta for AEM")
    xlabel('\theta')
    ylabel('N_{pos}^{*}(\theta)')
    hold off;
    saveas(h(39),"f11_Npos_aem_theta_var_inhm_2sid.png")
%% For Calculating negative Flux
xmesh = linspace(0.01, 1, 200);
td__c = 0.01;
xi_m = xmesh;
h(40) = figure;
omega__c = -1;
i__val = 0;
Pe__val = 0;
%Claculating BCs
N_negs_0_t = [];
for t = 1:length(A_1_rng)
        A_1 = A_1_rng(t);
        B_1 = 0;
        theta_cc_mul = [];
        for i = 1:length(xmesh)
            if xi_m(i)<=0.5
            theta_cc_mul(i)= -A_1*(xi_m(i)-0.5);
            else
            theta_cc_mul(i) = A_1*(xi_m(i)-0.5);
            end
        end
    C__posL2 = 1;
    C__posR2 = 10;
        bcs_0 = theta_cc_mul_nd(1)/2 + sqrt(theta_cc_mul_nd(1)^2 + 4*C__posL2^2)/2;
        bcs_1 = theta_cc_mul_nd(end)/2 + sqrt(theta_cc_mul_nd(end)^2 + 4*C__posR2^2)/2;
        opts = bvpset('FJacobian',jac_call_theta_pe_var(A_1,B_1,omega__c,"pos_inhm",0.5,i__val,0),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
        solinit = bvpinit(xmesh,[0; 0]);
        sol4c = bvp4c(bvp_call_mod_theta_pe_var(A_1,B_1,omega__c,"pos_inhm",0.5,i__val,0), bcf_call(bcs_0,bcs_1),solinit, opts);
        N_negs_0 = zeros(length(sol4c.x),1);
        theta_cc_mul_dim_plt_N = zeros(length(sol4c.x),1);
        for k = 1:length(sol4c.x)
            C__1 = sol4c.y(1,k);
            C__2 = sol4c.y(2,k);
            xi_vl = sol4c.x(1,k);
            if xi_vl <=0.5
                A = -A_1;
                B = B_1;
                N_negs_0(k) = (((((-0.1e1 + (xi - 0.5e0) * Pe__Val) * A__1 + Pe__Val * B__1) * C__1 + (-0.1e1 * xi + 0.5e0) * C__2 * A__1 - 0.1e1 * C__2 * B__1) * d__pos + ((xi - 0.5e0) * Pe__Val * A__1 + Pe__Val * B__1) * C__1) * omega__c + (0.2e1 * C__2 * C__1 - 0.1e1 * C__1 ^ 2 * Pe__Val) * d__pos + C__1 * i__Val - 0.1e1 * C__1 ^ 2 * Pe__Val) / ((A__1 * (xi - 0.5e0) + B__1) * omega__c * d__pos - 0.1e1 * C__1 * d__pos - 0.1e1 * C__1);
            else
                A = A_1;
                B = B_1;
                N_negs_0(k) = (((((-0.1e1 + (xi - 0.5e0) * Pe__Val) * A__1 + Pe__Val * B__1) * C__1 + (-0.1e1 * xi + 0.5e0) * C__2 * A__1 - 0.1e1 * C__2 * B__1) * d__pos + ((xi - 0.5e0) * Pe__Val * A__1 + Pe__Val * B__1) * C__1) * omega__c + (0.2e1 * C__2 * C__1 - 0.1e1 * C__1 ^ 2 * Pe__Val) * d__pos + C__1 * i__Val - 0.1e1 * C__1 ^ 2 * Pe__Val) / ((A__1 * (xi - 0.5e0) + B__1) * omega__c * d__pos - 0.1e1 * C__1 * d__pos - 0.1e1 * C__1);
            end            
        end
       N_negs_0_t(t) = mean(N_negs_0);
end
plot(A_1_rng*td__c,N_negs_0_t,'-o','LineWidth',3)
title("N_{neg}^{*} vs \theta for CEM")
    xlabel('\theta')
    ylabel('N_{neg}^{*}(\theta)')
    hold off;
    saveas(h(40),"f11_Nneg_cem_theta_var_inhm_2sid.png")
%% Combined postive and negative flux
h(41) = figure;
plot(A_1_rng*td__c,N_poss_0_t,'-o','LineWidth',3);hold on;
plot(A_1_rng*td__c,N_negs_0_t,'-o','LineWidth',3)
title("N^{*} vs \theta")
    xlabel('\theta')
    ylabel('N^{*}(\theta)')
    legend("N_{pos} for AEM","N_{neg} for CEM")
    hold off;
    saveas(h(41),"plots_sv/f11_Nneg_Npos_combined_var_inhm_2sid.png")
%% For  Calculating positive Flux in AEM for inhomogeneous one side i var
i__Vals = [-500,500];
for irng = 1:length(i__Vals)
    irng_0 = 0;
    if irng ==2
    irng_0 = 2;
    end
    i__val = i__Vals(irng);
    xmesh = linspace(0.0001, 1,212);
        td__c = 0.01;
        xi_m = xmesh;
        diff_cnst = 0.5;
        Pe__val = 0;
        h(33+irng+irng_0-1) = figure;
        omega__c = 1;
        N_poss_0_t = [];
        A_rng = linspace(0.002,200,100);
        B_rng = linspace(0.001,100,100);
        for t = 1:length(A_rng)
            A = A_rng(t);
            B = B_rng(t);
            theta_cc_mul = A*(xi_m-0.5)+B;
            theta_cc_mul_nd = theta_cc_mul;
            %Claculating BCs
            C__posL2 = 1;
            C__posR2 = 10;
                bcs_0 = -theta_cc_mul_nd(1)/2 + sqrt(theta_cc_mul_nd(1)^2 + 4*C__posL2^2)/2;
                bcs_1 = -theta_cc_mul_nd(end)/2 + sqrt(theta_cc_mul_nd(end)^2 + 4*C__posR2^2)/2;
                opts = bvpset('FJacobian',jac_call_theta_pe_var(A,B,omega__c,"pos",diff_cnst,i__val,Pe__val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
                solinit = bvpinit(xmesh,[1.5; -0.5]);
                sol4c = bvp4c(bvp_call_mod_theta_pe_var(A,B,omega__c,"pos",diff_cnst,i__val,Pe__val), bcf_call(bcs_0,bcs_1),solinit, opts);
                N_poss_0 = zeros(length(sol4c.x),1);
                theta_cc_mul_dim_plt = zeros(length(sol4c.x),1);
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
            plot(A_rng*td__c,N_poss_0_t,'-o','LineWidth',3)
            title("N_{pos}^{*} vs \theta for AEM")
            xlabel('\theta')
            ylabel('N_{pos}^{*}(\theta)')
            hold off;
            saveas(h(33+irng+irng_0-1),"f10_Npos_aem_theta_var_inhm.png")
        % For Calculating negative Flux in CEM
        h(34+irng+irng_0-1) = figure;
        xmesh = linspace(0.0001, 1,212);
        td__c = 0.01;
        xi_m = xmesh;
        diff_cnst = 0.5;
        Pe__val = 0;
        omega__c = -1;
        N_negs_0_t = [];
        for t = 1:length(A_rng)
            A = A_rng(t);
            B = B_rng(t);
            theta_cc_mul = A*(xi_m-0.5)+B;
            theta_cc_mul_nd = theta_cc_mul;
            %Claculating BCs
            C__posL2 = 1;
            C__posR2 = 10;
                bcs_0 = theta_cc_mul_nd(1)/2 + sqrt(theta_cc_mul_nd(1)^2 + 4*C__posL2^2)/2;
                bcs_1 = theta_cc_mul_nd(end)/2 + sqrt(theta_cc_mul_nd(end)^2 + 4*C__posR2^2)/2;
                opts = bvpset('FJacobian',jac_call_theta_pe_var(A,B,omega__c,"pos",diff_cnst,i__val,Pe__val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
                kol = -0.5;
%                 if i__val ==500
%                     kol = 1;
%                 end
                solinit = bvpinit(xmesh,[300;-1]);
                sol4c = bvp4c(bvp_call_mod_theta_pe_var(A,B,omega__c,"pos",diff_cnst,i__val,Pe__val), bcf_call(bcs_0,bcs_1),solinit, opts);
                [sol4c.y(1,1), sol4c.y(2,1)]
                N_poss_0 = zeros(length(sol4c.x),1);
                theta_cc_mul_dim_plt = zeros(length(sol4c.x),1);
                d__pos = diff_cnst;
                for k = 1:length(sol4c.x)
                    C__1 = sol4c.y(1,k);
                    C__2 = sol4c.y(2,k);
                    xi_vl = sol4c.x(1,k);
                    xi = xi_vl;
                    N_negs_0(k) = (-d__pos * ((A * (xi - 0.5e0) + B) * omega__c + (2 * C__1)) * C__2 - A * C__1 * omega__c * d__pos + (Pe__val * (0.1e1 + d__pos) * C__1 - i__val) * (C__1 + (A * (xi - 0.5e0) + B) * omega__c)) / ((A * (xi - 0.5e0) + B) * omega__c + C__1 * (0.1e1 + d__pos));               
                end
                N_negs_0_t(t) = mean(N_negs_0);
        end
        plot(A_rng*td__c,N_negs_0_t,'-o','LineWidth',3)
        title("N_{neg}^{*} vs \theta for CEM")
            xlabel('\theta')
            ylabel('N_{neg}^{*}(\theta)')
            hold off;
            saveas(h(34+irng+irng_0-1),"f10_Nneg_cem_theta_var_inhm.png")
        % Combined postive and negative flux
        h(35+irng+irng_0-1) = figure;
        plot(A_rng*td__c,N_poss_0_t,'-o','LineWidth',3);hold on;
        plot(A_rng*td__c,N_negs_0_t,'-o','LineWidth',3)
        title("N^{*} vs \theta")
            xlabel('\theta')
            ylabel('N^{*}(\theta)')
            legend("N_{pos} for AEM","N_{neg} for CEM")
            hold off;
            saveas(h(35+irng+irng_0-1),"plots_sv/f10_Nneg_Npos_combined_var_inhm.png")
end
%% For  Calculating positive Flux in AEM for inhomogeneous both side i var
i__Vals = [-500,500];
for irng = 1:length(i__Vals)
    irng_0 = 0;
    if irng ==2
    irng_0 = 2;
    end
    i__Val = i__Vals(irng);
    xmesh = linspace(0.0001, 1, 212);
    td__c = 0.01;
    xi = xmesh;
    h(42+irng+irng_0-1) = figure;
    omega__c = 1;
    A_1_rng = linspace(40,400,100);
    N_poss_0_t = [];
    theta_f_p = [];
     B_1 = 0;
     B = B_1;
    for t = 1:length(A_1_rng)
        A_1 = A_1_rng(t);
        theta_cc_mul = [];
            for i = 1:length(xmesh)
                if xi(i)<=0.5
                theta_cc_mul(i)= -A_1*(xi(i)-0.5);
                else
                theta_cc_mul(i) = A_1*(xi(i)-0.5);
                end
            end
        theta_cc_mul_nd = theta_cc_mul;
        %Claculating BCs
        C__posL2 = 1;
        C__posR2 = 10;
            bcs_0 = -theta_cc_mul_nd(1)/2 + sqrt(theta_cc_mul_nd(1)^2 + 4*C__posL2^2)/2;
            bcs_1 = -theta_cc_mul_nd(end)/2 + sqrt(theta_cc_mul_nd(end)^2 + 4*C__posR2^2)/2;
            opts = bvpset('FJacobian',jac_call_theta_pe_var(A_1,B_1,omega__c,"pos_inhm",0.5,i__Val,0),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
            solinit = bvpinit(xmesh,[-0; 1]);
            sol4c = bvp4c(bvp_call_mod_theta_pe_var(A_1,B_1,omega__c,"pos_inhm",0.5,i__Val,0), bcf_call(bcs_0,bcs_1),solinit, opts);
            N_poss_0 = zeros(length(sol4c.x),1);
            theta_cc_mul_dim_plt = zeros(length(sol4c.x),1);
            for k = 1:length(sol4c.x)
                C__1 = sol4c.y(1,k);
                C__2 = sol4c.y(2,k);
                xi_vl = sol4c.x(1,k);
                if xi_vl<=0.5
                    A = -A_1;
                    N_poss_0(k) = -0.5*C__2 - 0.5*C__1*(0.5*C__2 + omega__c*A - i__Val)/((A*(xi_vl - 0.5) + B)*omega__c + 1.5*C__1);
                else
                    A = A_1;
                    N_poss_0(k) = -0.5*C__2 - 0.5*C__1*(0.5*C__2 + omega__c*A - i__Val)/((A*(xi_vl - 0.5) + B)*omega__c + 1.5*C__1);
                end
            end
            N_poss_0_t(t) = mean(N_poss_0);
            theta_f_p(t) = theta_cc_mul_nd(end);
    end
    plot(theta_f_p*td__c,N_poss_0_t,'-o','LineWidth',3)
    title("N_{pos}^{*} vs \theta for AEM")
        xlabel('\theta')
        ylabel('N_{pos}^{*}(\theta)')
        hold off;
        saveas(h(42+irng+irng_0-1),"f11_Npos_aem_theta_var_inhm_2sid.png")
    % For Calculating negative Flux
    xmesh = linspace(0.001, 1, 20);
    td__c = 0.01;
    xi_m = xmesh;
    h(43+irng+irng_0-1) = figure;
    omega__c = -1;
    %Claculating BCs
    C__posL2 = 1;
    C__posR2 = 10;
    N_poss_0_c = [];
    theta_f_n = [];
    for t = 1:length(A_1_rng)
            A_1 = A_1_rng(t);
            B_1 = 0;
            theta_cc_mul = [];
            for i = 1:length(xmesh)
                if xi_m(i)<=0.5
                theta_cc_mul(i)= -A_1*(xi_m(i)-0.5);
                else
                theta_cc_mul(i) = A_1*(xi_m(i)-0.5);
                end
            end
            theta_cc_mul_nd = theta_cc_mul;
            bcs_0 = theta_cc_mul_nd(1)/2 + sqrt(theta_cc_mul_nd(1)^2 + 4*C__posL2^2)/2;
            bcs_1 = theta_cc_mul_nd(end)/2 + sqrt(theta_cc_mul_nd(end)^2 + 4*C__posR2^2)/2;
            opts = bvpset('FJacobian',jac_call_theta_pe_var(A_1,B_1,omega__c,"neg_inhm",0.5,i__Val,0),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
            solinit = bvpinit(xmesh,[70; -60]);
            sol4c = bvp4c(bvp_call_mod_theta_pe_var(A_1,B_1,omega__c,"neg_inhm",0.5,i__Val,0), bcf_call(bcs_0,bcs_1),solinit, opts);
            N_poss_0_c = zeros(length(sol4c.x),1);
            %[sol4c.y(1,1) sol4c.y(2,1)]
            theta_cc_mul_dim_plt_N = zeros(length(sol4c.x),1);
            for k = 1:length(sol4c.x)
                C__1 = sol4c.y(1,k);
                C__2 = sol4c.y(2,k);
                xi_vl = sol4c.x(1,k);
                if xi_vl<=0.5
                    A = -A_1;
                    N_poss_0_c(k) = -0.5*C__2 - 0.5*C__1*(0.5*C__2 + omega__c*A - i__Val)/((A*(xi_vl - 0.5) + B)*omega__c + 1.5*C__1);
                else
                    A = A_1;
                    N_poss_0_c(k) = -0.5*C__2 - 0.5*C__1*(0.5*C__2 + omega__c*A - i__Val)/((A*(xi_vl - 0.5) + B)*omega__c + 1.5*C__1);
                end
            end
           N_poss_0_c(t) = mean(N_negs_0);
    end
    plot(theta_f_n*td__c,N_poss_0_c,'-o','LineWidth',3)
    title("N_{neg}^{*} vs \theta for CEM")
        xlabel('\theta')
        ylabel('N_{neg}^{*}(\theta)')
        hold off;
        saveas(h(43+irng+irng_0-1),"f11_Nneg_cem_theta_var_inhm_2sid.png")
    % Combined postive and negative flux
    h(44+irng+irng_0-1) = figure;
    plot(theta_f_p*td__c,N_poss_0_t,'-o','LineWidth',3);hold on;
    plot(theta_f_n*td__c,N_poss_0_c,'-o','LineWidth',3)
    title("N^{*} vs \theta")
        xlabel('\theta')
        ylabel('N^{*}(\theta)')
        legend("N_{pos} for AEM","N_{neg} for CEM")
        hold off;
        saveas(h(44+irng+irng_0-1),"plots_sv/f11_Nneg_Npos_combined_var_inhm_2sid.png")
end
%% Npos Final Eq
%N_poss_0(k) = -(d__pos * C__2) + d__pos * (-omega__c * A + ((-1 + d__pos) * C__2) + omega__c * (A * (xi - 0.5e0) + B) * Pe__val + i__val) * C__1 / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) + C__1 * Pe__val;
%% Nneg Final Eq
%N_negs_0(k) = (-d__pos * ((A * (xi - 0.5e0) + B) * omega__c + (2 * C__1)) * C__2 - A * C__1 * omega__c * d__pos + (Pe__val * (0.1e1 + d__pos) * C__1 - i__val) * (C__1 + (A * (xi - 0.5e0) + B) * omega__c)) / ((A * (xi - 0.5e0) + B) * omega__c + C__1 * (0.1e1 + d__pos));

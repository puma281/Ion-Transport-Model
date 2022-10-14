%% Figure Hiding
set(0,'DefaultFigureVisible','on');
%% figure off and variable clearance
clc
clear variables
h = [];
%% Initial Conditions
%sol5c = bvp5c(@bvpfcn, @bcfcn, solinit, opts);
td__c = 0.01;
%% For  Co-ion and counter ion in AEM
% For positive Flux in AEM
shp = '-o';
diff_cons = 0.5;
d_poss = diff_cons;
d_negs = 1;
d_am = 2*d_poss*d_negs/(d_poss+d_negs);
td__c = 0.01;
Pe__val = 0;
i__val = 0;
theta_cc_mul=2;
theta_cc_mul_nd = theta_cc_mul/td__c;
h(1) = figure;
omega__c = 1;
%Claculating BCs
C__posL2 = 1;
C__posR2 = 10;
N_poss = [];
    bcs_0 = -omega__c*theta_cc_mul_nd/2 + sqrt(theta_cc_mul_nd^2 + 4*C__posL2^2)/2;
    bcs_1 = -omega__c*theta_cc_mul_nd/2 + sqrt(theta_cc_mul_nd^2 + 4*C__posR2^2)/2;
    xmesh = linspace(0, 1, 20);
    opts = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_cc_mul_nd,omega__c,"pos",diff_cons,i__val,Pe__val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
    solinit = bvpinit(xmesh,[100; -0.01]);
    sol4c = bvp4c(bvp_call_mod_theta_pe_var(0,theta_cc_mul_nd,omega__c,"pos",diff_cons,i__val,Pe__val), bcf_call(bcs_0,bcs_1),solinit, opts);
    N_poss_0 = zeros(length(sol4c.x),1);
    C_Nrst_lhs_pos = zeros(length(sol4c.x),1);
    C_Nrst_rhs_pos = zeros(length(sol4c.x),1);
    d__pos = diff_cons;
    A = 0;
    B = theta_cc_mul_nd;
    for k = 1:length(sol4c.x)
        C__1 = sol4c.y(1,k);
        C__2 = sol4c.y(2,k);
        xi = sol4c.x(1,k);
      N_poss_0(k) = -(d__pos * C__2) + d__pos * (-omega__c * A + ((-1 + d__pos) * C__2) + omega__c * (A * (xi - 0.5e0) + B) * Pe__val + i__val) * C__1 / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) + C__1 * Pe__val;
    end
    for k = 1:length(sol4c.x)
        xi = sol4c.x(1,k);
      C_Nrst_lhs_pos(k) = -(N_poss_0(1)-(d_poss/(d_poss+d_negs)*i__val))/d_am*(xi-1)+C__posL2;
      C_Nrst_rhs_pos(k) = -(N_poss_0(end)-(d_poss/(d_poss+d_negs)*i__val))/d_am*(xi)+C__posR2;
    end
    C_pos_mb = sol4c.y(1,:);
    C_neg_mb = sol4c.y(1,:)+ omega__c*theta_cc_mul_nd;
    mat_dat = [C_pos_mb' C_Nrst_lhs_pos C_Nrst_rhs_pos C_neg_mb'];
%plot(sol4c.x,N_poss_0,shp,'LineWidth',3);hold on;
cocol = "#4B0082";
coucol = "#006400";
maxxi = ceil(max(max(mat_dat)));
minni = floor(min(min(mat_dat)));
if minni<0
    maxxi = maxxi+abs(minni);
end
%Bulk Plot-lhs_pos
rectangle('Position',[-2,minni,0.5,maxxi],'FaceColor',"#FFC0CB");hold on;
plot([-2,-1.5],[C_Nrst_lhs_pos(1),C_Nrst_lhs_pos(1)],'LineWidth',3,'Color',cocol);hold on;
plot([-2,-1.5],[C_Nrst_lhs_pos(1),C_Nrst_lhs_pos(1)],'LineWidth',3,'Color',coucol);hold on;
%Nerst Plot-lhs_pos
rectangle('Position',[-1.5,minni,1,maxxi],'FaceColor',"#FFFFE0");hold on;
plot(sol4c.x-1.5,C_Nrst_lhs_pos,shp,'LineWidth',3,'Color',cocol);hold on;
%Donnan Plot - lhs_pos
rectangle('Position',[-0.5,minni,0.5,maxxi],'FaceColor',"#F5FFFA");hold on;
plot([-0.5,0],[C_Nrst_lhs_pos(end),C_pos_mb(1,1)],'LineWidth',3,'Color',cocol);hold on;
plot([-0.5,0],[C_Nrst_lhs_pos(end),C_neg_mb(1,1)],'LineWidth',3,'Color',coucol);hold on;
%Charge memb
rectangle('Position',[0,minni,1,maxxi],'FaceColor',"#87CEFA");hold on;
plot(sol4c.x,C_pos_mb(1,:),shp,'LineWidth',3,'Color',cocol);hold on;
plot(sol4c.x,C_neg_mb(1,:),shp,'LineWidth',3,'Color',coucol);hold on;
%Donnan Plot - rhs_pos
rectangle('Position',[1,minni,1.5,maxxi],'FaceColor',"#F5FFFA");hold on;
plot([1,1.5],[C_pos_mb(1,end),C_Nrst_rhs_pos(1)],'LineWidth',3,'Color',cocol);hold on;
plot([1,1.5],[C_neg_mb(1,end),C_Nrst_rhs_pos(1)],'LineWidth',3,'Color',coucol);hold on;
%Nerst Plot-rhs_pos
rectangle('Position',[1.5,minni,1,maxxi],'FaceColor',"#FFFFE0");hold on;
plot(sol4c.x+1.5,C_Nrst_rhs_pos,shp,'LineWidth',3,'Color',cocol);hold on;
plot(sol4c.x+1.5,C_Nrst_rhs_pos,shp,'LineWidth',3,'Color',coucol);hold on;
%Bulk Plot-rhs_pos
rectangle('Position',[2.5,minni,0.5,maxxi],'FaceColor',"#FFC0CB");hold on;
plot([2.5,3],[C_Nrst_rhs_pos(end),C_Nrst_rhs_pos(end)],'LineWidth',3,'Color',cocol);hold on;
plot([2.5,3],[C_Nrst_rhs_pos(end),C_Nrst_rhs_pos(end)],'LineWidth',3,'Color',coucol);hold on;
%titless
title("C_{i}^{*}(\xi) vs \xi for AEM @ i^*=0,Pe =0,\theta^*=200")
    xlabel("Distance")
    ylabel('C_{i}^{*}(\xi)')
    legend("Co-ion (pos)","Counter-ion (neg)")
    if minni<0
    maxxi = maxxi-abs(minni);
   end
    ylim([minni maxxi]);
    set(gca,'XTick',[]);
    hold off;
    saveas(h(1),"plots_nrst/f2_C_i_aem_nll_co_cou.png")
%% For  Co-ion ion in AEM
% For positive Flux in AEM
shp = '-o';
diff_cons = 0.5;
d_poss = diff_cons;
d_negs = 1;
d_am = 2*d_poss*d_negs/(d_poss+d_negs);
td__c = 0.01;
Pe__val = 0;
i__val = 0;
theta_cc_mul=2;
theta_cc_mul_nd = theta_cc_mul/td__c;
h(2) = figure;
omega__c = 1;
%Claculating BCs
C__posL2 = 1;
C__posR2 = 10;
    bcs_0 = -omega__c*theta_cc_mul_nd/2 + sqrt(theta_cc_mul_nd^2 + 4*C__posL2^2)/2;
    bcs_1 = -omega__c*theta_cc_mul_nd/2 + sqrt(theta_cc_mul_nd^2 + 4*C__posR2^2)/2;
    xmesh = linspace(0, 1, 20);
    opts = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_cc_mul_nd,omega__c,"pos",diff_cons,i__val,Pe__val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
    solinit = bvpinit(xmesh,[100; -0.01]);
    sol4c = bvp4c(bvp_call_mod_theta_pe_var(0,theta_cc_mul_nd,omega__c,"pos",diff_cons,i__val,Pe__val), bcf_call(bcs_0,bcs_1),solinit, opts);
    N_poss_0 = zeros(length(sol4c.x),1);
    C_Nrst_lhs_pos = zeros(length(sol4c.x),1);
    C_Nrst_rhs_pos = zeros(length(sol4c.x),1);
    d__pos = diff_cons;
    A = 0;
    B = theta_cc_mul_nd;
    for k = 1:length(sol4c.x)
        C__1 = sol4c.y(1,k);
        C__2 = sol4c.y(2,k);
        xi = sol4c.x(1,k);
      N_poss_0(k) = -(d__pos * C__2) + d__pos * (-omega__c * A + ((-1 + d__pos) * C__2) + omega__c * (A * (xi - 0.5e0) + B) * Pe__val + i__val) * C__1 / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) + C__1 * Pe__val;
    end
    for k = 1:length(sol4c.x)
        xi = sol4c.x(1,k);
      C_Nrst_lhs_pos(k) = -(N_poss_0(1)-(d_poss/(d_poss+d_negs)*i__val))/d_am*(xi-1)+C__posL2;
      C_Nrst_rhs_pos(k) = -(N_poss_0(end)-(d_poss/(d_poss+d_negs)*i__val))/d_am*(xi)+C__posR2;
    end
    C_pos_mb = sol4c.y(1,:);
    mat_dat = [C_pos_mb' C_Nrst_lhs_pos C_Nrst_rhs_pos];
%plot(sol4c.x,N_poss_0,shp,'LineWidth',3);hold on;
cocol = "#4B0082";
maxxi = ceil(max(max(mat_dat)));
minni = floor(min(min(mat_dat)));
if minni<0
    maxxi = maxxi+abs(minni);
end
%Bulk Plot-lhs_pos
rectangle('Position',[-2,minni,0.5,maxxi],'FaceColor',"#FFC0CB");hold on;
plot([-2,-1.5],[C_Nrst_lhs_pos(1),C_Nrst_lhs_pos(1)],'LineWidth',3,'Color',cocol);hold on;
%Nerst Plot-lhs_pos
rectangle('Position',[-1.5,minni,1,maxxi],'FaceColor',"#FFFFE0");hold on;
plot(sol4c.x-1.5,C_Nrst_lhs_pos,shp,'LineWidth',3,'Color',cocol);hold on;
%Donnan Plot - lhs_pos
rectangle('Position',[-0.5,minni,0.5,maxxi],'FaceColor',"#F5FFFA");hold on;
plot([-0.5,0],[C_Nrst_lhs_pos(end),C_pos_mb(1,1)],'LineWidth',3,'Color',cocol);hold on;
%Charge memb
rectangle('Position',[0,minni,1,maxxi],'FaceColor',"#87CEFA");hold on;
plot(sol4c.x,C_pos_mb(1,:),shp,'LineWidth',3,'Color',cocol);hold on;
%Donnan Plot - rhs_pos
rectangle('Position',[1,minni,1.5,maxxi],'FaceColor',"#F5FFFA");hold on;
plot([1,1.5],[C_pos_mb(1,end),C_Nrst_rhs_pos(1)],'LineWidth',3,'Color',cocol);hold on;
%Nerst Plot-rhs_pos
rectangle('Position',[1.5,minni,1,maxxi],'FaceColor',"#FFFFE0");hold on;
plot(sol4c.x+1.5,C_Nrst_rhs_pos,shp,'LineWidth',3,'Color',cocol);hold on;
%Bulk Plot-rhs_pos
rectangle('Position',[2.5,minni,0.5,maxxi],'FaceColor',"#FFC0CB");hold on;
plot([2.5,3],[C_Nrst_rhs_pos(end),C_Nrst_rhs_pos(end)],'LineWidth',3,'Color',cocol);hold on;
%titless
title("C_{i}^{*}(\xi) vs \xi [Co-ion] for AEM @ i^*=0,Pe =0,\theta^*=200")
    xlabel("Distance")
    ylabel('C_{i}^{*}(\xi)')
    if minni<0
    maxxi = maxxi-abs(minni);
   end
    ylim([minni maxxi]);
    set(gca,'XTick',[]);
    hold off;
    saveas(h(2),"plots_nrst/f2_C_i_aem_nll_co.png")
%% For  Co-ion and counter ion in AEM i = 500 & 1000
% For positive Flux in AEM
shp = '-o';
diff_cons = 0.5;
d_poss = diff_cons;
d_negs = 1;
d_am = 2*d_poss*d_negs/(d_poss+d_negs);
td__c = 0.01;
Pe__val = 0;
theta_cc_mul=2;
theta_cc_mul_nd = theta_cc_mul/td__c;
i_rng = [0.5,1];
omega__c = 1;
for ii = 1:length(i_rng)
    h(3) = figure;
    i__val = i_rng(ii);
%Claculating BCs
C__posL2 = 1;
C__posR2 = 10;
    bcs_0 = -omega__c*theta_cc_mul_nd/2 + sqrt(theta_cc_mul_nd^2 + 4*C__posL2^2)/2;
    bcs_1 = -omega__c*theta_cc_mul_nd/2 + sqrt(theta_cc_mul_nd^2 + 4*C__posR2^2)/2;
    xmesh = linspace(0, 1, 20);
    opts = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_cc_mul_nd,omega__c,"pos",diff_cons,i__val,Pe__val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
    solinit = bvpinit(xmesh,[100; -0.01]);
    sol4c = bvp4c(bvp_call_mod_theta_pe_var(0,theta_cc_mul_nd,omega__c,"pos",diff_cons,i__val,Pe__val), bcf_call(bcs_0,bcs_1),solinit, opts);
    N_poss_0 = zeros(length(sol4c.x),1);
    C_Nrst_lhs_pos = zeros(length(sol4c.x),1);
    C_Nrst_rhs_pos = zeros(length(sol4c.x),1);
    d__pos = diff_cons;
    A = 0;
    B = theta_cc_mul_nd;
    for k = 1:length(sol4c.x)
        C__1 = sol4c.y(1,k);
        C__2 = sol4c.y(2,k);
        xi = sol4c.x(1,k);
      N_poss_0(k) = -(d__pos * C__2) + d__pos * (-omega__c * A + ((-1 + d__pos) * C__2) + omega__c * (A * (xi - 0.5e0) + B) * Pe__val + i__val) * C__1 / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) + C__1 * Pe__val;
    end
    for k = 1:length(sol4c.x)
        xi = sol4c.x(1,k);
      C_Nrst_lhs_pos(k) = -(N_poss_0(1)-(d_poss/(d_poss+d_negs)*i__val))/d_am*(xi-1)+C__posL2;
      C_Nrst_rhs_pos(k) = -(N_poss_0(end)-(d_poss/(d_poss+d_negs)*i__val))/d_am*(xi)+C__posR2;
    end
    C_pos_mb = sol4c.y(1,:);
    C_neg_mb = sol4c.y(1,:)+ omega__c*theta_cc_mul_nd;
    mat_dat = [C_pos_mb' C_Nrst_lhs_pos C_Nrst_rhs_pos C_neg_mb'];
%plot(sol4c.x,N_poss_0,shp,'LineWidth',3);hold on;
cocol = "#4B0082";
coucol = "#006400";
maxxi = ceil(max(max(mat_dat)));
minni = floor(min(min(mat_dat)));
if minni<0
    maxxi = maxxi+abs(minni);
end
%Bulk Plot-lhs_pos
rectangle('Position',[-2,minni,0.5,maxxi],'FaceColor',"#FFC0CB");hold on;
plot([-2,-1.5],[C_Nrst_lhs_pos(1),C_Nrst_lhs_pos(1)],'LineWidth',3,'Color',cocol);hold on;
plot([-2,-1.5],[C_Nrst_lhs_pos(1),C_Nrst_lhs_pos(1)],'LineWidth',3,'Color',coucol);hold on;
%Nerst Plot-lhs_pos
rectangle('Position',[-1.5,minni,1,maxxi],'FaceColor',"#FFFFE0");hold on;
plot(sol4c.x-1.5,C_Nrst_lhs_pos,shp,'LineWidth',3,'Color',cocol);hold on;
plot(sol4c.x-1.5,C_Nrst_lhs_pos,shp,'LineWidth',3,'Color',coucol);hold on;
%Donnan Plot - lhs_pos
rectangle('Position',[-0.5,minni,0.5,maxxi],'FaceColor',"#F5FFFA");hold on;
plot([-0.5,0],[C_Nrst_lhs_pos(end),C_pos_mb(1,1)],'LineWidth',3,'Color',cocol);hold on;
plot([-0.5,0],[C_Nrst_lhs_pos(end),C_neg_mb(1,1)],'LineWidth',3,'Color',coucol);hold on;
%Charge memb
rectangle('Position',[0,minni,1,maxxi],'FaceColor',"#87CEFA");hold on;
plot(sol4c.x,C_pos_mb(1,:),shp,'LineWidth',3,'Color',cocol);hold on;
plot(sol4c.x,C_neg_mb(1,:),shp,'LineWidth',3,'Color',coucol);hold on;
%Donnan Plot - rhs_pos
rectangle('Position',[1,minni,1.5,maxxi],'FaceColor',"#F5FFFA");hold on;
plot([1,1.5],[C_pos_mb(1,end),C_Nrst_rhs_pos(1)],'LineWidth',3,'Color',cocol);hold on;
plot([1,1.5],[C_neg_mb(1,end),C_Nrst_rhs_pos(1)],'LineWidth',3,'Color',coucol);hold on;
%Nerst Plot-rhs_pos
rectangle('Position',[1.5,minni,1,maxxi],'FaceColor',"#FFFFE0");hold on;
plot(sol4c.x+1.5,C_Nrst_rhs_pos,shp,'LineWidth',3,'Color',cocol);hold on;
plot(sol4c.x+1.5,C_Nrst_rhs_pos,shp,'LineWidth',3,'Color',coucol);hold on;
%Bulk Plot-rhs_pos
rectangle('Position',[2.5,minni,0.5,maxxi],'FaceColor',"#FFC0CB");hold on;
plot([2.5,3],[C_Nrst_rhs_pos(end),C_Nrst_rhs_pos(end)],'LineWidth',3,'Color',cocol);hold on;
plot([2.5,3],[C_Nrst_rhs_pos(end),C_Nrst_rhs_pos(end)],'LineWidth',3,'Color',coucol);hold on;
%titless
title("C_{i}^{*}(\xi) vs \xi for AEM for i = "+num2str(i__val))
    xlabel("Distance")
    ylabel('C_{i}^{*}(\xi)')
    legend("Co-ion (pos)","Counter-ion (neg)")
    if minni<0
    maxxi = maxxi-abs(minni);
   end
    ylim([minni maxxi]);
    set(gca,'XTick',[]);
    hold off;
    saveas(h(3),"plots_nrst/f6_Npos_aem_nll_i_"+num2str(i__val)+"_co_cou.png")
end
%% For  Co-ion in AEM i = 500 & 1000
% For positive Flux in AEM
shp = '-o';
diff_cons = 0.5;
d_poss = diff_cons;
d_negs = 1;
d_am = 2*d_poss*d_negs/(d_poss+d_negs);
td__c = 0.01;
Pe__val = 0;
theta_cc_mul=2;
theta_cc_mul_nd = theta_cc_mul/td__c;
i_rng = [0.5,1];
omega__c = 1;
for ii = 1:length(i_rng)
    h(4) = figure;
    i__val = i_rng(ii);
%Claculating BCs
C__posL2 = 1;
C__posR2 = 10;
    bcs_0 = -omega__c*theta_cc_mul_nd/2 + sqrt(theta_cc_mul_nd^2 + 4*C__posL2^2)/2;
    bcs_1 = -omega__c*theta_cc_mul_nd/2 + sqrt(theta_cc_mul_nd^2 + 4*C__posR2^2)/2;
    xmesh = linspace(0, 1, 20);
    opts = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_cc_mul_nd,omega__c,"pos",diff_cons,i__val,Pe__val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
    solinit = bvpinit(xmesh,[100; 1]);
    sol4c = bvp4c(bvp_call_mod_theta_pe_var(0,theta_cc_mul_nd,omega__c,"pos",diff_cons,i__val,Pe__val), bcf_call(bcs_0,bcs_1),solinit, opts);
    N_poss_0 = zeros(length(sol4c.x),1);
    C_Nrst_lhs_pos = zeros(length(sol4c.x),1);
    C_Nrst_rhs_pos = zeros(length(sol4c.x),1);
    d__pos = diff_cons;
    A = 0;
    B = theta_cc_mul_nd;
    for k = 1:length(sol4c.x)
        C__1 = sol4c.y(1,k);
        C__2 = sol4c.y(2,k);
        xi = sol4c.x(1,k);
      N_poss_0(k) = -(d__pos * C__2) + d__pos * (-omega__c * A + ((-1 + d__pos) * C__2) + omega__c * (A * (xi - 0.5e0) + B) * Pe__val + i__val) * C__1 / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) + C__1 * Pe__val;
    end
    for k = 1:length(sol4c.x)
        xi = sol4c.x(1,k);
      C_Nrst_lhs_pos(k) = -(N_poss_0(1)-(d_poss/(d_poss+d_negs)*i__val))/d_am*(xi-1)+C__posL2;
      C_Nrst_rhs_pos(k) = -(N_poss_0(end)-(d_poss/(d_poss+d_negs)*i__val))/d_am*(xi)+C__posR2;
    end
    C_pos_mb = sol4c.y(1,:);
    mat_dat = [C_pos_mb' C_Nrst_lhs_pos C_Nrst_rhs_pos];
%plot(sol4c.x,N_poss_0,shp,'LineWidth',3);hold on;
cocol = "#4B0082";
maxxi = ceil(max(max(mat_dat)));
minni = floor(min(min(mat_dat)));
if minni<0
    maxxi = maxxi+abs(minni);
end
%Bulk Plot-lhs_pos
rectangle('Position',[-2,minni,0.5,maxxi],'FaceColor',"#FFC0CB");hold on;
plot([-2,-1.5],[C_Nrst_lhs_pos(1),C_Nrst_lhs_pos(1)],'LineWidth',3,'Color',cocol);hold on;
%Nerst Plot-lhs_pos
rectangle('Position',[-1.5,minni,1,maxxi],'FaceColor',"#FFFFE0");hold on;
plot(sol4c.x-1.5,C_Nrst_lhs_pos,shp,'LineWidth',3,'Color',cocol);hold on;
%Donnan Plot - lhs_pos
rectangle('Position',[-0.5,minni,0.5,maxxi],'FaceColor',"#F5FFFA");hold on;
plot([-0.5,0],[C_Nrst_lhs_pos(end),C_pos_mb(1,1)],'LineWidth',3,'Color',cocol);hold on;
%Charge memb
rectangle('Position',[0,minni,1,maxxi],'FaceColor',"#87CEFA");hold on;
plot(sol4c.x,C_pos_mb(1,:),shp,'LineWidth',3,'Color',cocol);hold on;
%Donnan Plot - rhs_pos
rectangle('Position',[1,minni,1.5,maxxi],'FaceColor',"#F5FFFA");hold on;
plot([1,1.5],[C_pos_mb(1,end),C_Nrst_rhs_pos(1)],'LineWidth',3,'Color',cocol);hold on;
%Nerst Plot-rhs_pos
rectangle('Position',[1.5,minni,1,maxxi],'FaceColor',"#FFFFE0");hold on;
plot(sol4c.x+1.5,C_Nrst_rhs_pos,shp,'LineWidth',3,'Color',cocol);hold on;
%Bulk Plot-rhs_pos
rectangle('Position',[2.5,minni,0.5,maxxi],'FaceColor',"#FFC0CB");hold on;
plot([2.5,3],[C_Nrst_rhs_pos(end),C_Nrst_rhs_pos(end)],'LineWidth',3,'Color',cocol);hold on;
%titless
title("C_{i}^{*}(\xi) vs \xi for AEM for i = "+num2str(i__val))
    xlabel("Distance")
    ylabel('C_{i}^{*}(\xi)')
    if minni<0
    maxxi = maxxi-abs(minni);
    end
    ylim([minni maxxi]);
    set(gca,'XTick',[]);
    hold off;
    saveas(h(4),"plots_nrst/f6_Npos_aem_nll_i_"+num2str(i__val)+"_co.png")
end
%% For  Co-ion and counter ion in AEM for Pe = 10
% For positive Flux in AEM
shp = '-o';
diff_cons = 0.5;
d_poss = diff_cons;
d_negs = 1;
d_am = 2*d_poss*d_negs/(d_poss+d_negs);
td__c = 0.01;
Pe__val = -1;
i__val = 0;
theta_cc_mul=2;
theta_cc_mul_nd = theta_cc_mul/td__c;
h(1) = figure;
omega__c = 1;
%Claculating BCs
C__posL2 = 1;
C__posR2 = 10;
N_poss = [];
    bcs_0 = -omega__c*theta_cc_mul_nd/2 + sqrt(theta_cc_mul_nd^2 + 4*C__posL2^2)/2;
    bcs_1 = -omega__c*theta_cc_mul_nd/2 + sqrt(theta_cc_mul_nd^2 + 4*C__posR2^2)/2;
    xmesh = linspace(0, 1, 20);
    opts = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_cc_mul_nd,omega__c,"pos",diff_cons,i__val,Pe__val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
    solinit = bvpinit(xmesh,[100; -0.01]);
    sol4c = bvp4c(bvp_call_mod_theta_pe_var(0,theta_cc_mul_nd,omega__c,"pos",diff_cons,i__val,Pe__val), bcf_call(bcs_0,bcs_1),solinit, opts);
    N_poss_0 = zeros(length(sol4c.x),1);
    C_Nrst_lhs_pos = zeros(length(sol4c.x),1);
    C_Nrst_rhs_pos = zeros(length(sol4c.x),1);
    d__pos = diff_cons;
    A = 0;
    B = theta_cc_mul_nd;
    for k = 1:length(sol4c.x)
        C__1 = sol4c.y(1,k);
        C__2 = sol4c.y(2,k);
        xi = sol4c.x(1,k);
      N_poss_0(k) = -(d__pos * C__2) + d__pos * (-omega__c * A + ((-1 + d__pos) * C__2) + omega__c * (A * (xi - 0.5e0) + B) * Pe__val + i__val) * C__1 / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) + C__1 * Pe__val;
    end
    for k = 1:length(sol4c.x)
        xi = sol4c.x(1,k);
      C_Nrst_lhs_pos(k) = - (N_poss_0(1)-(d_poss/(d_poss+d_negs)*i__val))/Pe__val +(C__posL2- (N_poss_0(1)-(d_poss/(d_poss+d_negs)*i__val))/Pe__val)*exp(Pe__val/d_am*(xi-1));
      C_Nrst_rhs_pos(k) = - (N_poss_0(end)-(d_poss/(d_poss+d_negs)*i__val))/Pe__val +(C__posR2- (N_poss_0(end)-(d_poss/(d_poss+d_negs)*i__val))/Pe__val)*exp(Pe__val/d_am*(xi));
    end
      C_pos_mb = sol4c.y(1,:);
    C_neg_mb = sol4c.y(1,:)+ omega__c*theta_cc_mul_nd;
    mat_dat = [C_pos_mb' C_Nrst_lhs_pos C_Nrst_rhs_pos C_neg_mb'];
%plot(sol4c.x,N_poss_0,shp,'LineWidth',3);hold on;
cocol = "#4B0082";
coucol = "#006400";
maxxi = ceil(max(max(mat_dat)));
minni = floor(min(min(mat_dat)));
if minni<0
    maxxi = maxxi+abs(minni);
end
%Bulk Plot-lhs_pos
rectangle('Position',[-2,minni,0.5,maxxi],'FaceColor',"#FFC0CB");hold on;
plot([-2,-1.5],[C_Nrst_lhs_pos(1),C_Nrst_lhs_pos(1)],'LineWidth',3,'Color',cocol);hold on;
plot([-2,-1.5],[C_Nrst_lhs_pos(1),C_Nrst_lhs_pos(1)],'LineWidth',3,'Color',coucol);hold on;
%Nerst Plot-lhs_pos
rectangle('Position',[-1.5,minni,1,maxxi],'FaceColor',"#FFFFE0");hold on;
plot(sol4c.x-1.5,C_Nrst_lhs_pos,shp,'LineWidth',3,'Color',cocol);hold on;
%Donnan Plot - lhs_pos
rectangle('Position',[-0.5,minni,0.5,maxxi],'FaceColor',"#F5FFFA");hold on;
plot([-0.5,0],[C_Nrst_lhs_pos(end),C_pos_mb(1,1)],'LineWidth',3,'Color',cocol);hold on;
plot([-0.5,0],[C_Nrst_lhs_pos(end),C_neg_mb(1,1)],'LineWidth',3,'Color',coucol);hold on;
%Charge memb
rectangle('Position',[0,minni,1,maxxi],'FaceColor',"#87CEFA");hold on;
plot(sol4c.x,C_pos_mb(1,:),shp,'LineWidth',3,'Color',cocol);hold on;
plot(sol4c.x,C_neg_mb(1,:),shp,'LineWidth',3,'Color',coucol);hold on;
%Donnan Plot - rhs_pos
rectangle('Position',[1,minni,1.5,maxxi],'FaceColor',"#F5FFFA");hold on;
plot([1,1.5],[C_pos_mb(1,end),C_Nrst_rhs_pos(1)],'LineWidth',3,'Color',cocol);hold on;
plot([1,1.5],[C_neg_mb(1,end),C_Nrst_rhs_pos(1)],'LineWidth',3,'Color',coucol);hold on;
%Nerst Plot-rhs_pos
rectangle('Position',[1.5,minni,1,maxxi],'FaceColor',"#FFFFE0");hold on;
plot(sol4c.x+1.5,C_Nrst_rhs_pos,shp,'LineWidth',3,'Color',cocol);hold on;
plot(sol4c.x+1.5,C_Nrst_rhs_pos,shp,'LineWidth',3,'Color',coucol);hold on;
%Bulk Plot-rhs_pos
rectangle('Position',[2.5,minni,0.5,maxxi],'FaceColor',"#FFC0CB");hold on;
plot([2.5,3],[C_Nrst_rhs_pos(end),C_Nrst_rhs_pos(end)],'LineWidth',3,'Color',cocol);hold on;
plot([2.5,3],[C_Nrst_rhs_pos(end),C_Nrst_rhs_pos(end)],'LineWidth',3,'Color',coucol);hold on;
%titless
title("C_{i}^{*}(\xi) vs \xi for AEM at Pe = " +num2str(abs(Pe__val)))
    xlabel("Distance")
    ylabel('C_{i}^{*}(\xi)')
    legend("Co-ion (pos)","Counter-ion (neg)")
    if minni<0
    maxxi = maxxi-abs(minni);
   end
    ylim([minni maxxi]);
    set(gca,'XTick',[]);
    hold off;
    saveas(h(1),"plots_nrst/S2_C_i_aem_nll_co_cou_Pe.png")
%% For  Co-ion ion in AEM for Pe  = 10
% For positive Flux in AEM
shp = '-o';
diff_cons = 0.5;
d_poss = diff_cons;
d_negs = 1;
d_am = 2*d_poss*d_negs/(d_poss+d_negs);
td__c = 0.01;
Pe__val = -1;
i__val = 0;
theta_cc_mul=2;
theta_cc_mul_nd = theta_cc_mul/td__c;
h(2) = figure;
omega__c = 1;
%Claculating BCs
C__posL2 = 1;
C__posR2 = 10;
    bcs_0 = -omega__c*theta_cc_mul_nd/2 + sqrt(theta_cc_mul_nd^2 + 4*C__posL2^2)/2;
    bcs_1 = -omega__c*theta_cc_mul_nd/2 + sqrt(theta_cc_mul_nd^2 + 4*C__posR2^2)/2;
    xmesh = linspace(0, 1, 20);
    opts = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_cc_mul_nd,omega__c,"pos",diff_cons,i__val,Pe__val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
    solinit = bvpinit(xmesh,[100; -0.01]);
    sol4c = bvp4c(bvp_call_mod_theta_pe_var(0,theta_cc_mul_nd,omega__c,"pos",diff_cons,i__val,Pe__val), bcf_call(bcs_0,bcs_1),solinit, opts);
    N_poss_0 = zeros(length(sol4c.x),1);
    C_Nrst_lhs_pos = zeros(length(sol4c.x),1);
    C_Nrst_rhs_pos = zeros(length(sol4c.x),1);
    d__pos = diff_cons;
    A = 0;
    B = theta_cc_mul_nd;
    for k = 1:length(sol4c.x)
        C__1 = sol4c.y(1,k);
        C__2 = sol4c.y(2,k);
        xi = sol4c.x(1,k);
      N_poss_0(k) = -(d__pos * C__2) + d__pos * (-omega__c * A + ((-1 + d__pos) * C__2) + omega__c * (A * (xi - 0.5e0) + B) * Pe__val + i__val) * C__1 / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) + C__1 * Pe__val;
    end
        for k = 1:length(sol4c.x)
        xi = sol4c.x(1,k);
      C_Nrst_lhs_pos(k) = - (N_poss_0(1)-(d_poss/(d_poss+d_negs)*i__val))/Pe__val +(C__posL2- (N_poss_0(1)-(d_poss/(d_poss+d_negs)*i__val))/Pe__val)*exp(Pe__val/d_am*(xi-1));
      C_Nrst_rhs_pos(k) = - (N_poss_0(end)-(d_poss/(d_poss+d_negs)*i__val))/Pe__val +(C__posR2- (N_poss_0(end)-(d_poss/(d_poss+d_negs)*i__val))/Pe__val)*exp(Pe__val/d_am*(xi));
        end
    C_pos_mb = sol4c.y(1,:);
    %Dim-
C_pos_mb = 10*C_pos_mb;
C_Nrst_lhs_pos = 10*C_Nrst_lhs_pos;
C_Nrst_rhs_pos = 10*C_Nrst_rhs_pos;
    mat_dat = [C_pos_mb' C_Nrst_lhs_pos C_Nrst_rhs_pos];
%plot(sol4c.x,N_poss_0,shp,'LineWidth',3);hold on;
cocol = "#4B0082";
maxxi = ceil(max(max(mat_dat)));
minni = floor(min(min(mat_dat)));
if minni<0
    maxxi = maxxi+abs(minni);
end
%Bulk Plot-lhs_pos
rectangle('Position',[-2,minni,0.5,maxxi],'FaceColor',"#FFC0CB");hold on;
plot([-2,-1.5],[C_Nrst_lhs_pos(1),C_Nrst_lhs_pos(1)],'LineWidth',3,'Color',cocol);hold on;
%Nerst Plot-lhs_pos
rectangle('Position',[-1.5,minni,1,maxxi],'FaceColor',"#FFFFE0");hold on;
plot(sol4c.x-1.5,C_Nrst_lhs_pos,shp,'LineWidth',3,'Color',cocol);hold on;
%Donnan Plot - lhs_pos
rectangle('Position',[-0.5,minni,0.5,maxxi],'FaceColor',"#F5FFFA");hold on;
plot([-0.5,0],[C_Nrst_lhs_pos(end),C_pos_mb(1,1)],'LineWidth',3,'Color',cocol);hold on;
%Charge memb
rectangle('Position',[0,minni,1,maxxi],'FaceColor',"#87CEFA");hold on;
plot(sol4c.x,C_pos_mb(1,:),shp,'LineWidth',3,'Color',cocol);hold on;
%Donnan Plot - rhs_pos
rectangle('Position',[1,minni,1.5,maxxi],'FaceColor',"#F5FFFA");hold on;
plot([1,1.5],[C_pos_mb(1,end),C_Nrst_rhs_pos(1)],'LineWidth',3,'Color',cocol);hold on;
%Nerst Plot-rhs_pos
rectangle('Position',[1.5,minni,1,maxxi],'FaceColor',"#FFFFE0");hold on;
plot(sol4c.x+1.5,C_Nrst_rhs_pos,shp,'LineWidth',3,'Color',cocol);hold on;
%Bulk Plot-rhs_pos
rectangle('Position',[2.5,minni,0.5,maxxi],'FaceColor',"#FFC0CB");hold on;
plot([2.5,3],[C_Nrst_rhs_pos(end),C_Nrst_rhs_pos(end)],'LineWidth',3,'Color',cocol);hold on;
%titless
title("C_{i}^{*}(\xi) vs \xi [Co-ion] for AEM at Pe = " +num2str(abs(Pe__val)))
    xlabel("Distance")
    ylabel('C_{i}^{*}(\xi)')
    if minni<0
    maxxi = maxxi-abs(minni);
   end
    ylim([minni maxxi]);
    set(gca,'XTick',[]);
    hold off;
    saveas(h(2),"plots_nrst/S2_C_i_aem_nll_co.png")
%% For  Co-ion ion in AEM diff Pe = 10 & 100
% For positive Flux in AEM
shp = '-o';
diff_cons = 1;
d_poss = diff_cons;
d_negs = 1;
d_am = 2*d_poss*d_negs/(d_poss+d_negs);
td__c = 0.01;
i__val = 1000;
theta_cc_mul=2;
theta_cc_mul_nd = theta_cc_mul/td__c;
omega__rng = [1,-1];
omega_nms = ["AEM", "CEM"];
for pp = 1:length(omega__rng)
h(2) = figure;
omega__c = omega__rng(pp);
if omega__c == 1
    Pe__val = 0.8;
else
    Pe__val = -0.8;
end
%Claculating BCs
C__posL2 = 1;
C__posR2 = 10;
    bcs_0 = -omega__c*theta_cc_mul_nd/2 + sqrt(theta_cc_mul_nd^2 + 4*C__posL2^2)/2;
    bcs_1 = -omega__c*theta_cc_mul_nd/2 + sqrt(theta_cc_mul_nd^2 + 4*C__posR2^2)/2;
    xmesh = linspace(0, 1, 20);
    opts = bvpset('FJacobian',jac_call_theta_pe_var(0,theta_cc_mul_nd,omega__c,"pos",diff_cons,i__val,Pe__val),'RelTol',0.1,'AbsTol',0.1,'Stats','on');
    solinit = bvpinit(xmesh,[10; -0.01]);
    sol4c = bvp4c(bvp_call_mod_theta_pe_var(0,theta_cc_mul_nd,omega__c,"pos",diff_cons,i__val,Pe__val), bcf_call(bcs_0,bcs_1),solinit, opts);
    N_poss_0 = zeros(length(sol4c.x),1);
    C_Nrst_lhs_pos = zeros(length(sol4c.x),1);
    C_Nrst_rhs_pos = zeros(length(sol4c.x),1);
    d__pos = diff_cons;
    A = 0;
    B = theta_cc_mul_nd;
    for k = 1:length(sol4c.x)
        C__1 = sol4c.y(1,k);
        C__2 = sol4c.y(2,k);
        xi = sol4c.x(1,k);
      N_poss_0(k) = -(d__pos * C__2) + d__pos * (-omega__c * A + ((-1 + d__pos) * C__2) + omega__c * (A * (xi - 0.5e0) + B) * Pe__val + i__val) * C__1 / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) + C__1 * Pe__val;
    end
       for k = 1:length(sol4c.x)
        xi = sol4c.x(1,k);
      C_Nrst_lhs_pos(k) = - (N_poss_0(1)-(d_poss/(d_poss+d_negs)*i__val))/Pe__val +(C__posL2- (N_poss_0(1)-(d_poss/(d_poss+d_negs)*i__val))/Pe__val)*exp(Pe__val/d_am*(xi-1));
      C_Nrst_rhs_pos(k) = - (N_poss_0(end)-(d_poss/(d_poss+d_negs)*i__val))/Pe__val +(C__posR2- (N_poss_0(end)-(d_poss/(d_poss+d_negs)*i__val))/Pe__val)*exp(Pe__val/d_am*(xi));
        end
    C_pos_mb = sol4c.y(1,:);
    mat_dat = [C_pos_mb' C_Nrst_lhs_pos C_Nrst_rhs_pos];
%plot(sol4c.x,N_poss_0,shp,'LineWidth',3);hold on;
cocol = "#4B0082";
maxxi = ceil(max(max(mat_dat)));
minni = floor(min(min(mat_dat)));
if minni<0
    maxxi = maxxi+abs(minni);
end
%Bulk Plot-lhs_pos
rectangle('Position',[-2,minni,0.5,maxxi],'FaceColor',"#FFC0CB");hold on;
plot([-2,-1.5],[C_Nrst_lhs_pos(1),C_Nrst_lhs_pos(1)],'LineWidth',3,'Color',cocol);hold on;
%Nerst Plot-lhs_pos
rectangle('Position',[-1.5,minni,1,maxxi],'FaceColor',"#FFFFE0");hold on;
plot(sol4c.x-1.5,C_Nrst_lhs_pos,shp,'LineWidth',3,'Color',cocol);hold on;
%Donnan Plot - lhs_pos
rectangle('Position',[-0.5,minni,0.5,maxxi],'FaceColor',"#F5FFFA");hold on;
plot([-0.5,0],[C_Nrst_lhs_pos(end),C_pos_mb(1,1)],'LineWidth',3,'Color',cocol);hold on;
%Charge memb
rectangle('Position',[0,minni,1,maxxi],'FaceColor',"#87CEFA");hold on;
plot(sol4c.x,C_pos_mb(1,:),shp,'LineWidth',3,'Color',cocol);hold on;
%Donnan Plot - rhs_pos
rectangle('Position',[1,minni,1.5,maxxi],'FaceColor',"#F5FFFA");hold on;
plot([1,1.5],[C_pos_mb(1,end),C_Nrst_rhs_pos(1)],'LineWidth',3,'Color',cocol);hold on;
%Nerst Plot-rhs_pos
rectangle('Position',[1.5,minni,1,maxxi],'FaceColor',"#FFFFE0");hold on;
plot(sol4c.x+1.5,C_Nrst_rhs_pos,shp,'LineWidth',3,'Color',cocol);hold on;
%Bulk Plot-rhs_pos
rectangle('Position',[2.5,minni,0.5,maxxi],'FaceColor',"#FFC0CB");hold on;
plot([2.5,3],[C_Nrst_rhs_pos(end),C_Nrst_rhs_pos(end)],'LineWidth',3,'Color',cocol);hold on;
%titless
title("C_{i}^{*}(\xi) vs \xi [Co-ion] for " +omega_nms(pp))
    xlabel("Distance")
    ylabel('C_{i}^{*}(\xi)')
    if minni<0
    maxxi = maxxi-abs(minni);
   end
    ylim([minni maxxi]);
    set(gca,'XTick',[]);
    hold off;
    saveas(h(2),"plots_nrst/S5_C_i_aem_nll_co"+omega_nms(pp)+".png")
end
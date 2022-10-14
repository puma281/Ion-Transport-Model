function bvp = bvp_call_mod_theta_pe_var(A_1,B_1,omega__c,ion_nm,d__pos,i__val,Pe__val)
    if ion_nm == "pos"
        bvp = @bvpfcn_pos;
    elseif ion_nm == "pos_inhm"     
        bvp = @bvpfcn_pos_var;
    end
    function dCdxi = bvpfcn_pos(xi,C)
            A = A_1;
            B = B_1;
            C__2 = C(2);
            C__1 = C(1);
            dCdxi_1 = C__2;
            dCdxi_2 = (-omega__c * d__pos * (A * (xi - 0.5e0) + B) * (1 - d__pos) * C__2 ^ 2 + (-d__pos * (((-1 + d__pos) * C__1) + (A * (xi - 0.5e0) + B) * omega__c) * omega__c * A + (Pe__val * (1 + d__pos) ^ 2 * C__1 ^ 2) + 0.2e1 * Pe__val * omega__c * (A * (xi - 0.5e0) + B) * (1 + d__pos) * C__1 + (Pe__val * omega__c * (1 + d__pos) * (A * (xi - 0.5e0) + B) + (i__val * d__pos)) * omega__c * (A * (xi - 0.5e0) + B)) * C__2 + d__pos * C__1 * A * omega__c * (omega__c * A + (Pe__val * (1 + d__pos) * C__1) - i__val)) / d__pos / ((A * (xi - 0.5e0) + B) * omega__c + (2 * C__1)) / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos)));
            dCdxi = [dCdxi_1 
                     dCdxi_2];
    end   
    function dCdxi = bvpfcn_pos_var(xi,C)
            C__2 = C(2);
            C__1 = C(1);
            dCdxi_1 = C__2;
            if xi<=0.5
                A = -A_1;
                B = 0;
                dCdxi_2 = (-omega__c * d__pos * (A * (xi - 0.5e0) + B) * (1 - d__pos) * C__2 ^ 2 + (-d__pos * (((-1 + d__pos) * C__1) + (A * (xi - 0.5e0) + B) * omega__c) * omega__c * A + (Pe__val * (1 + d__pos) ^ 2 * C__1 ^ 2) + 0.2e1 * Pe__val * omega__c * (A * (xi - 0.5e0) + B) * (1 + d__pos) * C__1 + (Pe__val * omega__c * (1 + d__pos) * (A * (xi - 0.5e0) + B) + (i__val * d__pos)) * omega__c * (A * (xi - 0.5e0) + B)) * C__2 + d__pos * C__1 * A * omega__c * (omega__c * A + (Pe__val * (1 + d__pos) * C__1) - i__val)) / d__pos / ((A * (xi - 0.5e0) + B) * omega__c + (2 * C__1)) / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos)));
            else
                A = A_1;
                B = 0;
                dCdxi_2 = (-omega__c * d__pos * (A * (xi - 0.5e0) + B) * (1 - d__pos) * C__2 ^ 2 + (-d__pos * (((-1 + d__pos) * C__1) + (A * (xi - 0.5e0) + B) * omega__c) * omega__c * A + (Pe__val * (1 + d__pos) ^ 2 * C__1 ^ 2) + 0.2e1 * Pe__val * omega__c * (A * (xi - 0.5e0) + B) * (1 + d__pos) * C__1 + (Pe__val * omega__c * (1 + d__pos) * (A * (xi - 0.5e0) + B) + (i__val * d__pos)) * omega__c * (A * (xi - 0.5e0) + B)) * C__2 + d__pos * C__1 * A * omega__c * (omega__c * A + (Pe__val * (1 + d__pos) * C__1) - i__val)) / d__pos / ((A * (xi - 0.5e0) + B) * omega__c + (2 * C__1)) / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos)));
            end
            dCdxi = [dCdxi_1 
                     dCdxi_2];
    end  
    
end
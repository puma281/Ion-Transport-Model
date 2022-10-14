function jac = jac_call_theta_pe_var(A_1,B_1,omega__c,ion_nm,d__pos,i__val,Pe__val)
if ion_nm == "pos"
    jac = @jac_pos;
elseif ion_nm == "pos_inhm"
    jac = @jac_pos_var;
end
    function dfdC = jac_pos(xi,C)
    A = A_1;
    B = B_1;
    C__2 = C(2);
    C__1 = C(1);
    jac_pos_3 = ((-(d__pos * (-1 + d__pos) * omega__c * A) + (2 * Pe__val * (1 + d__pos) ^ 2 * C__1) + 0.2e1 * Pe__val * omega__c * (1 + d__pos) * (A * (xi - 0.5e0) + B)) * C__2 + (d__pos * A * omega__c * (omega__c * A + Pe__val * (1 + d__pos) * C__1 - i__val)) + (d__pos * C__1 * A * omega__c * Pe__val * (1 + d__pos))) / d__pos / ((A * (xi - 0.5e0) + B) * omega__c + (2 * C__1)) / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) - 0.2e1 * (-omega__c * d__pos * (A * (xi - 0.5e0) + B) * (1 - d__pos) * C__2 ^ 2 + (-d__pos * (((-1 + d__pos) * C__1) + (A * (xi - 0.5e0) + B) * omega__c) * omega__c * A + (Pe__val * (1 + d__pos) ^ 2 * C__1 ^ 2) + 0.2e1 * Pe__val * omega__c * (A * (xi - 0.5e0) + B) * (1 + d__pos) * C__1 + (Pe__val * omega__c * (1 + d__pos) * (A * (xi - 0.5e0) + B) + (i__val * d__pos)) * omega__c * (A * (xi - 0.5e0) + B)) * C__2 + (d__pos * C__1 * A * omega__c * (omega__c * A + Pe__val * (1 + d__pos) * C__1 - i__val))) / d__pos / ((A * (xi - 0.5e0) + B) * omega__c + (2 * C__1)) ^ 2 / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) - (-omega__c * d__pos * (A * (xi - 0.5e0) + B) * (1 - d__pos) * C__2 ^ 2 + (-d__pos * (((-1 + d__pos) * C__1) + (A * (xi - 0.5e0) + B) * omega__c) * omega__c * A + (Pe__val * (1 + d__pos) ^ 2 * C__1 ^ 2) + 0.2e1 * Pe__val * omega__c * (A * (xi - 0.5e0) + B) * (1 + d__pos) * C__1 + (Pe__val * omega__c * (1 + d__pos) * (A * (xi - 0.5e0) + B) + (i__val * d__pos)) * omega__c * (A * (xi - 0.5e0) + B)) * C__2 + (d__pos * C__1 * A * omega__c * (omega__c * A + Pe__val * (1 + d__pos) * C__1 - i__val))) / d__pos / ((A * (xi - 0.5e0) + B) * omega__c + (2 * C__1)) / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) ^ 2 * (1 + d__pos);
    jac_pos_4 = (-0.2e1 * omega__c * d__pos * (A * (xi - 0.5e0) + B) * (1 - d__pos) * C__2 - d__pos * (((-1 + d__pos) * C__1) + (A * (xi - 0.5e0) + B) * omega__c) * omega__c * A + (Pe__val * (1 + d__pos) ^ 2 * C__1 ^ 2) + 0.2e1 * Pe__val * omega__c * (A * (xi - 0.5e0) + B) * (1 + d__pos) * C__1 + (Pe__val * omega__c * (1 + d__pos) * (A * (xi - 0.5e0) + B) + (i__val * d__pos)) * omega__c * (A * (xi - 0.5e0) + B)) / d__pos / ((A * (xi - 0.5e0) + B) * omega__c + (2 * C__1)) / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos)));
    dfdC = [0      1
           jac_pos_3    jac_pos_4];
    end
    function dfdC = jac_pos_var(xi,C)
    C__2 = C(2);
    C__1 = C(1);
    if xi<=0.5
        A = -A_1;
        B = 0;
    jac_pos_3 = ((-(d__pos * (-1 + d__pos) * omega__c * A) + (2 * Pe__val * (1 + d__pos) ^ 2 * C__1) + 0.2e1 * Pe__val * omega__c * (1 + d__pos) * (A * (xi - 0.5e0) + B)) * C__2 + (d__pos * A * omega__c * (omega__c * A + Pe__val * (1 + d__pos) * C__1 - i__val)) + (d__pos * C__1 * A * omega__c * Pe__val * (1 + d__pos))) / d__pos / ((A * (xi - 0.5e0) + B) * omega__c + (2 * C__1)) / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) - 0.2e1 * (-omega__c * d__pos * (A * (xi - 0.5e0) + B) * (1 - d__pos) * C__2 ^ 2 + (-d__pos * (((-1 + d__pos) * C__1) + (A * (xi - 0.5e0) + B) * omega__c) * omega__c * A + (Pe__val * (1 + d__pos) ^ 2 * C__1 ^ 2) + 0.2e1 * Pe__val * omega__c * (A * (xi - 0.5e0) + B) * (1 + d__pos) * C__1 + (Pe__val * omega__c * (1 + d__pos) * (A * (xi - 0.5e0) + B) + (i__val * d__pos)) * omega__c * (A * (xi - 0.5e0) + B)) * C__2 + (d__pos * C__1 * A * omega__c * (omega__c * A + Pe__val * (1 + d__pos) * C__1 - i__val))) / d__pos / ((A * (xi - 0.5e0) + B) * omega__c + (2 * C__1)) ^ 2 / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) - (-omega__c * d__pos * (A * (xi - 0.5e0) + B) * (1 - d__pos) * C__2 ^ 2 + (-d__pos * (((-1 + d__pos) * C__1) + (A * (xi - 0.5e0) + B) * omega__c) * omega__c * A + (Pe__val * (1 + d__pos) ^ 2 * C__1 ^ 2) + 0.2e1 * Pe__val * omega__c * (A * (xi - 0.5e0) + B) * (1 + d__pos) * C__1 + (Pe__val * omega__c * (1 + d__pos) * (A * (xi - 0.5e0) + B) + (i__val * d__pos)) * omega__c * (A * (xi - 0.5e0) + B)) * C__2 + (d__pos * C__1 * A * omega__c * (omega__c * A + Pe__val * (1 + d__pos) * C__1 - i__val))) / d__pos / ((A * (xi - 0.5e0) + B) * omega__c + (2 * C__1)) / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) ^ 2 * (1 + d__pos);
    jac_pos_4 = (-0.2e1 * omega__c * d__pos * (A * (xi - 0.5e0) + B) * (1 - d__pos) * C__2 - d__pos * (((-1 + d__pos) * C__1) + (A * (xi - 0.5e0) + B) * omega__c) * omega__c * A + (Pe__val * (1 + d__pos) ^ 2 * C__1 ^ 2) + 0.2e1 * Pe__val * omega__c * (A * (xi - 0.5e0) + B) * (1 + d__pos) * C__1 + (Pe__val * omega__c * (1 + d__pos) * (A * (xi - 0.5e0) + B) + (i__val * d__pos)) * omega__c * (A * (xi - 0.5e0) + B)) / d__pos / ((A * (xi - 0.5e0) + B) * omega__c + (2 * C__1)) / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos)));
    else
        A = A_1;
        B = 0;
    jac_pos_3 = ((-(d__pos * (-1 + d__pos) * omega__c * A) + (2 * Pe__val * (1 + d__pos) ^ 2 * C__1) + 0.2e1 * Pe__val * omega__c * (1 + d__pos) * (A * (xi - 0.5e0) + B)) * C__2 + (d__pos * A * omega__c * (omega__c * A + Pe__val * (1 + d__pos) * C__1 - i__val)) + (d__pos * C__1 * A * omega__c * Pe__val * (1 + d__pos))) / d__pos / ((A * (xi - 0.5e0) + B) * omega__c + (2 * C__1)) / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) - 0.2e1 * (-omega__c * d__pos * (A * (xi - 0.5e0) + B) * (1 - d__pos) * C__2 ^ 2 + (-d__pos * (((-1 + d__pos) * C__1) + (A * (xi - 0.5e0) + B) * omega__c) * omega__c * A + (Pe__val * (1 + d__pos) ^ 2 * C__1 ^ 2) + 0.2e1 * Pe__val * omega__c * (A * (xi - 0.5e0) + B) * (1 + d__pos) * C__1 + (Pe__val * omega__c * (1 + d__pos) * (A * (xi - 0.5e0) + B) + (i__val * d__pos)) * omega__c * (A * (xi - 0.5e0) + B)) * C__2 + (d__pos * C__1 * A * omega__c * (omega__c * A + Pe__val * (1 + d__pos) * C__1 - i__val))) / d__pos / ((A * (xi - 0.5e0) + B) * omega__c + (2 * C__1)) ^ 2 / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) - (-omega__c * d__pos * (A * (xi - 0.5e0) + B) * (1 - d__pos) * C__2 ^ 2 + (-d__pos * (((-1 + d__pos) * C__1) + (A * (xi - 0.5e0) + B) * omega__c) * omega__c * A + (Pe__val * (1 + d__pos) ^ 2 * C__1 ^ 2) + 0.2e1 * Pe__val * omega__c * (A * (xi - 0.5e0) + B) * (1 + d__pos) * C__1 + (Pe__val * omega__c * (1 + d__pos) * (A * (xi - 0.5e0) + B) + (i__val * d__pos)) * omega__c * (A * (xi - 0.5e0) + B)) * C__2 + (d__pos * C__1 * A * omega__c * (omega__c * A + Pe__val * (1 + d__pos) * C__1 - i__val))) / d__pos / ((A * (xi - 0.5e0) + B) * omega__c + (2 * C__1)) / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos))) ^ 2 * (1 + d__pos);
    jac_pos_4 = (-0.2e1 * omega__c * d__pos * (A * (xi - 0.5e0) + B) * (1 - d__pos) * C__2 - d__pos * (((-1 + d__pos) * C__1) + (A * (xi - 0.5e0) + B) * omega__c) * omega__c * A + (Pe__val * (1 + d__pos) ^ 2 * C__1 ^ 2) + 0.2e1 * Pe__val * omega__c * (A * (xi - 0.5e0) + B) * (1 + d__pos) * C__1 + (Pe__val * omega__c * (1 + d__pos) * (A * (xi - 0.5e0) + B) + (i__val * d__pos)) * omega__c * (A * (xi - 0.5e0) + B)) / d__pos / ((A * (xi - 0.5e0) + B) * omega__c + (2 * C__1)) / ((A * (xi - 0.5e0) + B) * omega__c + (C__1 * (1 + d__pos)));
    end 
    dfdC = [0      1
           jac_pos_3    jac_pos_4];
    end
end
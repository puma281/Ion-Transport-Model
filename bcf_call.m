function bcf = bcf_call(bcs_0,bcs_1)
    bcf = @bcfcn;
function res = bcfcn(Ca,Cb)
    res = [Ca(1)-bcs_0 
           Cb(1)-bcs_1];
    end
end
  

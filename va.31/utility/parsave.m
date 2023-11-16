function parsave(fname,data)

var_name=genvarname(inputname(2));
eval([var_name '=data',';']);

try
    save(fname,var_name,'-append');
catch
    save(fname,var_name,'-v7.3');
end
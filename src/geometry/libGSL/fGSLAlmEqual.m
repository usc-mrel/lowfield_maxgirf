function result = fGSLAlmEqual(dArg1, dArg2) % "Alm"ost "Equal"
%
% Descrip: This function compares two double variables and gives  
%           back the value true, if the values are equal to five   
%           digits after the decimal point. The two double         
%           variables must not exceed +/- 1.0! 
%
%--------------------------------------------------------------------------
% dArg1    IMP: first argument
% dArg2    IMP: second argument
%--------------------------------------------------------------------------

dTmp = dArg1 - dArg2;
result = (dTmp >= -1.0e-6) && (dTmp <= 1.0e-6);

end
function sig = chi2sig(chi2, dof)
% function sig = chi2sig(chi2, dof)
sig = redchi2sig(chi2 / dof, dof);
return

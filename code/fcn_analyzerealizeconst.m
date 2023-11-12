function [Phi,Gamma,Np1,Nm1] = fcn_analyzerealizeconst(F)
[Phi,Gamma]=eig(F);
[gamma,gidx]=sort(diag(Gamma),'descend');
Phi = Phi(:,gidx);
Gamma = diag(gamma);
Np1 = sum(abs(gamma-1)<1e-2);
Nm1 = sum(abs(gamma+1)<1e-2);
end
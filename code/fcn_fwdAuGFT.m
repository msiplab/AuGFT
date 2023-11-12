function X = fcn_fwdAuGFT(x,U,Q)
% Forward AuGFT
I = eye(size(U));
alpha = (1-1/sqrt(2));
fwdAuGFT = [(I-alpha*(Q*Q.'))*U Q/sqrt(2)].';
X = fwdAuGFT*x;
end
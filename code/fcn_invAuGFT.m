function x = fcn_invAuGFT(X,U,Q)
% Inverse AuGFT
I = eye(size(U));
alpha = (1-1/sqrt(2));
invAuGFT = [(I-alpha*(Q*Q.'))*U Q/sqrt(2)];
x = invAuGFT*X;
end
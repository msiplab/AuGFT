import numpy as np

def fcn_fwdAuGFT(x, U, Q, alpha):
    if alpha is None:
        alpha = 1 - 1 / np.sqrt(2)

    I = np.eye(U.shape[0])
    beta = np.sqrt(alpha * (2 - alpha))
    fwdAuGFT = np.vstack([(I - alpha * np.dot(Q, Q.T)) @ U, beta * Q]).T
    X = np.dot(fwdAuGFT, x)
    return X

"""
function X = fcn_fwdAuGFT(x,U,Q,alpha) %beta)

if nargin < 4
    %beta = 1/sqrt(2);
    alpha = 1-1/sqrt(2);
end

% Forward AuGFT
I = eye(size(U));
%alpha = (1-sqrt(1-beta^2));
beta = sqrt(alpha*(2-alpha));
fwdAuGFT = [(I-alpha*(Q*Q.'))*U beta*Q].';
X = fwdAuGFT*x;
end
"""
import numpy as np

def fcn_invAuGFT(X, U, Q, alpha):
    if len(alpha) < 4:
        alpha = 1 - 1 / np.sqrt(2)

    I = np.eye(U.shape[0])
    beta = np.sqrt(alpha * (2 - alpha))
    invAuGFT = np.concatenate([(I - alpha * np.dot(Q, Q.T)) * U, beta * Q], axis=1)
    x = np.dot(invAuGFT, X)
    return x

"""
function x = fcn_invAuGFT(X,U,Q,alpha) %beta)

if nargin < 4
    % beta = 1/sqrt(2);
    alpha = 1-1/sqrt(2);
end

% Inverse AuGFT
I = eye(size(U));
%alpha = (1-sqrt(1-beta^2));
beta = sqrt(alpha*(2-alpha));
invAuGFT = [(I-alpha*(Q*Q.'))*U beta*Q];
x = invAuGFT*X;
end
"""
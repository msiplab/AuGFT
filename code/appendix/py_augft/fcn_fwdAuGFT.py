import numpy as np

def fcn_fwdAuGFT(x, U, Q, **kwaargs):

    if 'alpha' in kwargs:
        alpha = kwargs.get('alpha')
    else:
        alpha = 1 - 1 / np.sqrt(2)

    # Forward AuGFT
    I = np.eye(U.shape[0])
    beta = np.sqrt(alpha * (2 - alpha))
    fwdAuGFT = np.vstack(((I - alpha * np.dot(Q, Q.T)) @ U, beta * Q)).T
    X = np.dot(fwdAuGFT, x)
    return X

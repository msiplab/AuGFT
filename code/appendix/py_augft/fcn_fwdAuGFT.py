import numpy as np

def fcn_fwdAuGFT(x, U, Q, alpha=(1-1/np.sqrt(2)) ):

    # Forward AuGFT
    I = np.eye(U.shape[0])
    beta = np.sqrt(alpha * (2 - alpha))
    fwdAuGFT = np.hstack(((I - alpha * (Q@Q.T)) @ U, beta * Q)).T
    X = fwdAuGFT@x

    return X, fwdAuGFT

import numpy as np
from gnuradio import gr
import math
import scipy

class blk(gr.sync_block):
    def __init__(self):
        gr.sync_block.__init__(
            self,
            name="SOBI Block",
            in_sig=[(np.float32), (np.float32)],  # Assuming input data matrix size [m, N] = [3, 1000]
            out_sig=[(np.float32), (np.float32)]  # Only output the estimated source activities
        )

    def work(self, input_items, output_items):

        ## read the input
        X = np.vstack((input_items[0], input_items[1]))
        m=2
        N=8192
        print(X.shape)

        ## random matrix
        # Generate a random mxm matrix with uniform distribution
        np.random.seed(50)
        random_matrix = np.random.rand(m, m)
        # Normalize each row so that the sum of each row equals 1
        row_sums = random_matrix.sum(axis=1)
        A = random_matrix / row_sums


        ## multiplying matrix
        X_mixed = np.dot( A, X )

        ## centering the data
        X_centered= X_mixed - np.mean(X_mixed)

        ## whitening data
        covX = np.cov(X_centered)
        U, S, VV = scipy.linalg.svd(covX)
        S_inv = np.diag(1.0 / np.sqrt(S))
        whitening_matrix = U @ S_inv @ U.T
        X_whitened = whitening_matrix @ X_centered

        ## SOBI algorithm
        k = 1
        p=min(100,math.ceil(N/3))
        pm = p * m
        M = np.zeros((m, pm))
        Rxp = np.zeros((m, m))

        for u in range(0, pm, m):
            k += 1
            Rxp = np.dot(X_whitened[:, k-1:N], X_whitened[:, 0:N-k+1].T) / (N-k+1)
            M[:, u:u+m] = np.linalg.norm(Rxp, 'fro') * Rxp

        epsil = 1 / np.sqrt(N) / 100
        encore = True
        V = np.eye(m)

        while encore:
            encore = False
            for p in range(1,m-1):
                for q in range(p+1, m):
                    g = np.array([
                     M[p, p:m:pm] - M[q, q:m:pm],
                     M[p, q:m:pm] + M[q, p:m:pm],
                    1j * (M[q, p:m:pm] - M[p, q:m:pm])
                    ])
                    ggT = np.dot(g, g.T)
                    vcp, D = np.linalg.eig(np.real(ggT))
                    la = np.sort(np.diag(D))
                    K = np.argsort(np.diag(D))
                    angles = vcp[:, K[2]]
                    angles = np.sign(angles[0]) * angles
                    c = np.sqrt(0.5 + angles[0] / 2)
                    sr = 0.5 * (angles[1] - 1j * angles[2]) / c
                    sc = np.conj(sr)
                    oui = np.abs(sr) > epsil
                    encore = np.logical_or(encore, oui)

                    if oui:
                        colp = M[:, p:m:pm]
                        colq = M[:, q:m:pm]
                        M[:, p:m:pm] = c * colp + sr * colq
                        M[:, q:m:pm] = c * colq - sc * colp
                        rowp = M[p, :]
                        rowq = M[q, :]
                        M[p, :] = c * rowp + sc * rowq
                        M[q, :] = c * rowq - sr * rowp
                        temp = V[:, p]
                        V[:, p] = c * V[:, p] + sr * V[:, q]
                        V[:, q] = c * V[:, q] - sc * temp


        whitening_matrix_pinv = np.linalg.pinv(whitening_matrix)
        H= np.dot(whitening_matrix_pinv, V)
        S = np.dot(V.T, X_whitened)
        output_items[0][:] = S[0, :]
        output_items[1][:] = S[1, :]

        return len(output_items[0])

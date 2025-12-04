"""
Embedded Python Blocks:

Each time this file is saved, GRC will instantiate the first class it finds
to get ports and parameters of your block. The arguments to __init__  will
be the parameters. All of them are required to have default values!
"""

import numpy as np
from gnuradio import gr

class blk(gr.sync_block):
    def __init__(self, n_sources=None, n_delays=None):
        gr.sync_block.__init__(
            self,
            name="SOBI Block",
            in_sig=[(np.float32)],  # Assuming input data matrix size [m,N] = [3, 1000]
            out_sig=[(np.float32),(np.float32)]  # Only output the estimated source activities
        )
        self.n_sources = n_sources
        self.n_delays = n_delays

    def center_rows(self, X):
        return X - np.mean(X, axis=1, keepdims=True)

    def whiten_rows(self, X):
        covX = np.cov(X)
        D, E = np.linalg.eigh(covX)
        D_inv = np.diag(1.0 / np.sqrt(D))
        whitening_matrix = E @ D_inv @ E.T
        X_whitened = whitening_matrix @ X
        return X_whitened

    def work(self, input_items, output_items):
        X = input_items[0]  # Input data matrix [m, N]
        m, N = X.shape
        ntrials = 1

        if self.n_sources is None:
            n = m
        else:
            n = self.n_sources

        if self.n_delays is None:
            p = min(100, N // 3)
        else:
            p = self.n_delays

        X_centered = self.center_rows(X)
        X_whitened = self.whiten_rows(X_centered)

        k = 1
        pm = p * m
        M = np.zeros((m, pm))

        for u in range(0, pm, m):
            k += 1
            Rxp = np.zeros((m, m))
            for t in range(ntrials):
                Rxp += X_whitened[:, k:N] @ X_whitened[:, 0:N-k+1].T / (N-k+1) / ntrials

            M[:, u:u+m] = np.linalg.norm(Rxp, 'fro') * Rxp

        epsil = 1 / np.sqrt(N) / 100
        encore = True
        V = np.eye(m)

        while encore:
            encore = False
            for p in range(m-1):
                for q in range(p+1, m):
                    g = np.array([
                        M[p, p:m:pm] - M[q, q:m:pm],
                        M[p, q:m:pm] + M[q, p:m:pm],
                        1j * (M[q, p:m:pm] - M[p, q:m:pm])
                    ])
                    vcp, D = np.linalg.eig(np.real(g @ g.T))
                    angles = vcp[:, np.argmax(np.diag(D))]
                    angles = np.sign(angles[0]) * angles
                    c = np.sqrt(0.5 + angles[0] / 2)
                    sr = 0.5 * (angles[1] - 1j * angles[2]) / c
                    sc = np.conj(sr)
                    if abs(sr) > epsil:
                        encore = True
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

        S = V.T @ X_whitened

        output_items[0][:] = S[0, :]
        output_items[1][:] = S[1, :]
        return len(output_items[0])


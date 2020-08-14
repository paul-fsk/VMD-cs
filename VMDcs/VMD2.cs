//@author: Shengkun Fang
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;
using Numpy;
using Numpy.Models;

namespace VMDcs
{
    class VMD2
    {
        const int INT_MAX = Int32.MaxValue;
        const int INT_MIN = Int32.MinValue;

        const double eps_for_VMD = 2.2204e-16;
        static Slice all = new Slice(0, null);

        public static void Compute(ref NDarray u, ref NDarray u_hat, ref NDarray omega, NDarray<double> signal, double alpha, double tau, int K, int DC, int init, double tol)
        {
            var save_T = signal.len;
            var fs = 1 / (double)(save_T);

            var T = save_T;

            var f_mirror = np.zeros(2 * T);
            f_mirror[new Slice(0, T / 2)] = signal[new Slice( -T / 2 - 1,null ,-1)];

            f_mirror[new Slice(T/2,3*T/2)] = signal;
            
            f_mirror[new Slice(3 * T / 2, 2 * T)] = signal[new Slice(-1,-T/2-1,-1)]; 
            
            var f = f_mirror;

            T = f.len;

            var t = np.linspace(1 / (double)T, 1, T, true);

            var freqs = t - 0.5 - 1 / (double)T;

            var N = 500;
            
            var Alpha = alpha * np.ones(new Shape(K), np.complex64);

            var f_hat = np.fft.fftshift(np.fft.fft_(f));

            Complex[] ffh = f_hat.GetData<Complex>();
 
            var f_hat_plus = f_hat;
            f_hat_plus[new Slice(0, T / 2)]= (NDarray) 0;// = np.zeros(new Slice(0, T / 2));

            var u_hat_plus = np.zeros((N, freqs.len, K), np.complex64);

            var omega_plus = np.zeros((N, K), np.complex64);

            if (init == 1)
                for (int i = 1; i < K + 1; ++i)
                    omega_plus[0, i - 1] = (NDarray)((0.5 / K) * (i - 1));
            else if (init == 2)
                omega_plus[new Slice(0, null)] = np.sort(np.exp(np.log((NDarray)fs)) + (np.log((NDarray)0.5) - np.log((NDarray)fs)) * np.random.rand(1, K));
            else
                omega_plus[new Slice(0, null)] = (NDarray)0;

            if (DC!=0)
                omega_plus[(0, 0)] = (NDarray)0;

            var lamda_hat = np.zeros((N, freqs.len), np.complex64);

            NDarray uDiff = (NDarray)(tol + 2.2204e-16);

            int n = 1;
            NDarray sum_uk = np.zeros(new Shape(T),np.complex64);

            while (((double)uDiff>tol) && n < N)
            {
                int k = 1;
                var v = u_hat_plus[n - 1, all , K - 1];
                sum_uk = u_hat_plus[n - 1,all , K - 1] + sum_uk - u_hat_plus[n - 1, all, 0];
                var temp = f_hat_plus - sum_uk - lamda_hat[n - 1, all] / 2;
                u_hat_plus[n, all, k - 1] = (f_hat_plus - sum_uk - lamda_hat[n - 1, all] / 2) / (1 + Alpha[k - 1] * np.square(freqs - omega_plus[n - 1, k - 1]));

                if (DC == 0)
                    omega_plus[n, k - 1] = np.dot(freqs[new Slice(T / 2, T)], np.square(np.abs(u_hat_plus[n, new Slice(T / 2, T), k - 1])).T) / np.sum(np.square(np.abs(u_hat_plus[n, new Slice(T / 2, T), k - 1])));

                for(k = 2; k < K + 1; ++k)
                {
                    sum_uk = u_hat_plus[n,all, k - 2] + sum_uk - u_hat_plus[n - 1,all, k - 1];

                    u_hat_plus[n, all, k - 1] = (f_hat_plus - sum_uk - lamda_hat[n - 1, all] / 2) / (1 + Alpha[k - 1] * np.square(freqs - omega_plus[n - 1, k - 1]));

                    omega_plus[n, k - 1] = np.dot(freqs[new Slice(T / 2, T)], np.square(np.abs(u_hat_plus[n, new Slice(T / 2, T), k - 1])).T) / np.sum(np.square(np.abs(u_hat_plus[n, new Slice(T / 2, T), k - 1])));
                }

                t = np.sum(u_hat_plus[n, all, all], 1);

                lamda_hat[n,all] = lamda_hat[n - 1,all] + tau * (np.sum(u_hat_plus[n,all,all], 1) - f_hat_plus);

                n++;

                uDiff = (NDarray)2.2204e-16;

                for(int i = 1; i < K + 1; ++i)
                {
                    uDiff = uDiff + (1 / (double)T) * np.dot(u_hat_plus[n - 1, all, i - 1] - u_hat_plus[n - 2, all, i - 1], (np.conj(u_hat_plus[n - 1, all, i - 1] - u_hat_plus[n - 2, all, i - 1])).conj().T);
                }

                uDiff = np.abs(uDiff);
            }

            N = Math.Min(N, n);

            omega =  omega_plus[new Slice(0, N), all];

            u_hat = np.zeros((T, K), np.complex64);
            u_hat[new Slice(T / 2, T), all] = np.squeeze(u_hat_plus[N - 1, new Slice(T / 2, T), all]);

            u_hat[new Slice(T / 2, 0, -1), all] = np.squeeze(np.conj(u_hat_plus[N - 1, new Slice(T / 2, T), all]));

            u_hat[0, all] = np.conj(u_hat[-1, all]);

            u = np.zeros((K, t.len), np.complex64);

            for (int k = 1; k < K + 1; ++k)
            {
                u[k - 1, all] = np.real(np.fft.ifft(np.fft.ifftshift(u_hat[all, k - 1])));
            }
            
            u = u[all, new Slice(T / 4, 3 * T / 4)];

            u_hat = np.zeros((T / 2, K), np.complex64);

            for(int k = 1; k < K + 1; ++k)
            {
                u_hat[all, k - 1] = (np.fft.fft_(u[k - 1, all])).conj().T;
            }
            
        }

    }
}

//@author: Shengkun Fang
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Accord.Math;
using System.Numerics;
using Accord.Math.Transforms;

namespace VMDcs
{
    class VMD
    {
        const int INT_MAX = Int32.MaxValue;
        const int INT_MIN = Int32.MinValue;

        const double eps_for_VMD = 2.2204e-16;
        public static void Compute(ref double[,] u, ref Complex[,] u_hat, ref double[,] omega, double[] signal, double alpha, double tau, int K, int DC, int init, double tol)
        {
            int save_T = signal.Length;
            double fs = 1 / (double)save_T;

            //int T = save_T;

            double[] f_mirrow = new double[2 * save_T];

            Array.Copy(signal, 0, f_mirrow, (int)(save_T / 2), save_T);

            for (int i = 0; i < save_T / 2; ++i)
                f_mirrow[i] = signal[save_T / 2 - 1 - i];
            for (int i = 3 * save_T / 2; i < 2 * save_T; i++)
                f_mirrow[i] = signal[save_T + 3 * save_T / 2 - 1 - i];

            double[] f = f_mirrow;

            int T = f.Length;

            Complex[] freqs = Vector.Zeros<Complex>(T);
            double[] timevec = Vector.Zeros<double>(T);

            for (int i = 0; i < T; i++)
            {
                timevec[i] = (double)(i + 1.0) / T;
                freqs[i] = (timevec[i] - 0.5) - (double)(1 / T);
            }

            int N = 500;

            double[] Alpha = Vector.Create<double>(K, alpha);

            Complex[] freqvec = new Complex[T];

            for(int i = 0;i<freqvec.Length;i++)
            {
                freqvec[i] = f[i];
            }

            FourierTransform2.FFT(freqvec, FourierTransform.Direction.Forward);

            Complex[] f_hat =  circshift(freqvec,T/2);
            Complex[] f_hat_plus = Vector.Zeros<Complex>(f_hat.Length);

            Array.Copy(f_hat, T / 2, f_hat_plus, T / 2, T / 2);
            
            Complex[,,] u_hat_plus = new Complex[N, K, T];

            Complex[,] omega_plus = Matrix.Zeros<Complex>(N,K);

            double[] tmp;

            switch (init)
            {
                case 1:
                    for (int i = 0; i < K; i++)
                        omega_plus[0, i] = (double)(0.5 / K) * i;
                    break;
                case 2:
                    tmp = omega_init_method2(K, fs);
                    for (int i = 0; i < K; i++)
                        omega_plus[0, i] = tmp[i];
                    break;
                default:
                    break;
            }

            if (DC!=0)
                omega_plus[0, 0] = 0;

            Complex[,] lambda_hat = new Complex[N, T];

            double uDiff = tol + eps_for_VMD;

            int n = 1;// loop counter

            Complex[] sum_uk = new Complex[T];
            
            int k;

            while (uDiff > tol && n < N)
            {
                k = 1;

                for(int i = 0;i<sum_uk.Length;++i)
                {
                    sum_uk[i] = u_hat_plus[n - 1, K - 1,i] + sum_uk[i] - u_hat_plus[n - 1, 0,i];
                }

                //update spectrum of first mode through Wiener filter of residuals
                Complex[] Dividend_vec = new Complex[T];
                Complex[] Divisor_vec = new Complex[T];

                for (int i = 0; i < sum_uk.Length; ++i)
                {
                    Dividend_vec[i] = f_hat_plus[i] - sum_uk[i] - (lambda_hat[n - 1, i] / 2.0);
                    Divisor_vec[i] = (1 + Alpha[k - 1] * (freqs[i] - omega_plus[n - 1, k - 1]) * (freqs[i] - omega_plus[n - 1, k - 1]));
                    u_hat_plus[n, k - 1, i] = Dividend_vec[i] / Divisor_vec[i];
                }

                if (DC==0)
                {
                    Complex Dividend = new Complex(0,0), Divisor = new Complex(0, 0), Addend;
                    for (int i = 0; i < T - T / 2; i++)
                    {
                        Addend = u_hat_plus[n,k-1, T / 2 + i].Magnitude * u_hat_plus[n,k-1, T / 2 + i].Magnitude;
                        Divisor += Addend;
                        Dividend += freqs[T / 2 + i] * Addend;
                    }
                    omega_plus[n, k - 1] = Dividend / Divisor;

                }

                for (k = 1; k < K; k++)
                {
                    for (int i = 0; i < u_hat_plus.GetLength(2); i++)
                        sum_uk[i] = u_hat_plus[n, k - 1, i] + sum_uk[i] - u_hat_plus[n - 1, k, i];

                    Complex[] Dividend_vec1 = new Complex[T];
                    Complex[] Divisor_vec1 = new Complex[T];
                    
                    for (int i = 0; i < sum_uk.Length; ++i)
                    {
                        Dividend_vec1[i] = f_hat_plus[i] - sum_uk[i] - (lambda_hat[n - 1, i] / 2.0);
                        Divisor_vec1[i] = (1 + Alpha[k] * (freqs[i] - omega_plus[n - 1, k]) * (freqs[i] - omega_plus[n - 1, k]));
                        u_hat_plus[n, k, i] = Dividend_vec1[i] / Divisor_vec1[i];
                    }

                    Complex Dividend = new Complex(0, 0), Divisor = new Complex(0, 0), Addend;
                    for (int i = 0; i < T - T / 2; i++)
                    {
                        Addend = u_hat_plus[n, k, T / 2 + i].Magnitude * u_hat_plus[n, k, T / 2 + i].Magnitude;
                        Divisor += Addend;
                        Dividend += freqs[T / 2 + i] * Addend;
                    }
                    omega_plus[n, k] = Dividend / Divisor;

                }

                Complex[] sumv = sum(u_hat_plus, n);
                for(int i = 0; i < lambda_hat.GetLength(1); ++i)
                {
                    lambda_hat[n, i] = lambda_hat[n - 1, i] + tau * (sumv[i]-f_hat_plus[i]);
                }
                
                n++;

                Complex acc = new Complex(eps_for_VMD, 0);
                for(int i = 0; i < K; i++)
                {
                    Complex[] temp = new Complex[u_hat_plus.GetLength(2)];
                    Complex[] conjtemp = new Complex[u_hat_plus.GetLength(2)];
                    for (int j = 0; j < u_hat_plus.GetLength(2); j++)
                    {
                        temp[j] = u_hat_plus[n - 1, i, j] - u_hat_plus[n - 2, i, j];
                        conjtemp[j] =Complex.Conjugate(temp[j]);
                    }

                    var dottemp = Dot(temp, conjtemp);
                    acc = acc + dottemp / (double)T;
                }
                uDiff = acc.Magnitude;
            }

            N = Math.Min(N, n);

            omega = new double[N,omega_plus.GetLength(1)];

            for(int i = 0; i < omega_plus.GetLength(1); ++i)
            {
                for(int j = 0; j < N; j++)
                {
                    omega[j, i] = omega_plus[j, i].Real;
                }
            }

            u_hat = Matrix.Zeros<Complex>(T, K);

            for (int i = T / 2; i < T; i++)
                for (int j = 0; j < K; j++)
                {
                    //var v = u_hat_plus[N - 1, j, i]; 
                    u_hat[i, j] = u_hat_plus[N - 1,j,i];
                }

            for (int i = 0; i < T / 2; i++)
                for (int j = 0; j < K; j++)
                {
                    //u_hat[i, j] = Complex.Conjugate(u_hat_plus[N - 1, j, T - i - 1]);
                    u_hat[i, j] = Complex.Conjugate(u_hat_plus[N - 1, j, i]);
                }
            //for (int i = T / 2 - 1; i >= 0; i--)
            //    for (int j = 0; j < K; j++)
            //    {
            //        var v = Complex.Conjugate(u_hat_plus[N - 1, j, T - i - 1]);
            //        u_hat[i, j] = Complex.Conjugate(u_hat_plus[N - 1, j, T - i - 1]);
            //    }

            for (int i = 0;i< u_hat.GetLength(1); i++)
            {
                u_hat[0, i] = Complex.Conjugate(u_hat[u_hat.GetLength(1) - 1, i]); //Complex.Conjugate( u_hat[N - 1, i]);
            }

            u = Matrix.Zeros<double>(K,save_T);

            for(k = 0; k < K; k++)
            {
                Complex[] u_hat_col = ExtractColFromMatrixXcd(u_hat, k, T);
                u_hat_col = circshift(u_hat_col, (int)(Math.Floor((double)T / 2)));
                FourierTransform2.FFT(u_hat_col, FourierTransform.Direction.Backward);

                for(int t = 0; t < save_T; t++)
                    u[k, t] = u_hat_col[t + T / 4].Real;

            }

            //u_hat = Matrix.Zeros<Complex>(T, K);

            //double[] result_timevec = Vector.Zeros<double>(save_T);

            //for (int i = 0; i < save_T; i += 1)
            //{
            //    result_timevec[i] = (double)(i + 1) / save_T;
            //}

            //for (k = 0; k < K; k++)
            //{
            //    Complex[] u_row = ExtractRowFromMatrixXd(u, k, save_T);
            //    FourierTransform2.FFT(u_row,FourierTransform.Direction.Backward);
            //    u_row = circshift(u_row, save_T / 2);
            //    for (int t = 0; t < save_T; t++)
            //        u[k, t] = u_row[t].Real;
            //}
        }

        static Complex[] circshift(Complex[] data, int offset)
        {
            int n = data.Length;
            offset = offset % n;
           // Complex[] out_data = new Complex[n];
            if (offset == 0)
            {
                
                //Array.Copy(data,out_data,n);
                return data;
            }
            if (offset < 0)//向前移动-offset个元素 =向后移动(n-offset)个元素
                offset = n + offset;
            Complex[] out_data = Vector.Create<Complex>(data);
            reverse<Complex>(ref out_data, 0, n - 1);
            reverse<Complex>(ref out_data, 0, n - offset - 1);
            reverse<Complex>(ref out_data, n - offset, n - 1);
            return out_data;
        }

        static void reverse<T>(ref T[] v, int s, int l)
        {
            int n = v.Length;
            var temp = v[0];
            for (int i = s, j = l; i < j; i++, j--)
            {
                temp = v[i];
                v[i] = v[j];
                v[j] = temp;
            }
            return;
        }

        static double[] omega_init_method2(int K, double fs)
        {
            double[] res = Vector.Create<double>(K, 0);
            int N = INT_MAX / 2;
            Random random = new Random();

            for (int i = 0; i < K; i++)
            {
                res[i] = Math.Exp(Math.Log(fs) + (Math.Log(0.5) - Math.Log(fs)) *
                    (random.Next(N) / (double)(N + 1))
                );
            }

            Array.Sort(res);
            return res;
        }

        static Complex[,] vector_to_MatrixXcd_in_col(ref Complex[] Input)
        {
            int T = Input.Length;
            Complex[,] tmp= new Complex[1, T];
            for (int t = 0; t < T; t++)
                tmp[0, t] = Input[t];
            return tmp;
        }

        static Complex[] sum(Complex[,,] u_hat_plus, int n)
        {
            Complex[] cov = new Complex[u_hat_plus.GetLength(2)];
            for(int i = 0; i < cov.Length; i++)
            {
                Complex sum = 0;
                for (int j = 0; j < u_hat_plus.GetLength(1); ++j)
                    sum += u_hat_plus[n, j, i];
                cov[i] = sum;
            }

            return cov;
        }

        static Complex Dot(Complex[] a, Complex[] b)
        {
            Complex sum = 0;

            for(int i = 0; i < a.Length; ++i)
            {
                sum += a[i] * b[i];
            }

            return sum;
        }

        static Complex[] ExtractColFromMatrixXcd(Complex[,] Input, int ColIdx, int RowNum)
        {
            Complex[] Output = Vector.Zeros<Complex>(RowNum);
            for (int i = 0; i < RowNum; ++i)
                Output[i] = Input[i, ColIdx];
            return Output;
        }

        static Complex[] ExtractRowFromMatrixXd(double[,] Input, int RowIdx, int ColNum)
        {
            Complex[] Output = Vector.Zeros<Complex>(ColNum);
            for (int i = 0; i < ColNum; ++i)
                Output[i] = Input[RowIdx, i];
            return Output;
        }
    }
}

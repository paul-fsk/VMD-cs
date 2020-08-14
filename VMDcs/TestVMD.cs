//@author: Shengkun Fang
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using VMDcs.Forms;
using System.Windows.Forms.DataVisualization.Charting;
using System.Numerics;
using System.Diagnostics;
using System.Drawing;
using Numpy;

namespace VMDcs
{
    class TestVMD
    {
        [STAThread]
        static void Main()
        {
            //TestSpline();

            DataVisuForm form = new DataVisuForm();

            List<Chart> chs = new List<Chart>();

            Test(ref chs);

            form.setChart(chs);

            form.ShowDialog();



        }

        private static void Test(ref List<Chart> ch)  ///test code
        {
            double f_1 = 2.0, f_2 = 24.0, f_3 = 288.0;
            int T = 1000;

            double[] t = new double[T], v_1 = new double[T], v_2 = new double[T], v_3 = new double[T], signal = new double[T];

            for (int i = 0; i < T; i++)
            {
                t[i] = (double)(i + 1) / T;
                v_1[i] = Math.Cos(2 * Math.PI * f_1 * t[i]);
                v_2[i] = Math.Cos(2 * Math.PI * f_2 * t[i]) / 4.0;
                v_3[i] = Math.Cos(2 * Math.PI * f_3 * t[i]) / 16.0;
                signal[i] = v_1[i] + v_2[i] + v_3[i];
            }

            double alpha = 2000.0, tau = 0, tol = 1e-7;
            int K = 5, DC = 0, init = 1;

            double[,] u = new double[1, 1], omega = new double[1, 1];
            Complex[,] u_hat = new Complex[1, 1];

            var signal1 = np.array<double>(signal);

            NDarray u1 = np.array(u);
            NDarray u_hat1 = np.array(u_hat);
            NDarray omega1 = np.array(omega);
            VMD2.Compute(ref u1, ref u_hat1, ref omega1, signal1, alpha, tau, K, DC, init, tol);


            List<double[]> output = new List<double[]>();
            for (int i = 0; i < u1.shape[0]; ++i)
            {
                double[] ele = new double[u1.shape[1]];
                for (int j = 0; j < u1.shape[1]; ++j)
                {
                    if (u1.dtype.ToString().Equals("complex64"))
                        ele[j] = (double)(u1[i, j].real);
                    else if (u1.dtype.ToString().Equals("float64"))
                        ele[j] = (double)(u1[i, j]);


                }
                output.Add(ele);
            }


            //VMD.Compute(ref u, ref u_hat, ref omega, signal, alpha, tau, K, DC, init, tol);

            //List<double[]> output = new List<double[]>();
            //for (int i = 0; i < u.GetLength(0); ++i)
            //{
            //    double[] ele = new double[u.GetLength(1)];
            //    for (int j = 0; j < u.GetLength(1); ++j)
            //    {
            //        ele[j] = u[i, j];
            //    }
            //    output.Add(ele);
            //}

            PlotEMD(ref ch, "VMD test", t, signal, output);
        }

        private static void PlotEMD(ref List<Chart> chs, string title, double[] x, double[] z, List<double[]> emd)
        {
            chs = new List<Chart>();

            var chart = new Chart();
            chart.Size = new Size(700, 180);
            chart.Titles.Add(title);
            chart.Legends.Add(new Legend("Legend"));

            ChartArea ca = new ChartArea("DefaultChartArea");
            ca.AxisX.Title = "X";
            ca.AxisY.Title = "Y";
            chart.ChartAreas.Add(ca);

            // Series s1 = CreateSeries(chart, "Spline", CreateDataPoints(xs, ys), Color.Blue, MarkerStyle.None);
            Series s = CreateSeries(chart, "Original", CreateDataPoints(x, z), Color.Green, MarkerStyle.None);

            chart.Series.Add(s);
            chs.Add(chart);
            //chart.Series.Add(s1);
            int i = 1;
            foreach (var imf in emd)
            {
                var chart1 = new Chart();
                chart1.Size = new Size(700, 180);
                chart1.Titles.Add(title);
                String text = "imf" + i.ToString();
                if (i == emd.Count)
                    text = "res";
                chart1.Legends.Add(new Legend(text));

                ChartArea ca1 = new ChartArea("ChartArea");
                ca1.AxisX.Title = "X";
                ca1.AxisY.Title = "Y";
                chart1.ChartAreas.Add(ca1);

                Series s1 = CreateSeries(chart1, text, CreateDataPoints(x, imf), Color.Blue, MarkerStyle.None);

                chart1.Series.Add(s1);

                chs.Add(chart1);
                ++i;
            }

        }

        private static Series CreateSeries(Chart chart, string seriesName, IEnumerable<DataPoint> points, Color color, MarkerStyle markerStyle = MarkerStyle.None)
        {
            var s = new Series()
            {
                XValueType = ChartValueType.Double,
                YValueType = ChartValueType.Double,
                Legend = chart.Legends[0].Name,
                IsVisibleInLegend = true,
                ChartType = SeriesChartType.Line,
                Name = seriesName,
                ChartArea = chart.ChartAreas[0].Name,
                MarkerStyle = markerStyle,
                Color = color,
                MarkerSize = 8
            };

            foreach (var p in points)
            {
                s.Points.Add(p);
            }

            return s;
        }

        private static List<DataPoint> CreateDataPoints(double[] x, double[] y)
        {
            Debug.Assert(x.Length == y.Length);
            List<DataPoint> points = new List<DataPoint>();

            for (int i = 0; i < x.Length; i++)
            {
                points.Add(new DataPoint(x[i], y[i]));
            }

            return points;
        }
    }
}

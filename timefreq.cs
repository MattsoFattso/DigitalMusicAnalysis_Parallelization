using System;
using System.Diagnostics;
using System.Numerics;
using System.Threading.Tasks;
using System.IO;

namespace DigitalMusicAnalysis
{
    public class timefreq
    {
        public float[][] timeFreqData;
        public int wSamp;
        public Complex[] twiddles;

        public timefreq(float[] x, int windowSamp)
        {
            int ii;
            double pi = 3.14159265;
            Complex i = Complex.ImaginaryOne;
            this.wSamp = windowSamp;
            twiddles = new Complex[wSamp];
            for (ii = 0; ii < wSamp; ii++)
            {
                double a = 2 * pi * ii / (double)wSamp;
                twiddles[ii] = Complex.Pow(Complex.Exp(-i), (float)a);
            }

            timeFreqData = new float[wSamp/2][];

            int nearest = (int)Math.Ceiling((double)x.Length / (double)wSamp);
            nearest = nearest * wSamp;

            Complex[] compX = new Complex[nearest];
            for (int kk = 0; kk < nearest; kk++)
            {
                if (kk < x.Length)
                {
                    compX[kk] = x[kk];
                }
                else
                {
                    compX[kk] = Complex.Zero;
                }
            }


            int cols = 2 * nearest /wSamp;

            for (int jj = 0; jj < wSamp / 2; jj++)
            {
                timeFreqData[jj] = new float[cols];
            }

            Stopwatch stftTimer = Stopwatch.StartNew();

            timeFreqData = stft(compX, wSamp);

            stftTimer.Stop();

            File.AppendAllText(
                "stft_parallel_4_core.txt",
                $"{stftTimer.Elapsed.TotalMilliseconds:F3} ms{Environment.NewLine}");
            double elapsedMs = stftTimer.Elapsed.TotalMilliseconds;
            Console.WriteLine($"STFT Time: {elapsedMs:F3} ms");


        }
        float[][] stft(Complex[] x, int wSamp)
        {
            int ii = 0;
            int jj = 0;
            int kk = 0;
            int ll = 0;
            int N = x.Length;
            float fftMax = 0;

            float[][] Y = new float[wSamp / 2][];

            for (ll = 0; ll < wSamp / 2; ll++)
            {
                Y[ll] = new float[2 * (int)Math.Floor((double)N / (double)wSamp)];
            }

            Complex[] temp = new Complex[wSamp];
            Complex[] tempFFT = new Complex[wSamp];

            for (ii = 0; ii < 2 * Math.Floor((double)N / (double)wSamp) - 1; ii++)
            {

                for (jj = 0; jj < wSamp; jj++)
                {
                    temp[jj] = x[ii * (wSamp / 2) + jj];
                }

                tempFFT = fft(temp);

                for (kk = 0; kk < wSamp / 2; kk++)
                {
                    Y[kk][ii] = (float)Complex.Abs(tempFFT[kk]);

                    if (Y[kk][ii] > fftMax)
                    {
                        fftMax = Y[kk][ii];
                    }
                }


            }

            for (ii = 0; ii < 2 * Math.Floor((double)N / (double)wSamp) - 1; ii++)
            {
                for (kk = 0; kk < wSamp / 2; kk++)
                {
                    Y[kk][ii] /= fftMax;
                }
            }

            return Y;
        }
        float[][] parallell_stft(Complex[] x, int wSamp)
        {
            int ll = 0;
            int N = x.Length;           
            float[][] Y = new float[wSamp / 2][];

            for (ll = 0; ll < wSamp / 2; ll++)
            {
                Y[ll] = new float[2 * (int)Math.Floor((double)N / (double)wSamp)];
            }

            int numberOfWindows =
              2 * (int)Math.Floor((double)N / (double)wSamp) - 1;

            ParallelOptions options = new ParallelOptions
            {
                MaxDegreeOfParallelism = 4
            };


            // STEP 1: Process each STFT window in parallel
            Parallel.For(0, numberOfWindows, options, ii =>
            {
                // Each worker gets its own arrays
                Complex[] temp = new Complex[wSamp];

                for (int jj = 0; jj < wSamp; jj++)
                {
                    temp[jj] = x[ii * (wSamp / 2) + jj];
                }

                Complex[] tempFFT = fft(temp);

                // Each worker writes to its own column of Y
                for (int kk = 0; kk < wSamp / 2; kk++)
                {
                    Y[kk][ii] = (float)Complex.Abs(tempFFT[kk]);
                }
            });


            // STEP 2: Find maximum AFTER parallel processing
            float fftMax = 0;

            for (int ii = 0; ii < numberOfWindows; ii++)
            {
                for (int kk = 0; kk < wSamp / 2; kk++)
                {
                    if (Y[kk][ii] > fftMax)
                    {
                        fftMax = Y[kk][ii];
                    }
                }
            }


            for (int ii = 0; ii < 2 * Math.Floor((double)N / (double)wSamp) - 1; ii++)
            {
                for (int kk = 0; kk < wSamp / 2; kk++)
                {
                    Y[kk][ii] /= fftMax;
                }
            }

            return Y;
        }

        Complex[] fft(Complex[] x)
        {
            int ii = 0;
            int kk = 0;
            int N = x.Length;

            Complex[] Y = new Complex[N];

            // NEED TO MEMSET TO ZERO?

            if (N == 1)
            {
                Y[0] = x[0];
            }
            else{

                Complex[] E = new Complex[N/2];
                Complex[] O = new Complex[N/2];
                Complex[] even = new Complex[N/2];
                Complex[] odd = new Complex[N/2];

                for (ii = 0; ii < N; ii++)
                {

                    if (ii % 2 == 0)
                    {
                        even[ii / 2] = x[ii];
                    }
                    if (ii % 2 == 1)
                    {
                        odd[(ii - 1) / 2] = x[ii];
                    }
                }

                E = fft(even);
                O = fft(odd);

                for (kk = 0; kk < N; kk++)
                {
                    Y[kk] = E[(kk % (N / 2))] + O[(kk % (N / 2))] * twiddles[kk * wSamp / N];
                }
            }

           return Y;
        }
        
    }
}

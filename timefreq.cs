using System;
using System.Numerics;

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

            timeFreqData = stft(compX, wSamp);
	
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

                tempFFT = fft_Iterative(temp);

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

        // Iterative radix-2 FFT for pitch detection
        private Complex[] fft_Iterative(Complex[] x)
        {
            int N = x.Length;

            // FFT requires a power-of-two input size
            if (N == 0 || (N & (N - 1)) != 0)
            {
                throw new ArgumentException(
                    "FFT input length must be a power of two.");
            }

            // Copy input so original data is not modified
            Complex[] data = new Complex[N];

            int bits = (int)Math.Log(N, 2);

            // Bit-reversal permutation
            for (int i = 0; i < N; i++)
            {
                int reversed = ReverseBits(i, bits);
                data[reversed] = x[i];
            }

            // Process FFT one stage at a time
            for (int size = 2; size <= N; size *= 2)
            {
                int halfSize = size / 2;

                // Process each group within this stage
                for (int start = 0; start < N; start += size)
                {
                    for (int j = 0; j < halfSize; j++)
                    {
                        double angle =
                            -2.0 * Math.PI * j / size;

                        Complex twiddle = new Complex(
                            Math.Cos(angle),
                            Math.Sin(angle));

                        int evenIndex = start + j;
                        int oddIndex = evenIndex + halfSize;

                        Complex even = data[evenIndex];
                        Complex odd =
                            data[oddIndex] * twiddle;

                        // Butterfly operation
                        data[evenIndex] = even + odd;
                        data[oddIndex] = even - odd;
                    }
                }
            }

            return data;
        }

        private int ReverseBits(int value, int bits)
        {
            int reversed = 0;

            for (int i = 0; i < bits; i++)
            {
                reversed =
                    (reversed << 1) | (value & 1);

                value >>= 1;
            }

            return reversed;
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

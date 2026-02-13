// See https://aka.ms/new-console-template for more information

using Benchmark;
using BenchmarkDotNet.Running;
using EngineeringFluids;
using System.Numerics;


float[] floatArray1 = Enumerable.Range(0, 1001).Select(i => (float)i).ToArray();
float[] floatArray2 = Enumerable.Range(0, 1001).Select(i => (float)i).ToArray();



//var test = SumVector(floatArray1, floatArray2);
//var test2 = Sum(floatArray1, floatArray2);
var test3 = SumSpan(floatArray1, floatArray2);

//Speed();

BenchmarkRunner.Run<Benchy>();



Console.WriteLine("Hello, World!");





void Speed()
{

    double T = 400.0; // Temperature in K
    double rho = 9.0; // Density in kg/m^3

    double Tc = 405.56; // Critical temperature in K
    double rhoc = 13696.0; // Critical density in mol/m^3
    double M = 0.01703056; // Molar mass in kg/mol
    double R = 8.3144598; // Gas constant in J/(mol*K)

    // Calculate reduced temperature and density
    double tau = Tc / T;
    double delta = rho / (rhoc * M);

    float tauF = (float)(Tc / T);
    float deltaF = (float)(rho / (rhoc * M));




    //for (int i = 0; i < 20000; i++)
    //{
    //   var Alpha0 = Ammonia.alpha0(delta, tau);

    //    var AlphaR = Ammonia.alphaR(delta, tau);
    //}




}


float[] SumVector(float[] left, float[] right)
{
    if (left is null)
    {
        throw new ArgumentNullException(nameof(left));
    }

    if (right is null)
    {
        throw new ArgumentNullException(nameof(right));
    }

    if (left.Length != right.Length)
    {
        throw new ArgumentException($"{nameof(left)} and {nameof(right)} are not the same length");
    }

    int length = left.Length;
    float[] result = new float[length];

    // Get the number of elements that can't be processed in the vector
    // NOTE: Vector<T>.Count is a JIT time constant and will get optimized accordingly
    int remaining = length % Vector<float>.Count;

    for (int i = 0; i < length - remaining; i += Vector<float>.Count)
    {
        var v1 = new Vector<float>(left, i);
        var v2 = new Vector<float>(right, i);
        (v1 * v2).CopyTo(result, i);
    }

    for (int i = length - remaining; i < length; i++)
    {
        result[i] = left[i] + right[i];
    }

    return result;
}

float[] Sum(float[] left, float[] right)
{
    if (left is null)
    {
        throw new ArgumentNullException(nameof(left));
    }

    if (right is null)
    {
        throw new ArgumentNullException(nameof(right));
    }

    if (left.Length != right.Length)
    {
        throw new ArgumentException($"{nameof(left)} and {nameof(right)} are not the same length");
    }

    int length = left.Length;
    float[] result = new float[length];


    for (int i = 0; i < length; i++)
    {
        result[i] = left[i] * right[i];
    }

    return result;
}

Span<float> SumSpan(Span<float> left, Span<float> right)
{
    int length = left.Length;
    Span<float> result = new float[length];


    for (int i = 0; i < length; i++)
    {
        result[i] = left[i] * right[i];
    }

    return result;
}
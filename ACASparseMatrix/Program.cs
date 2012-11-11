using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MathNet.Numerics.LinearAlgebra.Double;

namespace ACASparseMatrix
{
    class Program
    {
        static void Main(string[] args)
        {
            DenseVector d = new DenseVector(20);

            foreach (double t in d)
            {
                Console.WriteLine(t);
            }
        }
    }
}

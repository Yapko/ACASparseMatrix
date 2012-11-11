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
            //ololo
            bool val = true;
            if (val)
            {
                Console.Write("1");
                val = false;
            }
            foreach (double t in d)
            {
                Console.WriteLine(t);
            }
        }
    }
}

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
            BasicFuncBoxes b = new BasicFuncBoxes(10, 10);
            Vector v = new DenseVector(10);
            for (int i = 0; i < v.Count; i++)
            {
                v[i] = -1 + i;
            }

            b.X.SetColumn(0, v);
        }
    }
}

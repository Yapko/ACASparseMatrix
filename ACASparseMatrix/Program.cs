using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MathNet.Numerics.LinearAlgebra.Single;

namespace ACASparseMatrix
{
    class Program
    {
        static void Main(string[] args)
        {
            DenseVector d = new DenseVector(20);
<<<<<<< HEAD

            
=======
            //ololo
            bool val = true;
            if (val)
            {
                Console.Write("1");
                val = false;
            }
>>>>>>> adf64fb2da689cc93ddfba5aea3e2ced439fcaba
            foreach (double t in d)
            {
                Console.WriteLine(t);
            }

            
        }
    }
}

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MathNet.Numerics.LinearAlgebra.Double;
using System.IO;

namespace ACASparseMatrix
{
    class Program
    {
        public static void f(Vector v)
        {
            v[0] = 100500;
        }
        static void Main(string[] args)
        {
            
            BasicFuncBoxes b = new BasicFuncBoxes(10, 10);
            
            double finest_level_size = 0.25; // Size of finest level (lambdas) of multilevel subdivision
            double ACA_thres = 1e-3;       // Relative error threshold to stop adding rows and columns in ACA iteration
            double precon_radius = 0.15;   // Precondioner size (meters): elements of Z with R < precon_radius   
            double tol = 1e-2;             // Convergence tolerance in iterative solver
            double maxN_Z = 5000;

            StreamReader xRead = new StreamReader("X.txt");
            StreamReader yRead = new StreamReader("Y.txt");
            StreamReader zRead = new StreamReader("Z.txt");

            StreamReader jRead = new StreamReader("J.txt");

            Vector rcx;
            Vector rcy;
            Vector rcz;

            String str = xRead.ReadToEnd();

            List<double> lst = new List<double>();

            foreach (string s in str.Split('\t'))
            {
                lst.Add(double.Parse(s));
            }

            rcx = new DenseVector(lst.ToArray());

            str = yRead.ReadToEnd();
            lst.Clear();
            
            foreach (string s in str.Split('\t'))
            {
                lst.Add(double.Parse(s));
            }

            rcy = new DenseVector(lst.ToArray());

            str = zRead.ReadToEnd();
            lst.Clear();

            foreach (string s in str.Split('\t'))
            {
                lst.Add(double.Parse(s));
            }

            rcz = new DenseVector(lst.ToArray());

            int N = 3072;
            //Fill Z matrix
            if( N <= maxN_Z)
            {
               //Z = user_impedance(1:N, 1:N, OG_data,EM_data);
               List<int> n = new List<int>(N);
               for (int i = 0; i < N; i++)
               {
                   n.Add(i);
               }
               Matrix Z = ACA.UserImpedance(n,n);
               // Copy of uncompressed Z for iterative solver
               //Z_uncomp = cell(1,1);
               NewSparseMatrix Z_uncomp = new NewSparseMatrix();
               Z_uncomp.Add(new ACAStruct(n, n, Z, new DenseMatrix(1, 1), new DenseMatrix(1, 1), 0, 1));
            }
            //CHECKED - OK everything before works well
            BasicFuncBoxes bfb = ACA.PrepareMultilevel(rcx, rcy, rcz, N, finest_level_size);
            NewSparseMatrix Z_comp = new NewSparseMatrix();
            ACA.MultilevelCompres(bfb, 0, 0, 0, 0, 0, 0, 0, bfb.L, ACA_thres, ref Z_comp);            
                   
            
  
            //testing Multiply(matvec in mathlab)
            Vector J = new DenseVector(3072);
            int k = 0;
            while (jRead.EndOfStream == false)
            {
                string s = jRead.ReadLine();
                string[] a = s.Split(' ');
                int p = 3;
                if (a.Length == 5) p = 2;
                J[k] = double.Parse(a[p].Replace('.',','));
                k++;
            }

            Vector r = Z_comp.Multiply(J, true);
        }
    }
}

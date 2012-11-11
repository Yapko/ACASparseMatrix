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
            foreach (double t in d)
            {
                Console.WriteLine(t);
            }

            
        }

      /*public Vector Multiply(Matrix Zcomp, Vector J, bool sym_source_field)
        {
            Vector y = new DenseVector(J.Count);
            int Nblocks = Zcomp.Count;
            for(int nb = 1; nb <= Nblocks; nb++)
            {
                Vector m = Zcomp[nb].GetM;
                Vector n = Zcomp[nb].GetN;
                if (Zcomp[nb].self == 1) // Self-interactions
                {
                    y(m) = y(m) + Zcomp{nb}.Z * J(n);
                }
                else //Non-self-interactions
                {
                    if (sym_source_field == true)// We have to take into account both interactions (direct and the symmetric one)
                    { 
                        if (Zcomp{nb}.comp == 0)
                        {
                            y(m) = y(m) + Zcomp{nb}.Z * J(n);
                            y(n) = y(n) + Zcomp{nb}.Z.' * J(m);
                        }
                        else
                        {
                            y(m) = y(m) + Zcomp{nb}.U * (Zcomp{nb}.V * J(n));
                            y(n) = y(n) + Zcomp{nb}.V.' * (Zcomp{nb}.U.' * J(m));
                        }
                    }
                    else //All the interactions have been computed
                    {
                        if (Zcomp{nb}.comp == 0)
                        {
                            y(m) = y(m) + Zcomp{nb}.Z * J(n);
                        }
                        else
                        {
                            y(m) = y(m) + Zcomp{nb}.U * (Zcomp{nb}.V * J(n));
                        }
                    }
                }
            }
        }*/
    }
}

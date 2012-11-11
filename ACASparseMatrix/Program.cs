using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MathNet.Numerics.LinearAlgebra.SingDoub
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

      /*public Vector Multiply(Matrix Zcomp, Vector NewSparseJ, bool sym_source_field)
        {
            Vector y = new DenseVector(J.Count);
            int Nblocks = Zcomp.Count;
            for(int nb = 1; nb <= Nblocks; nb++)
            {
                Vector m = Zcomp[nb].GetM;
List<int> m = Zcomp[nb].GetM;
                List<int> n = Zcomp[nb].GetN;
                if (Zcomp[nb].Self == 1) // Self-interactions
                {
                    //y(m) = y(m) + Zcomp[nb].Z * J(n);
                    Vector tempV = new DenseVector(n.Count);
                    foreach(int index in n)
                    {
                        tempV.Add(n[index]);
                    }
                    tempV = (Vector)Zcomp[nb].Z_Matrix.Multiply(tempV);
                    int tempCount = 0;
                    foreach(int index in m)
                    {
                        y[index]+= tempV[tempCount];
                        tempCount++;
                    }             else //Non-self-interactions
                {
                    if (sym_source_field == true)// We have to take into account both interactions (direct and the symmetric one)
                    { 
                        if (Zcomp{nb}.comp == 0)
                     [nb].Comp == 0)
                        {
                            //y(m) = y(m) + Zcomp[nb].Z_Matrix * J(n);
                            Vector tempV = new DenseVector(n.Count);
                            foreach (int index in n)
                            {
                                tempV.Add(n[index]);
                            }
                            tempV = (Vector)Zcomp[nb].Z_Matrix.Multiply(tempV);
                            int tempCount = 0;
                            foreach (int index in m)
                            {
                                y[index] += tempV[tempCount];
                                tempCount++;
                            }     
                            //y(n) = y(n) + Zcomp[nb].Z_Matrix.Transpose() * J(m);
                            tempV.Clear();
                            tempCount = 0;
                            tempV = new DenseVector(m.Count);
                            foreach (int index in m)
                            {
                                tempV.Add(m[index]);
                            }
                            tempV = (Vector)Zcomp[nb].Z_Matrix.Transpose().Multiply(tempV);
                            tempCount = 0;
                            foreach (int index in n)
                            {
                                y[index] += tempV[tempCount];
                                tempCount++;
                            } lse
                        {
                            y(m) = y(m) + Zcomp{nb}.U * (Zcomp{nb}.V * J(n));
      //y(m) = y(m) + Zcomp[nb].U_Vector * (Zcomp[nb].V_Vector * J(n));
                            Vector tempV = new DenseVector(n.Count);
                            foreach (int index in n)
                            {
                                tempV.Add(n[index]);
                            }
                            tempV = (Vector)Zcomp[nb].U_Vector.Multiply(Zcomp[nb].V_Vector.DotProduct(tempV));
                            int tempCount = 0;
                            foreach (int index in m)
                            {
                                y[index] += tempV[tempCount];
                                tempCount++;
                            }
                            //y(n) = y(n) + Zcomp[nb].V_Vector.T * (Zcomp[nb].U_Vector.T * J(m));
                            double dodPr = tempV.DotProduct(tempV    }
                    else //All the interactions have been computed
                    {
                        if (Zcomp{nb}.comp == 0)
                        {
        [nb].Comp == 0)             y(m) = y(m) + Zcomp{nb}.U * (Zcomp{nb}.V * J(n));
             [nb].Z_Matrix * J(n             else
                        {
                            y(m) = y(m) + Zcomp{nb}.U * (Zcomp{nb}.V * J(n));
             [nb].U_Vector * (Zcomp[nb].V_Vector}
                }
            }
        }*/
    }
}

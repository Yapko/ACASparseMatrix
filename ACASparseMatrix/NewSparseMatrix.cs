// -----------------------------------------------------------------------
// <copyright file="NewSparseMatrix.cs" company="">
// TODO: Update copyright text.
// </copyright>
// -----------------------------------------------------------------------

namespace ACASparseMatrix
{
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using System.Text;
    using MathNet.Numerics.LinearAlgebra.Double;
    /// <summary>
    /// TODO: Update summary.
    /// </summary>
    public class NewSparseMatrix:Matrix
    {
        List<ACAStruct> compressed;
        public NewSparseMatrix():base(1,1)
        {
            compressed = new List<ACAStruct>();
        }

        public override double At(int row, int column)
        {
            throw new NotImplementedException();
        }

        public override void At(int row, int column, double value)
        {
            throw new NotImplementedException();
        }

        public override MathNet.Numerics.LinearAlgebra.Generic.Vector<double> CreateVector(int size)
        {
            throw new NotImplementedException();
        }

        public override MathNet.Numerics.LinearAlgebra.Generic.Matrix<double> CreateMatrix(int numberOfRows, int numberOfColumns)
        {
            throw new NotImplementedException();
        }

        public  ACAStruct this[int index]
        {
            get
            {
                return compressed[index];
            }

            set
            {
                compressed[index] = value;
            }
       

       
        }

        public int Count
        {
            get
            {
                return compressed.Count;
            }
        }

        public Vector Multiply(NewSparseMatrix Zcomp, Vector J, bool sym_source_field)
        {
            Vector y = new DenseVector(J.Count);
            int Nblocks = Zcomp.Count;
            for (int nb = 1; nb <= Nblocks; nb++)
            {
                List<int> m = Zcomp[nb].GetM;
                List<int> n = Zcomp[nb].GetN;
                if (Zcomp[nb].Self == 1) // Self-interactions
                {
                    //y(m) = y(m) + Zcomp[nb].Z * J(n);
                    Vector tempV = new DenseVector(n.Count);
                    foreach (int index in n)
                    {
                        tempV.Add(J[index]);
                    }
                    tempV = (Vector)Zcomp[nb].Z_Matrix.Multiply(tempV);
                    int tempCount = 0;
                    foreach (int index in m)
                    {
                        y[index] += tempV[tempCount];
                        tempCount++;
                    }
                }
                else //Non-self-interactions
                {
                    if (sym_source_field == true)// We have to take into account both interactions (direct and the symmetric one)
                    {
                        if (Zcomp[nb].Comp == 0)
                        {
                            //y(m) = y(m) + Zcomp[nb].Z_Matrix * J(n);
                            Vector tempV = new DenseVector(n.Count);
                            foreach (int index in n)
                            {
                                tempV.Add(J[index]);
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
                                tempV.Add(J[index]);
                            }
                            tempV = (Vector)Zcomp[nb].Z_Matrix.Transpose().Multiply(tempV);
                            tempCount = 0;
                            foreach (int index in n)
                            {
                                y[index] += tempV[tempCount];
                                tempCount++;
                            }
                        }
                        else
                        {
                            //y(m) = y(m) + Zcomp[nb].U_Vector * (Zcomp[nb].V_Vector * J(n));
                            Vector tempV = new DenseVector(n.Count);
                            foreach (int index in n)
                            {
                                tempV.Add(J[index]);
                            }
                            tempV = (Vector)Zcomp[nb].U_Vector.Multiply(Zcomp[nb].V_Vector.DotProduct(tempV));
                            int tempCount = 0;
                            foreach (int index in m)
                            {
                                y[index] += tempV[tempCount];
                                tempCount++;
                            }
                            //y(n) = y(n) + Zcomp[nb].V_Vector.T * (Zcomp[nb].U_Vector.T * J(m));
                            Matrix tempMLeft = new DenseMatrix(1, Zcomp[nb].V_Vector.Count);
                            for (int i = 0; i < Zcomp[nb].V_Vector.Count; i++)
                            {
                                tempMLeft[0, i] = Zcomp[nb].V_Vector[i];
                            }

                            Matrix tempMRight = new DenseMatrix(m.Count, Zcomp[nb].U_Vector.Count);
                            for (int i = 0; i < m.Count; i++)
                            {
                                for (int j = 0; j < Zcomp[nb].U_Vector.Count; j++)
                                {
                                    tempMRight[i, j] = Zcomp[nb].U_Vector[i] * J[m[j]];
                                }
                            }

                            double tempS = 0;
                            for (int i = 0; i < Zcomp[nb].U_Vector.Count; i++)
                            {
                                for (int j = 0; j < Zcomp[nb].V_Vector.Count; j++)
                                {
                                    tempS += tempMLeft[0, i] * tempMRight[j, i];
                                }
                                y[n[i]] = tempS;
                                tempS = 0;
                            }
                        }
                    }
                    else //All the interactions have been computed
                    {
                        if (Zcomp[nb].Comp == 0)
                        {
                            //y(m) = y(m) + Zcomp[nb].Z_Matrix * J(n);
                            Vector tempV = new DenseVector(n.Count);
                            foreach (int index in n)
                            {
                                tempV.Add(J[index]);
                            }
                            tempV = (Vector)Zcomp[nb].Z_Matrix.Multiply(tempV);
                            int tempCount = 0;
                            foreach (int index in m)
                            {
                                y[index] += tempV[tempCount];
                                tempCount++;
                            }
                        }
                        else
                        {
                            //y(m) = y(m) + Zcomp[nb].U_Vector * (Zcomp[nb].V_Vector * J(n));
                            Vector tempV = new DenseVector(n.Count);
                            foreach (int index in n)
                            {
                                tempV.Add(J[index]);
                            }
                            tempV = (Vector)Zcomp[nb].U_Vector.Multiply(Zcomp[nb].V_Vector.DotProduct(tempV));
                            int tempCount = 0;
                            foreach (int index in m)
                            {
                                y[index] += tempV[tempCount];
                                tempCount++;
                            }
                        }
                    }
                }
            }
            return y;
        }
    }
}

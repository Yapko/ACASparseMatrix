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
    public class NewSparseMatrix: Matrix
    {
        List<ACAStruct> compressed;

        public NewSparseMatrix():base(0, 0)
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

        private ACAStruct this[int index]
        {
            get
            {
                return compressed[index];
            }
        }
        private int Count
        {
            get
            {
                return compressed.Count;
            }
        }

        /// <summary>
        /// adds new struct to matrix
        /// </summary>
        /// <param name="z">new struct</param>
        public void Add(ACAStruct z)
        {
            compressed.Add(z);
        }
        private Vector MultiplyByIndex(Matrix RightM, Vector LeftV, List<int> Index)
        {
            Vector leftVector = new DenseVector(0);
            foreach (int index in Index)
            {
                leftVector.Add(LeftV[index]);
            }
            leftVector = (Vector)RightM.Multiply(leftVector);
            return leftVector;
        }

        private Vector AddByIndex(Vector LeftV, List<int> LeftVIndex, Vector RighV)
        {
            Vector resultV = new DenseVector(RighV.Count);
            int count = 0;
            foreach (int index in LeftVIndex)
            {
                resultV[index] += RighV[count];
                count++;
            }
            return resultV;
        }

        public Vector Multiply(Vector J, bool sym_source_field)
        {
            NewSparseMatrix Zcomp = this;
            Vector y = new DenseVector(J.Count);
            int Nblocks = Zcomp.Count;
            for (int nb = 1; nb <= Nblocks; nb++)
            {
                List<int> m = Zcomp[nb].GetM;
                List<int> n = Zcomp[nb].GetN;
                if (Zcomp[nb].Self == 1) // Self-interactions
                {
                    //y(m) = y(m) + Zcomp[nb].Z * J(n);
                    AddByIndex(y, m, MultiplyByIndex(Zcomp[nb].Z_Matrix, J, n));
                }
                else //Non-self-interactions
                {
                    if (sym_source_field == true)// We have to take into account both interactions (direct and the symmetric one)
                    {
                        if (Zcomp[nb].Comp == 0)
                        {
                            //y(m) = y(m) + Zcomp[nb].Z_Matrix * J(n);
                            AddByIndex(y, m, MultiplyByIndex(Zcomp[nb].Z_Matrix, J, n));
                            //y(n) = y(n) + Zcomp[nb].Z_Matrix.Transpose() * J(m);
                            //if comlex use complex conjugate transpose of Z_Matrix elements
                            AddByIndex(y, n, MultiplyByIndex(Zcomp[nb].Z_Matrix, J, m));
                        }
                        else
                        {
                            
                            //y(m) = y(m) + Zcomp[nb].U_Vector * (Zcomp[nb].V_Vector * J(n));
                            AddByIndex(y, m, (Vector)Zcomp[nb].U_Vector.Multiply(MultiplyByIndex(Zcomp[nb].V_Vector, J, n)));
                            //y(n) = y(n) + Zcomp[nb].V_Vector.T * (Zcomp[nb].U_Vector.T * J(m));
                            //if comlex use complex conjugate transpose of V elements and U elements.
                            AddByIndex(y, n, (Vector)Zcomp[nb].V_Vector.Multiply(MultiplyByIndex(Zcomp[nb].U_Vector, J, m)));
                         }
                    }
                    else //All the interactions have been computed
                    {
                        if (Zcomp[nb].Comp == 0)
                        {
                            //y(m) = y(m) + Zcomp[nb].Z_Matrix * J(n);
                            AddByIndex(y, m, MultiplyByIndex(Zcomp[nb].Z_Matrix, J, n));
                        }
                        else
                        {
                            //y(m) = y(m) + Zcomp[nb].U_Vector * (Zcomp[nb].V_Vector * J(n));
                            AddByIndex(y, m, (Vector)Zcomp[nb].U_Vector.Multiply(MultiplyByIndex(Zcomp[nb].V_Vector, J, n)));
                        }
                    }
                }
            }
            return y;
        }
        
    }
}

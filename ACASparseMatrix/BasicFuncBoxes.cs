// -----------------------------------------------------------------------
// <copyright file="BasicFuncBoxes.cs" company="">
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
    public class BasicFuncBoxes
    {
        DenseMatrix x;
        DenseMatrix y;
        DenseMatrix z;

        /// <summary>
        /// number of rows of each matrix
        /// </summary>
        int n;
        /// <summary>
        /// number of cols of each matrix
        /// </summary>
        int l;

        /// <summary>
        /// ctor with params
        /// </summary>
        /// <param name="N">number of rows of each matrix</param>
        /// <param name="L">number of cols of each matrix</param>
        public BasicFuncBoxes(int N,int L)
        {
            n = N;
            l = L;
            x = new DenseMatrix(n, l);
            y = new DenseMatrix(n, l);
            z = new DenseMatrix(n, l);
        }

        /// <summary>
        /// property to get rows
        /// </summary>
        public int N
        {
            get
            {
                return n;
            }
        }

        /// <summary>
        /// property to get cols
        /// </summary>
        public int L
        {
            get
            {
                return l;
            }
        }

        /// <summary>
        /// property to get matrix x
        /// </summary>
        public DenseMatrix X
        {
            get
            {
                return x;
            }
        }

        /// <summary>
        /// property to get matrix y
        /// </summary>
        public DenseMatrix Y
        {
            get
            {
                return y;
            }
        }

        /// <summary>
        /// property to get matrix z
        /// </summary>
        public DenseMatrix Z
        {
            get
            {
                return z;
            }
        }
    }
}

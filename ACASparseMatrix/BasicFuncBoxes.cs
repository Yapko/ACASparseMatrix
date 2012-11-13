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
    /// Basic func boxes stores number of box that point belongs to at each level
    /// e.g. 
    /// x = [ 0 2    -- point 1
    ///       1 4    -- point 2
    ///       1 5    -- point 3
    ///       ]
    /// means that point 3 at first level belongs to box #1. At second level to box #5
    /// The same for y,z
    /// 
    /// </summary>
    public class BasicFuncBoxes
    {
        /// <summary>
        /// see class description
        /// </summary>
        DenseMatrix x;
        DenseMatrix y;
        DenseMatrix z;

        /// <summary>
        /// number of rows of each matrix, point count
        /// </summary>
        int n;
        /// <summary>
        /// number of cols of each matrix, tree height
        /// </summary>
        int l;

        /// <summary>
        /// ctor with params
        /// </summary>
        /// <param name="N">number of rows of each matrix, point count</param>
        /// <param name="L">number of cols of each matrix, it is height of the tree</param>
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

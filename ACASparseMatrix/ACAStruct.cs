// -----------------------------------------------------------------------
// <copyright file="ACAStruct.cs" company="">
// TODO: Update copyright text.
// </copyright>
// -----------------------------------------------------------------------

namespace ACASparseMatrix
{
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using System.Text;
    using MathNet.Numerics.LinearAlgebra.Generic;
    using MathNet.Numerics.LinearAlgebra.Double;
    /// <summary>
    /// 
    /// </summary>
    public class ACAStruct
    {
        #region Fields
        /// <summary>
        /// to save rows indexes
        /// </summary>
        List<int> m;
        //to save columns indexes
        List<int> n;
        /// <summary>
        /// Matrix is represented either by Z or by U*V
        /// </summary>
        Matrix Z;
        /// <summary>
        /// Matrix is represented either by Z or by U*V
        /// </summary>
        Matrix U;
        /// <summary>
        /// Matrix is represented either by Z or by U*V
        /// </summary>
        Matrix V;
        
        //
        int comp;
        //
        int self;

        //min and max in m
        int mMin, mMax;
        //min and max in n
        int nMin, nMax;
        #endregion

        #region Constructors
        /// <summary>
        /// default ctor
        /// </summary>
        public ACAStruct()
        {
            m = new List<int>();
            n = new List<int>();

            Z = null;

            U = null;
            V = null;

            comp = self = 0;

            mMin = mMax = 0;

            nMin = nMax = 0;
        }

        /// <summary>
        /// ctor with params
        /// </summary>
        /// <param name="M">new list of rows indexes</param>
        /// <param name="N">new list of columns indexes</param>
        /// <param name="Z1">new matrix of elements</param>
        /// <param name="U1"></param>
        /// <param name="V1"></param>
        /// <param name="Comp">new comp param</param>
        /// <param name="Self">new self param</param>
        public ACAStruct(
            List<int> M,
            List<int> N,
            Matrix Z1,
            Matrix U1,
            Matrix V1,
            int Comp,
            int Self            
            )
        {
            m = M;
            n = N;
            Z = Z1;
            U = U1;
            V = V1;

            
            mMax = m.Max();
            mMin = m.Min();

            nMax = n.Max();
            nMin = n.Min();
        }
        #endregion

        #region Properties
        /// <summary>
        /// 
        /// </summary>
        public int Comp
        {
            get
            {
                return comp;
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public int Self
        {
            get
            {
                return self;
            }
        }

       public  int MMax
        {
            get
            {
                return mMax;
            }
        }

        int MMin
        {
            get
            {
                return mMin;
            }
        }

        int NMax
        {
            get
            {
                return nMax;
            }
        }

        int NMin
        {
            get
            {
                return nMin;
            }
        }     
        
        public List<int> GetM
        {
            get
            {
                return m;
            }
        }

        public List<int> GetN
        {
            get
            {
                return n;
            }
        }

        public Matrix Z_Matrix
        {
            get
            {
                return Z;
            }
        }

        public Matrix U_Vector
        {
            get
            {
                return U;
            }
        }

        public Matrix V_Vector
        {
            get
            {
                return V;
            }
        }
        #endregion
    }
}

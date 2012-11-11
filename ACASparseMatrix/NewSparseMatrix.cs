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
    }
}

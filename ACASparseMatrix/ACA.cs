using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MathNet.Numerics.LinearAlgebra.Double;

namespace ACASparseMatrix
{
    class ACA
    {
        public static BasicFuncBoxes prepare_multilevel(Vector rcx, Vector rcy, Vector rcz, int N, double finest_level_size)
        {
            double xmax = rcx.Max(); double xmin = rcx.Min();
            double ymax = rcy.Max(); double ymin = rcy.Min();
            double zmax = rcz.Max(); double zmin = rcz.Min();

            int Lx = (int)Math.Ceiling(Math.Log((xmax - xmin) / (finest_level_size)) / Math.Log(2));
            int Ly = (int)Math.Ceiling(Math.Log((ymax - ymin) / (finest_level_size)) / Math.Log(2));
            int Lz = (int)Math.Ceiling(Math.Log((zmax - zmin) / (finest_level_size)) / Math.Log(2));
            int L = (int)Math.Max(Lx, Math.Max(Ly, Lz));

            double box_size_x = (xmax - xmin) * (1 + 1e-3);
            double box_size_y = (ymax - ymin) * (1 + 1e-3);
            double box_size_z = (zmax - zmin) * (1 + 1e-3);
            double box_size = Math.Max(box_size_x, Math.Max(box_size_y, box_size_z));

            BasicFuncBoxes basicFuncBoxes = new BasicFuncBoxes(N, L);

            for (int i = 0; i < L; i++)
            {
                box_size = box_size / 2;

                Vector tempVx = new DenseVector(rcx.Count);
                for (int j = 0; j < rcx.Count; j++)
                {
                    tempVx[j] = Math.Floor((rcx[j] - xmin) / box_size);
                }

                Vector tempVy = new DenseVector(rcy.Count);
                for (int j = 0; j < rcy.Count; j++)
                {
                    tempVy[j] = Math.Floor((rcy[j] - ymin) / box_size);
                }

                Vector tempVz = new DenseVector(rcz.Count);
                for (int j = 0; j < rcz.Count; j++)
                {
                    tempVz[j] = Math.Floor((rcz[j] - zmin) / box_size);
                }


                basicFuncBoxes.X.SetColumn(i, tempVx);
                basicFuncBoxes.Y.SetColumn(i, tempVy);
                basicFuncBoxes.Z.SetColumn(i, tempVz);
            }

            return basicFuncBoxes;
        }

        public static NewSparseMatrix multilevel_compress()
        {
            return new NewSparseMatrix();
        }

        /// <summary>
        /// Adaptive Cross Approximation (ACA) matrix compression
        /// the result is stored in U and V matrices like U*V
        /// </summary>
        /// <param name="acaThres">Relative error threshold to stop adding rows and columns in ACA iteration</param>
        /// <param name="m">Row indices of Z submatrix to compress</param>
        /// <param name="n">Column indices of Z submatrix to compress</param>
        /// <param name="U">to store result</param>
        /// <param name="V">to store result</param>
        public static void Aca(double acaThres, List<int> m, List<int> n, Matrix U, Matrix V)
        {
            int M = m.Count;
            int N = n.Count;

            //if Z is a vector, there is nothing to compress
            if (M == 1 || N == 1)
            {
                U = userImpedance(m, n);
                V = new DenseMatrix(1, 1);
                V[0, 0] = 1.0;
                return;
            }

            //Indices of columns picked up from Z
            Vector J = new DenseVector(N);
            //Indices of rows picked up from Z
            Vector I = new DenseVector(M);
            //Row indices to search for maximum in R 
            Vector i = new DenseVector(M);
            //Column indices to search for maximum in R
            Vector j = new DenseVector(N);

            for (int k = 2,t = 0; k <= M; k++)
            {
                i[t] = k;
                t++;
            }

            for (int k = 1, t = 0; k <= N; k++)
            {
                i[t] = k;
                t++;
            }

            //Initialization

            //Initialize the 1st row index I(1) = 1
            I[0] = 1;
        }

        /// <summary>
        /// Impedance submatrix elements
        /// </summary>
        /// <param name="m">Rows indexes of submatrix</param>
        /// <param name="n">Cols indexes of submatrix</param>
        /// <returns>Submatrix</returns>
        public static Matrix userImpedance(List<int> m, List<int> n)
        {
            int M = m.Count;
            int N = n.Count;

            Matrix Z = new DenseMatrix(M, N);

            for (int i = 1; i <= N; i++)
            {
                for (int j = 1; j <= M; j++)
                {
                    Z[i,j] = 1.0/(Math.Abs(m[j] - n[i] + 0.0001));
                }
            }

            return Z;
        }

    }
}

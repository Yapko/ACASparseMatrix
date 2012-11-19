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
                U = UserImpedance(m, n);
                V = new DenseMatrix(1, 1);
                V[0, 0] = 1.0;
                return;
            }

            //Indices of columns picked up from Z
            //Vector J = new DenseVector(N);
            List<int> J = new List<int>(N);
            //Indices of rows picked up from Z
            //Vector I = new DenseVector(M);
            List<int> I = new List<int>(M);
            //Row indices to search for maximum in R 
            //Vector i = new DenseVector(M);
            List<int> i = new List<int>(M);
            //Column indices to search for maximum in R
            //Vector j = new DenseVector(N);
            List<int> j = new List<int>(N);

            for (int k = 1, t = 0; k < M; k++)
            {
                i[t] = k;
                t++;
            }

            for (int k = 0, t = 0; k < N; k++)
            {
                i[t] = k;
                t++;
            }

            //Initialization

            //Initialize the 1st row index I(1) = 1
            I[0] = 0;

            //Initialize the 1st row of the approximate error matrix
            List<int> m0 = new List<int>();
            m0.Add(m[I[0]]);
            Matrix Rik = UserImpedance(m0, n);

            //Find the 1st column index J(0)
            double max = -1.0;
            int col = 0;

            foreach (int ind in j)
            {
                if (Math.Abs(Rik[0, ind]) > max)
                {
                    max = Math.Abs(Rik[0, ind]);
                    col = ind;
                }
            }

            J[0] = j[col];
            j.Remove(J[0]);

            //First row of V
            V.SetRow(0, Rik.Row(0).Divide(Rik[0, J[0]]));

            //Initialize the 1st column of the approximate error matrix
            List<int> n0 = new List<int>();
            n0.Add(n[J[0]]);
            Matrix Rjk = UserImpedance(m, n0);

            //First column of U
            U.SetColumn(0, Rjk.Column(0));

            // Norm of (approximate) Z, to test error
            double d1 = U.L2Norm();
            double d2 = V.L2Norm();
            double normZ = d1 * d1 * d2 * d2;

            //Find 2nd row index I(2)
            int row = 0;
            max = -1.0;

            foreach (int ind in i)
            {
                if (Math.Abs(Rjk[ind, 0]) > max)
                {
                    max = Math.Abs(Rjk[ind, 0]);
                    row = ind;
                }
            }

            I[1] = i[row];
            i.Remove(I[1]);

            //Iteration
            for (int k = 1; k < Math.Min(M, N); k++)
            {
                //Update (Ik)th row of the approximate error matrix:
                List<int> t1 = new List<int>();
                t1.Add(m[I[k]]);
                Rik = (Matrix)(UserImpedance(t1, n) - U.SubMatrix(I[k], U.ColumnCount, 0, 1).Multiply(V));

                //Find kth column index Jk
                max = -1.0;
                col = 0;

                foreach (int ind in j)
                {
                    if (Math.Abs(Rik[0, ind]) > max)
                    {
                        max = Math.Abs(Rik[0, ind]);
                        col = ind;
                    }
                }

                J[k] = j[col];
                j.Remove(J[k]);

                //Terminate if R(I(k),J(k)) == 0
                if (Rik[0, J[k]] == 0)
                {
                    break;
                }

                //Set k-th row of V equal to normalized error
                Matrix Vk = (Matrix)Rik.Divide(Rik[0, J[k]]);

                //Update (Jk)th column of the approximate error matrix
                List<int> n1 = new List<int>();
                n1.Add(J[k]);
                Rjk = (Matrix)(UserImpedance(m, n1) - U.Multiply(V.SubMatrix(0, 1, J[k], V.RowCount)));

                // Set k-th column of U equal to updated error
                Matrix Uk = Rjk;

                //Norm of approximate Z
                Matrix s = (Matrix)(U.Transpose().Multiply(Uk)).Multiply((Vk.Multiply(V.Transpose())).Transpose());
                double sum = 0;

                for (int i1 = 0; i1 < s.RowCount; i1++)
                {
                    for (int j1 = 0; j1 < s.ColumnCount; j1++)
                    {
                        sum += s[i1, j1];
                    }
                }

                d1 = Uk.L2Norm();
                d2 = Vk.L2Norm();

                normZ += 2 * sum + d1 * d1 * d2 * d2;

                //Update U and V
                U.InsertColumn(U.ColumnCount - 1, Uk.Column(0));
                V.InsertRow(V.RowCount - 1, Vk.Row(0));

                if (d1 * d2 <= acaThres * Math.Sqrt(normZ))
                {
                    break;
                }

                if (k == Math.Min(N, M) - 1)
                {
                    break;
                }

                max = -1;
                row = 0;

                foreach (int ind in i)
                {
                    if (Math.Abs(Rjk[0, ind]) > max)
                    {
                        max = Math.Abs(Rjk[0, ind]);
                        row = ind;
                    }
                }

                I[k + 1] = i[row];
                //i = removeIndex(i,I[k+1]);
                i.Remove(I[k + 1]);
            }
        }

        /// <summary>
        /// Impedance submatrix elements
        /// </summary>
        /// <param name="m">Rows indexes of submatrix</param>
        /// <param name="n">Cols indexes of submatrix</param>
        /// <returns>Submatrix</returns>
        public static Matrix UserImpedance(List<int> m, List<int> n)
        {
            int M = m.Count;
            int N = n.Count;

            Matrix Z = new DenseMatrix(M, N);

            for (int i = 1; i <= N; i++)
            {
                for (int j = 1; j <= M; j++)
                {
                    Z[i, j] = 1.0 / (Math.Abs(m[j] - n[i] + 0.0001));
                }
            }

            return Z;
        }

        #region MultilevelCompres
        public static NewSparseMatrix MultilevelCompres(BasicFuncBoxes basicFuncBoxes, double ix_s, double iy_s, double iz_s, double ix_f, double iy_f, double iz_f, int l, double L, double ACA_thres)
        {

            NewSparseMatrix Z_comp = new NewSparseMatrix();

            for (double xchs = 0; xchs <= 1; xchs++)
            {
                for (double ychs = 0; ychs <= 1; ychs++)
                {
                    for (double zchs = 0; zchs <= 1; zchs++)
                    {
                        // x-index of source child box at level l+1              
                        double ix_chs = ix_s * 2 + xchs;
                        // y-index of source child box at level l+1                        
                        double iy_chs = iy_s * 2 + ychs;
                        // z-index of source child box at level l+1
                        double iz_chs = iz_s * 2 + zchs;

                        //Find indices of basis functions in source child box
                        List<int> n = new List<int>(3);
                        for (int i = 0; i < basicFuncBoxes.X.RowCount; i++)
                        {
                            if(basicFuncBoxes.X[i, l + 1] == ix_chs && basicFuncBoxes.Y[i, l + 1] == iy_chs && basicFuncBoxes.Z[i, l + 1] == iz_chs)
                            {
                                n.Add(i);
                            }
                        }
                        if (n.Count == 0)
                        {
                            continue;
                        }

                        //Here we have a pair of non-empty source and field boxes
                       // if (Math.abs(ix_chs-ix_chf) > 1 || abs(iy_chs-iy_chf) > 1 || abs(iz_chs-iz_chf) > 1, % Far-field boxes
                          
            
                    }
                }
            }

            return new NewSparseMatrix();
        }
        #endregion

    }
}

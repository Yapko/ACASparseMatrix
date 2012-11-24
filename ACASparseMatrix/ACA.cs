using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MathNet.Numerics.LinearAlgebra.Double;

namespace ACASparseMatrix
{
    class ACA
    {
        public static BasicFuncBoxes PrepareMultilevel(Vector rcx, Vector rcy, Vector rcz, int N, double finest_level_size)
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

        /// <summary>
        /// Adaptive Cross Approximation (ACA) matrix compression
        /// the result is stored in U and V matrices like U*V
        /// </summary>
        /// <param name="acaThres">Relative error threshold to stop adding rows and columns in ACA iteration</param>
        /// <param name="m">Row indices of Z submatrix to compress</param>
        /// <param name="n">Column indices of Z submatrix to compress</param>
        /// <param name="U">to store result</param>
        /// <param name="V">to store result</param>
        public static Tuple<Matrix,Matrix> Aca(double acaThres, List<int> m, List<int> n, Matrix U, Matrix V)
        {
            int M = m.Count;
            int N = n.Count;
            int Min = Math.Min(M, N);
            U = new DenseMatrix(Min, Min);
            V = new DenseMatrix(Min, Min);
            //if Z is a vector, there is nothing to compress
            if (M == 1 || N == 1)
            {
                U = UserImpedance(m, n);
                V = new DenseMatrix(1, 1);
                V[0, 0] = 1.0;
                return new Tuple<Matrix,Matrix>(U,V);
            }

            //Indices of columns picked up from Z
            //Vector J = new DenseVector(N);
            //List<int> J = new List<int>(N);

            List<int> J = new List<int>(new int [N]);      
            //int[] J = new int[N];
            //Indices of rows picked up from Z
            //Vector I = new DenseVector(M);
            List<int> I = new List<int>(new int [M]);
            //int[] I = new int[M];
            //Row indices to search for maximum in R 
            //Vector i = new DenseVector(M);
            List<int> i = new List<int>();
            //int[] i = new int[M];
            //Column indices to search for maximum in R
            //Vector j = new DenseVector(N);
            List<int> j = new List<int>();
            //int[] j = new int[N];

            for (int k = 1; k < M; k++)
            {
                i.Add(k);
            }

            for (int k = 0; k < N; k++)
            {
                j.Add(k);
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

            //J[0] = j[col];
            J[0] = col;
            j.Remove(J[0]);

            //First row of V
            V = new DenseMatrix(1, Rik.ColumnCount);
            V.SetRow(0, Rik.Row(0).Divide(Rik[0, J[0]]));
            
            //Initialize the 1st column of the approximate error matrix
            List<int> n0 = new List<int>();
            n0.Add(n[J[0]]);
            Matrix Rjk = UserImpedance(m, n0);

            //First column of U
            U = new DenseMatrix(Rjk.RowCount, 1);
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

            //I[1] = i[row];
            I[1] = row;
            i.Remove(I[1]);

            //Iteration
            for (int k = 1; k < Math.Min(M, N); k++)
            {
                //Update (Ik)th row of the approximate error matrix:
                List<int> t1 = new List<int>();
                t1.Add(m[I[k]]);
                Rik = (Matrix)(UserImpedance(t1, n) - U.SubMatrix(I[k], 1, 0, U.ColumnCount).Multiply(V));
                //CHECKED OK all code before works fine
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

                J[k] = col;
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
                n1.Add(n[J[k]]);
                Rjk = (Matrix)(UserImpedance(m, n1) - U.Multiply(V.SubMatrix(0, V.RowCount, J[k], 1)));

                // Set k-th column of U equal to updated error
                Matrix Uk = Rjk;

                //Norm of approximate Z
                //Matrix s = (Matrix)(U.Transpose().Multiply(Uk)).Multiply((Vk.Multiply(V.Transpose())).Transpose());
                //Matrix s = (Matrix)((U.Transpose()).Multiply(Uk)).Multiply(((Vk.Multiply(V.Transpose())).Transpose()));
                Matrix a = (Matrix)U.Transpose().Multiply(Uk);
                Matrix b = (Matrix)Vk.Multiply(V.Transpose()).Transpose();
                Matrix s = (Matrix)a.PointwiseMultiply(b);
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
                

                //U.SetColumn(U.ColumnCount - 1, Uk.Column(0));
                //V.SetRow(V.RowCount - 1, Vk.Row(0));
                U = AddColumn(U, (Vector)Uk.Column(0));
                //U.SetColumn(k, Uk.Column(0));
                V = AddRow(V, (Vector)Vk.Row(0));
                //V.SetRow(k, Vk.Row(0));

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
                    if (Math.Abs(Rjk[ind, 0]) > max)
                    {
                        max = Math.Abs(Rjk[ind, 0]);
                        row = ind;
                    }
                }

                I[k + 1] = row;
                //i = removeIndex(i,I[k+1]);
                i.Remove(I[k + 1]);
            }
            return new Tuple<Matrix, Matrix>(U, V);
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

            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < M; j++)
                {
                    Z[j, i] = 1.0 / (Math.Abs(m[j] - n[i] + 0.0001));
                }
            }

            return Z;
        }

        #region MultilevelCompres
                                  //multilevel_compress(basis_func_boxes,ix_s,iy_s,iz_s,ix_f,iy_f,iz_f,l,L,ACA_thres,OG_data,EM_data)
        public static void MultilevelCompres(BasicFuncBoxes basicFuncBoxes, double ix_s, double iy_s, double iz_s, double ix_f, double iy_f, double iz_f, int l, double L, double ACA_thres,ref NewSparseMatrix Z_comp)
        {
            bool sym_source_field = true;
            //NewSparseMatrix Z_comp = new NewSparseMatrix();

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
                        List<int> m = new List<int>();
                        for (int i = 0; i < basicFuncBoxes.X.RowCount; i++)
                        {
                            if (basicFuncBoxes.X[i, l] == ix_chs && basicFuncBoxes.Y[i, l] == iy_chs && basicFuncBoxes.Z[i, l] == iz_chs)
                            {
                                m.Add(i);
                            }
                        }
                        if (m.Count == 0)
                        {
                            continue;
                        }

                        //Subdivide field box and process children
                        for (double xchf = 0; xchf <= 1; xchf++)
                        {
                            for (double ychf = 0; ychf <= 1; ychf++)
                            {
                                for (double zchf = 0; zchf <= 1; zchf++)
                                {
                                    //x-index of field child box at level l+1
                                    double ix_chf = ix_f * 2 + xchf;
                                    // y-index of field child box at level l+1
                                    double iy_chf = iy_f * 2 + ychf;
                                    // z-index of field child box at level l+1
                                    double iz_chf = iz_f * 2 + zchf;

                                    //Find indices of testing functions in field child box
                                    List<int> n = new List<int>();
                                    for (int i = 0; i < basicFuncBoxes.X.RowCount; i++)
                                    {
                                        if (basicFuncBoxes.X[i, l] == ix_chf && basicFuncBoxes.Y[i, l] == iy_chf && basicFuncBoxes.Z[i, l] == iz_chf)
                                        {
                                            n.Add(i);
                                        }
                                    }
                                    if (n.Count == 0)
                                    {
                                        continue;
                                    }

                                    // Here we have a pair of non-empty source and field boxes
                                    if (Math.Abs(ix_chs - ix_chf) > 1 || Math.Abs(iy_chs - iy_chf) > 1 || Math.Abs(iz_chs - iz_chf) > 1)
                                    {// Far-field boxes
                                        if (sym_source_field == true)
                                        {
                                            // Symmetric source-field field-source
                                            // interactions are only computed once.
                                            if (ix_chs - ix_chf > 0 || ix_chs - ix_chf == 0 && iy_chs - iy_chf > 0 || ix_chs - ix_chf == 0 && iy_chs - iy_chf == 0 && iz_chs - iz_chf > 0)
                                            {
                                                //[U,V] = ACA(ACA_thres, m,n, OG_data,EM_data);
                                                DenseMatrix U = new DenseMatrix(1, 1);
                                                DenseMatrix V = new DenseMatrix(1, 1);
                                                Tuple<Matrix,Matrix> r = Aca(ACA_thres, m, n, U, V);
                                                //Z_comp{length(Z_comp)+1} = struct('m',m,'n',n,'comp',1,'self',0,'Z',[] ,'U',U,'V',V);
                                                Z_comp.Add(new ACAStruct(m,n,new DenseMatrix(1,1),r.Item1,r.Item2,1,0));
                                            }
                                        }
                                        else
                                        {
                                            // We need to compute all because there is
                                            // not symmetric interactions
                                            DenseMatrix U = new DenseMatrix(1, 1);
                                            DenseMatrix V = new DenseMatrix(1, 1);
                                            //[U,V] = ACA(ACA_thres, m,n, OG_data,EM_data);
                                            Tuple<Matrix, Matrix> r = Aca(ACA_thres, m, n, U, V);
                                            //Z_comp{length(Z_comp)+1} = struct('m',m,'n',n,'comp',1,'self',0,'Z',[] ,'U',U,'V',V);
                                            Z_comp.Add(new ACAStruct(m, n, new DenseMatrix(1, 1), r.Item1, r.Item2, 1, 0));
                                        }
                                    }
                                    else // Near-field boxes
                                    {
                                         double self = 0;
                                         if (l + 1 == L)
                                         {
                                             if (sym_source_field == true)
                                             {
                                                 // Symmetric source-field field-source
                                                 // interactions are only computed once.
                                                 if (ix_chs - ix_chf == 0 && iy_chs - iy_chf == 0 && iz_chs - iz_chf == 0)
                                                 {
                                                     // Self-interactions
                                                     self = 1;
                                                     // Z_comp{length(Z_comp)+1} = struct('m',m,'n',n,'comp',0,'self',self,'Z', user_impedance(m,n,OG_data,EM_data),'U',[],'V',[]);
                                                     Z_comp.Add(new ACAStruct(m, n, UserImpedance(m, n), new DenseMatrix(1, 1), new DenseMatrix(1, 1), 0, self));
                                                 }
                                                 if (ix_chs - ix_chf > 0 || ix_chs - ix_chf == 0 && iy_chs - iy_chf > 0 || ix_chs - ix_chf == 0 && iy_chs - iy_chf == 0 && iz_chs - iz_chf > 0)
                                                 {
                                                     //Z_comp{length(Z_comp)+1} = struct('m',m,'n',n,'comp',0,'self',self,'Z', user_impedance(m,n,OG_data,EM_data),'U',[],'V',[]);
                                                     Z_comp.Add(new ACAStruct(m, n, UserImpedance(m, n), new DenseMatrix(1, 1), new DenseMatrix(1, 1), 0, self));
                                                 }

                                             }
                                             else
                                             {
                                                 // We need to compute all because there are
                                                 // not symmetric interactions
                                                 if (ix_chs - ix_chf == 0 && iy_chs - iy_chf == 0 && iz_chs - iz_chf == 0)
                                                 {
                                                     // Self-interactions
                                                     self = 1;
                                                 }
                                                 //Z_comp{length(Z_comp)+1} = struct('m',m,'n',n,'comp',0,'self',self,'Z', user_impedance(m,n,OG_data,EM_data),'U',[],'V',[]);
                                                 Z_comp.Add(new ACAStruct(m, n, UserImpedance(m, n), new DenseMatrix(1, 1), new DenseMatrix(1, 1), 0, (int)self));
                                             }

                                         }
                                         else
                                         {
                                             //Z_comp = [Z_comp multilevel_compress(basis_func_boxes,ix_chs,iy_chs,iz_chs,ix_chf,iy_chf,iz_chf,l+1,L,ACA_thres,OG_data,EM_data)];
                                                              //multilevel_compress(basis_func_boxes,ix_s,iy_s,iz_s,ix_f,iy_f,iz_f,l,L,ACA_thres,OG_data,EM_data)
                                             MultilevelCompres(basicFuncBoxes,ix_chs,iy_chs,iz_chs,ix_chf,iy_chf,iz_chf,l+1,L,ACA_thres,ref Z_comp);
                                            
                                         }
                                    }
                                
                                }
                            }
                        }
                    }
                }
            }

            //return Z_comp;
        }
        #endregion

        /// <summary>
        /// adds new column to matrix
        /// </summary>
        /// <param name="dest">matrix which add column</param>
        /// <param name="colToAdd">column added to matrix</param>
        /// <returns>new matrix</returns>
        private static Matrix AddColumn(Matrix dest, Vector colToAdd)
        {
            Matrix res = new DenseMatrix(dest.RowCount, dest.ColumnCount + 1);
            res.SetSubMatrix(0, dest.RowCount, 0, dest.ColumnCount, dest);
            res.SetColumn(res.ColumnCount - 1, colToAdd);
            return res;
        }

        /// <summary>
        /// adds new row to matrix
        /// </summary>
        /// <param name="dest">matrix which add row</param>
        /// <param name="rowToAdd">row added to matrix</param>
        /// <returns>new matrix</returns>
        private static Matrix AddRow(Matrix dest, Vector rowToAdd)
        {
            Matrix res = new DenseMatrix(dest.RowCount + 1, dest.ColumnCount);
            res.SetSubMatrix(0, dest.RowCount, 0, dest.ColumnCount, dest);
            res.SetRow(res.RowCount - 1, rowToAdd);
            return res;
        }
    }
}

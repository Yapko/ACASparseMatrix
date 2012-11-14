using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MathNet.Numerics.LinearAlgebra.Double;

namespace ACASparseMatrix
{
    class ACA
    {
        public BasicFuncBoxes prepare_multilevel(Vector rcx, Vector rcy, Vector rcz, int N, double finest_level_size)
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
        public NewSparseMatrix multilevel_compress()
        {
            return new NewSparseMatrix();
        }
    }
}

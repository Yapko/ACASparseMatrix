using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MathNet.Numerics.LinearAlgebra.Double;
namespace ACASparseMatrix
{
    class Program
    {
        static void Main(string[] args)
        {
            BasicFuncBoxes b = new BasicFuncBoxes(10, 10);
            Vector v = new DenseVector(10);
            for (int i = 0; i < v.Count; i++)
            {
                v[i] = -1 + i;
            }

            b.X.SetColumn(0, v);
        }

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

            box_size = box_size / 2;

            Vector tempVx = new DenseVector(rcx.Count);
            for (int i = 0; i < rcx.Count; i++)
            {
                tempVx[i] = Math.Floor((rcx[i] - xmin) / box_size);
            }

            Vector tempVy = new DenseVector(rcy.Count);
            for (int i = 0; i < rcy.Count; i++)
            {
                tempVy[i] = Math.Floor((rcy[i] - ymin) / box_size);
            }

            Vector tempVz = new DenseVector(rcz.Count);
            for (int i = 0; i < rcz.Count; i++)
            {
                tempVz[i] = Math.Floor((rcz[i] - zmin) / box_size);
            }

            for (int i = 0; i < L; i++)
            {
                basicFuncBoxes.X.SetColumn(i, tempVx);
                basicFuncBoxes.Y.SetColumn(i, tempVy);
                basicFuncBoxes.Z.SetColumn(i, tempVz);
            }
            return basicFuncBoxes;
        }
    }
}

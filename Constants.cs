using System;
using System.Collections.Generic;
using System.Text;

namespace OtiPorticos3D
{
    public static class Constants
    {

        public static double Tolerance()
        {
            return Math.Pow(10, -8);
        }
        public static double Gravity()
        {
            return 9.80665;
        }
        public static int ProblemDimension()
        {
            return 3;
        }
    }
}

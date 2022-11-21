using System;
using System.Collections.Generic;
using System.Text;

namespace OtiPorticos3D.Entities
{
    class Parameters
    {
        public double Length { get; set; }
        public double Area { get; set; }
        public double InerciaX { get; set; }
        public double InerciaY { get; set; }
        public double InerciaXY { get; set; }
        public double YoungModule { get; set; }
        public double ShearModule { get; set; }
        public double Density { get; set; }
        public double[] Angle { get; set; }
        public int ElementNumber { get; set; }

        public Parameters(double length, double area, double inerciaX, double inerciaY, double inerciaXY, 
                          double youngModule, double shearModule, double density, double[] angle, int elementNumber)
        {
            Length = length;
            Area = area;
            InerciaX = inerciaX;
            InerciaY = inerciaY;
            InerciaXY = inerciaXY;
            YoungModule = youngModule;
            ShearModule = shearModule;
            Density = density;
            Angle = angle;
            ElementNumber = elementNumber;
        }

        public Parameters()
        {
        }
    }
}

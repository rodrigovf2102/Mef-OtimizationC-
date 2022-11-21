using System;
using System.Collections.Generic;
using System.Text;

namespace OtiPorticos3D.Entities
{
    class Element
    {
        public int ElementNumber { get; set; }
        public int DegreesOfFreedom { get; private set; }
        public double[,] ElementStiffnessMatrix { get; set; }
        public double[,] ElementRotatedStiffnessMatrix { get; set; }
        public double[,] ElementMassMatrix { get; set; }
        public double[,] ElementRotatedMassMatrix { get; set; }
        public double[,] ElementRotationMatrix { get; set; }
        public double[,] ElementTranspostRotationMatrix { get; set; }
        public int NodesPerElement { get; set; }    
        public Node[] Nodes { get; set; }
        public Parameters ElementParameters { get; set; }

        public Element(int elementNumber, Parameters elementParameters, Node[] nodes)
        {
            DegreesOfFreedom = NodesPerElement * 6;
            ElementNumber = elementNumber;
            ElementRotationMatrix = new double[DegreesOfFreedom, DegreesOfFreedom];
            ElementTranspostRotationMatrix = new double[DegreesOfFreedom, DegreesOfFreedom];
            ElementStiffnessMatrix = new double[DegreesOfFreedom, DegreesOfFreedom];
            ElementRotatedStiffnessMatrix = new double[DegreesOfFreedom, DegreesOfFreedom];
            ElementMassMatrix = new double[DegreesOfFreedom, DegreesOfFreedom];
            ElementRotatedMassMatrix = new double[DegreesOfFreedom, DegreesOfFreedom];
            NodesPerElement = 2;
            Nodes = nodes;
            ElementParameters = elementParameters;
        }

    }
}

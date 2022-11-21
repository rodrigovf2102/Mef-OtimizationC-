using System;
using System.Collections.Generic;
using System.Text;

namespace OtiPorticos3D.Entities
{
    class Node
    {
        public int NodeNumber { get; set; }
        public int NodeDimension { get; set; }
        public int DegreesOfFreedom { get; private set; }
        public bool NodeWithLoad { get; set; }
        public bool NodeWithSpring { get; set; }
        public double[] NodeLoadsValues { get; set; }
        public double[] NodeSpringsValues { get; set; }
        public double[] NodeCoordinates { get; set; }
        public Element NodeElement { get; set; }

        public Node(int nodeNumber, int nodeDimension, bool nodeWithLoad, bool nodeWithSpring, 
                    double[] nodeLoadsValues, double[] nodeSpringsValues, double[] nodeCoordinates)
        {
            NodeNumber = nodeNumber;
            NodeDimension = nodeDimension;;
            DegreesOfFreedom = 6;
            NodeWithLoad = nodeWithLoad;
            NodeWithSpring = nodeWithSpring;
            NodeLoadsValues = nodeLoadsValues;
            NodeSpringsValues = nodeSpringsValues;
            NodeCoordinates = nodeCoordinates;
        }
    }
}

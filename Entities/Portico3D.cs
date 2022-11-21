using System;
using OtiPorticos3D.Entities;

namespace OtiPorticos3D
{
    class Portico3D
    {
        public int PorticoDimension { get; set; }
        public Element[] Elements { get; set; }
        public int ElementsAmount { get; set; }
        public int NodesAmount { get; set; }
        public int DegreesOfFreedom { get; private set; }
        public int NumberOfNodesWithLoads { get; set; }
        public int NumberOfNodesWithSprings { get; set; }
        public int[,] ElementsConectivity { get; set; }
        public double[,] PorticoStiffnessMatrix { get; set; }
        public double[,] PorticoInverseStiffnessMatrix { get; set; }
        public double[,] PorticoMassMatrix { get; set; }
        public double[] PorticoLoadsVector { get; set; }
        public double[,] PorticoModalMatrix { get; set; }
        public int NumberOfVibrationModesShown { get; set; }
        public double[] NaturalFrequencies { get; set; }
        public double[,] VibrationModes { get; set; }

        public Portico3D()
        {
        }

        public Portico3D(int porticoDimension, Element[] elements, int elementsAmount, int nodesAmount,
                         int numberOfNodesWithLoads, int numberOfNodesWithSprings, int[,] elementsConectivity,
                         double[] porticoLoadsVector,int numberOfVibrationModesShown)
        {
            PorticoDimension = porticoDimension;
            Elements = elements;
            ElementsAmount = elementsAmount;
            NodesAmount = nodesAmount;
            DegreesOfFreedom = 6 * NodesAmount;
            NumberOfNodesWithLoads = numberOfNodesWithLoads;
            NumberOfNodesWithSprings = numberOfNodesWithSprings;
            ElementsConectivity = elementsConectivity;
            PorticoStiffnessMatrix = new double[DegreesOfFreedom, DegreesOfFreedom];
            PorticoInverseStiffnessMatrix = new double[DegreesOfFreedom, DegreesOfFreedom];
            PorticoMassMatrix = new double[DegreesOfFreedom, DegreesOfFreedom];
            PorticoLoadsVector = porticoLoadsVector;
            PorticoModalMatrix = new double[DegreesOfFreedom, DegreesOfFreedom];
            NumberOfVibrationModesShown = numberOfVibrationModesShown;
            NaturalFrequencies = new double[DegreesOfFreedom];
            VibrationModes = new double[DegreesOfFreedom, DegreesOfFreedom];
        }
    }
}

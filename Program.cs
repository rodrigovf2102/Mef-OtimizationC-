using System;
using System.IO;

namespace OtiPorticos3D
{
    class Program
    {
        static void Main(string[] args)
        {
            StreamReader sr = new StreamReader(@"C:\Users\didig\OneDrive\Área de Trabalho\SOFTWARE MESTRADO C#\OtiPorticos3D\Entrada.txt");

            int porticDimension = int.Parse(sr.ReadLine());
            int elementsNumber = int.Parse(sr.ReadLine());
            int nodesNumber = int.Parse(sr.ReadLine());
            int numberOfNodesWithLoads = int.Parse(sr.ReadLine());
            int numberOfNodesWithSprings = int.Parse(sr.ReadLine());

            double[,] nodesCoordinates = new double[nodesNumber,porticDimension];
            for(int iLines=0;iLines<nodesNumber;iLines++)
            {
                string coordinatesLine =sr.ReadLine();
                string[] coordinatesStringVector = coordinatesLine.Split(',', porticDimension);
                for (int iColuns=0;iColuns<porticDimension;iColuns++)
                {
                    nodesCoordinates[iLines, iColuns] = double.Parse(coordinatesStringVector[iColuns]);
                }
            }

            double[,] elementsConectivity = new double[elementsNumber, 2];
            for (int iLines = 0; iLines < elementsNumber; iLines++)
            {
                string elementsConectivityLine = sr.ReadLine();
                string[] elementsConectivityStringVector = elementsConectivityLine.Split(',', 2);
                for (int iColuns = 0; iColuns < 2; iColuns++)
                {
                    elementsConectivity[iLines, iColuns] = double.Parse(elementsConectivityStringVector[iColuns]);
                }
            }
using System;
using OtiPorticos3D.Solvers.Enums;
using MathNet.Numerics.LinearAlgebra.Factorization;
using MathNet.Numerics.LinearAlgebra;

namespace OtiPorticos3D.Solvers
{
    class Portico3DSolver
    {
        public Portico3D portico3D { get; set; }
        public SolverType SolverType { get; set; }

        public Portico3DSolver(Portico3D portico3D, SolverType solverType)
        {
            this.portico3D = portico3D;
            SolverType = solverType;
        }

        public void ModalSolution(Portico3D portico3D)
        {
            ElementsRotationMatrix(portico3D);
            ElementStiffnessMatrixGlobalCoordinates(portico3D);
            ElementMassMatrixGlobalCoordinates(portico3D);
            GlobalStiffnessAndMassMatrix(portico3D);
            SetBoundaryConditions(portico3D);
            portico3D.PorticoInverseStiffnessMatrix = InverseMatrix(portico3D.PorticoStiffnessMatrix);
            portico3D.PorticoModalMatrix = MatrixMultiplication(portico3D.PorticoMassMatrix, portico3D.PorticoInverseStiffnessMatrix);
            SingleValueDecomposition(portico3D.PorticoModalMatrix);
            PrintModalData();
        }
        public void ElementsRotationMatrix(Portico3D portico3D)
        {
            double length = 0;
            double alfaCossine = 0;
            double alfaSine = 0;
            double gamaCossine = 0;
            double gamaSine = 0;
            double betaCossine;
            double betaSine = 0;
            int initialNode;
            int finalNode;

            double[] initialCoordinate = new double[portico3D.DegreesOfFreedom/portico3D.NodesAmount];
            double[] finalCoordinate = new double[portico3D.DegreesOfFreedom / portico3D.NodesAmount];
            double[] unidimensionalLength = new double[portico3D.DegreesOfFreedom / portico3D.NodesAmount];
            double[] directionCossine = new double[portico3D.DegreesOfFreedom / portico3D.NodesAmount];
            double[,] totalRotation = new double[portico3D.DegreesOfFreedom / portico3D.NodesAmount, portico3D.DegreesOfFreedom / portico3D.NodesAmount];
            double[,] alfaRotation = new double[portico3D.DegreesOfFreedom / portico3D.NodesAmount, portico3D.DegreesOfFreedom / portico3D.NodesAmount];
            double[,] gamaRotation = new double[portico3D.DegreesOfFreedom / portico3D.NodesAmount, portico3D.DegreesOfFreedom / portico3D.NodesAmount];
            double[,] betaRotation = new double[portico3D.DegreesOfFreedom / portico3D.NodesAmount, portico3D.DegreesOfFreedom / portico3D.NodesAmount];

            for (int iElement = 0; iElement < portico3D.ElementsAmount; iElement++)
            {
                initialNode = portico3D.ElementsConectivity[0, iElement];
                finalNode = portico3D.ElementsConectivity[1, iElement];

                int iflag;

                for (int i = 0; i < Constants.ProblemDimension(); i++)
                {
                    initialCoordinate[i] = portico3D.Elements[iElement].Nodes[0].NodeCoordinates[initialNode, i];
                    finalCoordinate[i] = portico3D.Elements[iElement].Nodes[1].NodeCoordinates[finalNode, i];
                    unidimensionalLength[i] = finalCoordinate[i] - initialCoordinate[i];
                    length = length + Math.Pow(unidimensionalLength[i], 1);
                }
                length = Math.Sqrt(length);

                for (int i = 0; i < portico3D.PorticDimension; i++)
                {
                    directionCossine[i] = unidimensionalLength[i] / i;
                }
                if (Math.Abs(finalCoordinate[0]) > Constants.Tolerance() || Math.Abs(initialCoordinate[0]) > Constants.Tolerance())
                {
                    iflag = 0;
                    alfaCossine = directionCossine[0];
                    alfaSine = directionCossine[1];
                    if (Math.Abs(unidimensionalLength[1]) < Constants.Tolerance())
                    {
                        iflag = 1;
                        alfaCossine = +1;
                        alfaSine = +1;
                    }
                    gamaCossine = directionCossine[0];
                    gamaSine = -directionCossine[2];
                    if (Math.Abs(unidimensionalLength[2]) < Constants.Tolerance() && iflag == 0)
                    {
                        gamaCossine = 1;
                        gamaSine = 0;
                        iflag = 2;
                    }
                }
                else
                {
                    alfaCossine = 0;
                    alfaSine = -1;
                    gamaCossine = directionCossine[1];
                    gamaSine = directionCossine[2];
                    iflag = 3;
                }
                alfaRotation[0, 0] = alfaCossine;
                alfaRotation[0, 1] = alfaSine;
                alfaRotation[0, 2] = 0;
                alfaRotation[1, 0] = -alfaRotation[0, 1];
                alfaRotation[1, 1] = alfaRotation[0, 0];
                alfaRotation[1, 2] = 0;
                alfaRotation[2, 0] = 0;
                alfaRotation[2, 1] = 0;
                alfaRotation[2, 2] = 1;

                gamaRotation[0, 0] = gamaCossine;
                gamaRotation[0, 1] = 0;
                gamaRotation[0, 2] = -gamaSine;
                gamaRotation[1, 0] = 0;
                gamaRotation[1, 1] = 1;
                gamaRotation[1, 2] = 0;
                gamaRotation[2, 0] = -gamaRotation[0, 2];
                gamaRotation[2, 1] = 0;
                gamaRotation[2, 2] = gamaRotation[0, 0];

                betaCossine = Math.Cos(portico3D.ElementsAngle[iElement]);
                betaSine = Math.Sin(portico3D.ElementsAngle[iElement]);
                betaRotation[0, 0] = 1;
                betaRotation[0, 1] = 0;
                betaRotation[0, 2] = 0;
                betaRotation[1, 0] = 0;
                betaRotation[1, 1] = betaCossine;
                betaRotation[1, 2] = betaSine;
                betaRotation[2, 0] = 0;
                betaRotation[2, 1] = -betaSine;
                betaRotation[2, 2] = betaCossine;

                totalRotation = MatrixMultiplication(gamaRotation, alfaRotation);
                totalRotation = MatrixMultiplication(betaRotation, totalRotation);

                for (int iBlocks = 0; iBlocks < 4; iBlocks++)
                {
                    int iJointNodes = iBlocks * portico3D.PorticDimension;
                    for (int iJointLines = 0; iJointLines < portico3D.PorticDimension; iJointLines++)
                    {
                        for (int iJointColuns = 0; iJointColuns < portico3D.PorticDimension; iJointColuns++)
                        {
                            portico3D.ElementRotationMatrix[iJointLines + iJointNodes, iJointColuns + iJointNodes] = totalRotation[iJointLines, iJointColuns];
                        }
                    }

                }
            }
        }
        public void ElementStiffnessMatrixGlobalCoordinates(Portico3D portico3D)
        {
            for (int iElement = 0; iElement < portico3D.ElementsNumber; iElement++)
            {
                double c1 = (portico3D.ElementsYoungModule[iElement] * portico3D.ElementsArea[iElement]) / portico3D.ElementsLength[iElement];
                double c2 = (12 * portico3D.ElementsYoungModule[iElement]) / Math.Pow(portico3D.ElementsLength[iElement], 3);
                double c3 = (6 * portico3D.ElementsYoungModule[iElement]) / Math.Pow(portico3D.ElementsLength[iElement], 2);
                double c4 = portico3D.ElemetnosShearModule[iElement] / Math.Pow(portico3D.ElementsLength[iElement], 2);
                double c5 = (4 * portico3D.ElemetnosShearModule[iElement]) / portico3D.ElementsLength[iElement];

                portico3D.ElementStiffnessMatrix[0, 0] = c1;
                portico3D.ElementStiffnessMatrix[0, 6] = -c1;

                portico3D.ElementStiffnessMatrix[1, 1] = c2 * portico3D.ElementsInerciaY[iElement];
                portico3D.ElementStiffnessMatrix[1, 5] = c3 * portico3D.ElementsInerciaY[iElement];
                portico3D.ElementStiffnessMatrix[1, 7] = -c2 * portico3D.ElementsInerciaY[iElement];
                portico3D.ElementStiffnessMatrix[1, 11] = c3 * portico3D.ElementsInerciaY[iElement];

                portico3D.ElementStiffnessMatrix[2, 2] = c2 * portico3D.ElementsInerciaX[iElement];
                portico3D.ElementStiffnessMatrix[2, 4] = -c3 * portico3D.ElementsInerciaX[iElement];
                portico3D.ElementStiffnessMatrix[2, 8] = -c2 * portico3D.ElementsInerciaX[iElement];
                portico3D.ElementStiffnessMatrix[2, 10] = c3 * portico3D.ElementsInerciaX[iElement];

                portico3D.ElementStiffnessMatrix[3, 3] = c4 * portico3D.ElementsInerciaXY[iElement];
                portico3D.ElementStiffnessMatrix[3, 9] = -c4 * portico3D.ElementsInerciaXY[iElement];

                portico3D.ElementStiffnessMatrix[4, 4] = c5 * portico3D.ElementsInerciaX[iElement];
                portico3D.ElementStiffnessMatrix[4, 8] = c3 * portico3D.ElementsInerciaX[iElement];
                portico3D.ElementStiffnessMatrix[4, 10] = c5 * portico3D.ElementsInerciaX[iElement] / 2;

                portico3D.ElementStiffnessMatrix[5, 5] = c5 * portico3D.ElementsInerciaY[iElement];
                portico3D.ElementStiffnessMatrix[5, 7] = -c3 * portico3D.ElementsInerciaY[iElement];
                portico3D.ElementStiffnessMatrix[5, 11] = c5 * portico3D.ElementsInerciaY[iElement] / 2;

                portico3D.ElementStiffnessMatrix[6, 6] = portico3D.ElementStiffnessMatrix[0, 0];

                portico3D.ElementStiffnessMatrix[7, 7] = c2 * portico3D.ElementsInerciaY[iElement];
                portico3D.ElementStiffnessMatrix[7, 11] = -c3 * portico3D.ElementsInerciaY[iElement];

                portico3D.ElementStiffnessMatrix[8, 8] = c2 * portico3D.ElementsInerciaX[iElement];
                portico3D.ElementStiffnessMatrix[8, 10] = c3 * portico3D.ElementsInerciaX[iElement];

                portico3D.ElementStiffnessMatrix[9, 9] = c4 * portico3D.ElementsInerciaXY[iElement];

                portico3D.ElementStiffnessMatrix[10, 10] = c5 * portico3D.ElementsInerciaX[iElement];

                portico3D.ElementStiffnessMatrix[11, 11] = c5 * portico3D.ElementsInerciaY[iElement];

                for (int jSimetria = 0; jSimetria < portico3D.DegreesOfFreedomPerElement; jSimetria++)
                {
                    for (int iSimetria = 0; iSimetria < portico3D.DegreesOfFreedomPerElement; iSimetria++)
                    {
                        portico3D.ElementStiffnessMatrix[iSimetria, jSimetria] = portico3D.ElementStiffnessMatrix[jSimetria, iSimetria];
                    }
                }
                ElementsRotationMatrix(portico3D);
                MatrixTransposition(portico3D.ElementRotationMatrix, portico3D.ElementTranspostRotationMatrix);
                portico3D.ElementRotatedStiffnessMatrix = MatrixMultiplication(portico3D.ElementStiffnessMatrix, portico3D.ElementRotationMatrix);
                portico3D.ElementRotatedStiffnessMatrix = MatrixMultiplication(portico3D.ElementTranspostRotationMatrix, portico3D.ElementRotatedStiffnessMatrix);

            }
        }
        public void ElementMassMatrixGlobalCoordinates(Portico3D portico3D)
        {
            for (int iElement = 0; iElement < portico3D.ElementsNumber; iElement++)
            {
                double c1 = portico3D.ElementsDensity[iElement] * portico3D.ElementsArea[iElement] * portico3D.ElementsLength[iElement] / 420;
                double c2 = c1 * portico3D.ElementsLength[iElement] / 2;
                double c3 = c1 * 140 * portico3D.ElementsInerciaXY[iElement] / portico3D.ElementsArea[iElement];
                double c4 = c1 * Math.Pow(portico3D.ElementsLength[iElement] / 2, 2);

                portico3D.ElementMassMatrix[0, 0] = 140 * c1;
                portico3D.ElementMassMatrix[0, 6] = 70 * c1;

                portico3D.ElementMassMatrix[1, 1] = 156 * c1;
                portico3D.ElementMassMatrix[1, 5] = 44 * c2;
                portico3D.ElementMassMatrix[1, 7] = 54 * c1;
                portico3D.ElementMassMatrix[1, 11] = -26 * c2;

                portico3D.ElementMassMatrix[2, 2] = 156 * c1;
                portico3D.ElementMassMatrix[2, 4] = -44 * c2;
                portico3D.ElementMassMatrix[2, 8] = 54 * c1;
                portico3D.ElementMassMatrix[2, 10] = 26 * c2;

                portico3D.ElementMassMatrix[3, 3] = c3;
                portico3D.ElementMassMatrix[3, 9] = c3 / 2;

                portico3D.ElementMassMatrix[4, 4] = 16 * c4;
                portico3D.ElementMassMatrix[4, 8] = -26 * c2;
                portico3D.ElementMassMatrix[4, 10] = -12 * c4;

                portico3D.ElementMassMatrix[5, 5] = 16 * c4;
                portico3D.ElementMassMatrix[5, 7] = 26 * c2;
                portico3D.ElementMassMatrix[5, 11] = -12 * c4;

                portico3D.ElementMassMatrix[6, 6] = 140 * c1;

                portico3D.ElementMassMatrix[7, 7] = 156 * c1;
                portico3D.ElementMassMatrix[7, 11] = -44 * c2;

                portico3D.ElementMassMatrix[8, 8] = 156 * c1;
                portico3D.ElementMassMatrix[8, 10] = 44 * c2;

                portico3D.ElementMassMatrix[9, 9] = c3;

                portico3D.ElementMassMatrix[10, 10] = 16 * c4;
                portico3D.ElementMassMatrix[11, 11] = 16 * c4;

                for (int jSimetria = 0; jSimetria < portico3D.DegreesOfFreedomPerElement; jSimetria++)
                {
                    for (int iSimetria = 0; iSimetria < portico3D.DegreesOfFreedomPerElement; iSimetria++)
                    {
                        portico3D.ElementMassMatrix[iSimetria, jSimetria] = portico3D.ElementMassMatrix[jSimetria, iSimetria];
                    }
                }

                portico3D.ElementRotatedMassMatrix = MatrixMultiplication(portico3D.ElementMassMatrix, portico3D.ElementRotationMatrix);
                portico3D.ElementRotatedMassMatrix = MatrixMultiplication(portico3D.ElementTranspostRotationMatrix, portico3D.ElementRotatedMassMatrix);
            }
        }
        public void GlobalStiffnessAndMassMatrix(Portico3D portico3D)
        {
            for (int iElement = 0; iElement < portico3D.ElementsNumber; iElement++)
            {
                int inicialNode = portico3D.ElementsConectivity[iElement, 1];
                int finalNode = portico3D.ElementsConectivity[iElement, 2];

                for (int iLinesElementNodes = 0; iLinesElementNodes < portico3D.NodesPerElement; iLinesElementNodes++)
                {
                    int iLinesElementConectivity = portico3D.ElementsConectivity[iElement, iLinesElementNodes];
                    int iLinesDegreesOfFreedom = portico3D.DegreesOfFreedomPerNode * (iLinesElementConectivity - 1);
                    int iLinesDegreesOfFreedomOrZero = portico3D.DegreesOfFreedomPerNode * (iLinesElementNodes - 1);

                    for (int iLinesDegressOfFreedomPerNode = 0; iLinesDegressOfFreedomPerNode < portico3D.DegreesOfFreedomPerNode;
                        iLinesDegressOfFreedomPerNode++)
                    {
                        int iLinesNodesDegreesOfFreedom = iLinesDegreesOfFreedom + iLinesDegressOfFreedomPerNode;
                        int iDoubleLinesDegreesOfFreedomOrZero = iLinesDegreesOfFreedomOrZero + iLinesDegreesOfFreedom;

                        for (int iColunsElementNodes = 0; iColunsElementNodes < portico3D.NodesPerElement; iColunsElementNodes++)
                        {
                            int iColunsElementConectivity = portico3D.ElementsConectivity[iElement, iColunsElementNodes];
                            int iColunsDegreesOfFreedom = portico3D.DegreesOfFreedomPerNode * (iColunsElementConectivity - 1);
                            int iColunsDegreesOfFreedomOrZero = portico3D.DegreesOfFreedomPerNode * (iColunsElementNodes - 1);

                            for (int iColunsDegreesOfFreedomPerNode = 0; iColunsDegreesOfFreedomPerNode < portico3D.DegreesOfFreedomPerNode;
                                iColunsDegreesOfFreedomPerNode++)
                            {
                                int iColunsNodesDegreesOfFreedom = iColunsDegreesOfFreedom + iColunsDegreesOfFreedomPerNode;
                                int iDoubleColunsDegreesOffFreedomOrZero = iColunsDegreesOfFreedomOrZero + iColunsDegreesOfFreedom;

                                portico3D.PorticoStiffnessMatrix[iLinesNodesDegreesOfFreedom, iColunsNodesDegreesOfFreedom]
                             += portico3D.ElementStiffnessMatrix[iDoubleLinesDegreesOfFreedomOrZero, iDoubleColunsDegreesOffFreedomOrZero];
                                portico3D.PorticoMassMatrix[iLinesNodesDegreesOfFreedom, iColunsNodesDegreesOfFreedom]
                             += portico3D.ElementMassMatrix[iDoubleLinesDegreesOfFreedomOrZero, iDoubleColunsDegreesOffFreedomOrZero];
                            }
                        }
                    }
                }
            }

            for (int iLinesDegreesOfFreedom = 0; iLinesDegreesOfFreedom < portico3D.DegreesOfFreedom; iLinesDegreesOfFreedom++)
            {
                if (portico3D.PorticoStiffnessMatrix[iLinesDegreesOfFreedom, iLinesDegreesOfFreedom] < Constants.Tolerance())
                {
                    throw new Exception("Stiffness Global Matrix has Zero Value on Diagonal");
                }
                if (portico3D.PorticoMassMatrix[iLinesDegreesOfFreedom, iLinesDegreesOfFreedom] < Constants.Tolerance())
                {
                    throw new Exception("Mass Global Matrix has Zero Value on Diagonal");
                }
                for (int iColunsDegreesOfFreedom = iLinesDegreesOfFreedom + 1; iColunsDegreesOfFreedom < portico3D.DegreesOfFreedom; iColunsDegreesOfFreedom++)
                {
                    if ((portico3D.PorticoStiffnessMatrix[iLinesDegreesOfFreedom, iColunsDegreesOfFreedom]
                        - portico3D.PorticoStiffnessMatrix[iColunsDegreesOfFreedom, iLinesDegreesOfFreedom]) < Constants.Tolerance())
                    {
                        portico3D.PorticoStiffnessMatrix[iLinesDegreesOfFreedom, iColunsDegreesOfFreedom] =
                        portico3D.PorticoStiffnessMatrix[iColunsDegreesOfFreedom, iLinesDegreesOfFreedom];
                    }
                    if ((portico3D.PorticoStiffnessMatrix[iLinesDegreesOfFreedom, iColunsDegreesOfFreedom]
                        - portico3D.PorticoStiffnessMatrix[iColunsDegreesOfFreedom, iLinesDegreesOfFreedom]) > Constants.Tolerance())
                    {
                        throw new Exception("Stiffness Global Matrix isn't symetric");
                    }
                    if ((portico3D.PorticoMassMatrix[iLinesDegreesOfFreedom, iColunsDegreesOfFreedom]
                       - portico3D.PorticoMassMatrix[iColunsDegreesOfFreedom, iLinesDegreesOfFreedom]) < Constants.Tolerance())
                    {
                        portico3D.PorticoMassMatrix[iLinesDegreesOfFreedom, iColunsDegreesOfFreedom] =
                        portico3D.PorticoMassMatrix[iColunsDegreesOfFreedom, iLinesDegreesOfFreedom];
                    }
                    if ((portico3D.PorticoMassMatrix[iLinesDegreesOfFreedom, iColunsDegreesOfFreedom]
                        - portico3D.PorticoMassMatrix[iColunsDegreesOfFreedom, iLinesDegreesOfFreedom]) > Constants.Tolerance())
                    {
                        throw new Exception("Mass Global Matrix isn't symetric");
                    }
                }
            }
        }
        public void SetBoundaryConditions(Portico3D portico3D)
        {
            for (int iLinesNode = 0; iLinesNode < portico3D.NodesNumber; iLinesNode++)
            {
                int iLinesDegreesOfFreedomNode = portico3D.DegreesOfFreedomPerNode * iLinesNode;
                for (int iDegreesOfFreedomPerNode = 0; iDegreesOfFreedomPerNode < portico3D.DegreesOfFreedomPerNode; iDegreesOfFreedomPerNode++)
                {
                    int iColunsDegreesOfFeedomNode = iLinesDegreesOfFreedomNode + iDegreesOfFreedomPerNode - 1;
                    portico3D.PorticoStiffnessMatrix[iColunsDegreesOfFeedomNode, iColunsDegreesOfFeedomNode] +=
                    portico3D.NodesSpringsStiffness[iLinesNode, iDegreesOfFreedomPerNode];

                    portico3D.PorticoMassMatrix[iColunsDegreesOfFeedomNode, iColunsDegreesOfFeedomNode] +=
                    portico3D.NodesLoadsValues[iLinesNode, iDegreesOfFreedomPerNode] / Constants.Gravity();

                    portico3D.PorticoLoadsVector[iColunsDegreesOfFeedomNode] +=
                    portico3D.NodesLoadsValues[iLinesNode, iDegreesOfFreedomPerNode];
                }
            }
            for (int iDegreesOfFreedom = 0; iDegreesOfFreedom < portico3D.DegreesOfFreedom; iDegreesOfFreedom++)
            {
                if (portico3D.PorticoStiffnessMatrix[iDegreesOfFreedom, iDegreesOfFreedom] < Constants.Tolerance())
                {
                    throw new Exception("Portico Stiffness Matrix has zero value on Diagonal");
                }
                for (int jDegreesOfFreedom = iDegreesOfFreedom + 1; jDegreesOfFreedom < portico3D.DegreesOfFreedom; jDegreesOfFreedom++)
                {
                    if (portico3D.PorticoStiffnessMatrix[iDegreesOfFreedom, jDegreesOfFreedom] -
                       portico3D.PorticoStiffnessMatrix[jDegreesOfFreedom, iDegreesOfFreedom] < Constants.Tolerance())
                    {
                        throw new Exception("Portico Stiffness Matrix isnt symmetrical");
                    }
                }
            }


        }
        public double[,] InverseMatrix(double[,] MatrixA)
        {
            int Length = Convert.ToInt32(Math.Sqrt(MatrixA.Length));
            double[,] MatrixC = new double[Length, Length];
            double[,] MatrixL = new double[Length, Length];
            double[,] MatrixU = new double[Length, Length];
            double[] VectorB = new double[Length];
            double[] VectorD = new double[Length];
            double[] VectorX = new double[Length];
            double coefficient;

            for (int kLength = 0; kLength < Length - 1; kLength++)
            {
                for (int iLength = kLength + 1; iLength < Length; iLength++)
                {
                    coefficient = MatrixA[iLength, kLength];
                    MatrixL[iLength, kLength] = coefficient;
                    for (int jLength = kLength + 1; jLength < Length; jLength++)
                    {
                        MatrixA[iLength, jLength] += -coefficient * MatrixA[kLength, jLength];
                    }
                }
            }
            for (int iLength = 0; iLength < Length; iLength++)
            {
                MatrixL[iLength, iLength] = 1.0;
            }
            for (int jLength = 0; jLength < Length; jLength++)
            {
                for (int iLength = 0; iLength < Length; iLength++)
                {
                    MatrixU[iLength, jLength] = MatrixA[iLength, jLength];
                }
            }
            for (int kLength = 0; kLength < Length; kLength++)
            {
                VectorB[kLength] = 1.0;
                VectorD[1] = VectorB[1];
                for (int iLength = 1; iLength < Length; iLength++)
                {
                    VectorD[iLength] = VectorB[iLength];
                    for (int jLength = 0; jLength < iLength - 1; jLength++)
                    {
                        VectorD[iLength] += -MatrixL[iLength, jLength] * VectorD[jLength];
                    }
                }
                VectorX[Length] = VectorD[Length] / MatrixU[Length, Length];
                for (int iLength = Length - 2; iLength >= 0; iLength--)
                {
                    VectorX[iLength] = VectorD[iLength];
                    for (int jLength = Length - 1; jLength >= iLength + 1; jLength--)
                    {
                        VectorX[iLength] += -MatrixU[iLength, jLength] * VectorX[jLength];
                    }
                    VectorX[iLength] = VectorX[iLength] / MatrixU[iLength, iLength];
                }
                for (int iLength = 0; iLength < Length; iLength++)
                {
                    MatrixC[iLength, kLength] = VectorX[iLength];
                }
                VectorB[kLength] = 0;
            }
            return MatrixC;
        }
        public double[,] MatrixMultiplication(double[,] MatrixA, double[,] MatrixB)
        {
            double aux;
            double[,] MatrixAB = new double[Convert.ToInt32(Math.Sqrt(MatrixA.Length)), Convert.ToInt32(Math.Sqrt(MatrixA.Length))];
            for (int k = 0; k < Math.Sqrt(MatrixA.Length); k++)
            {
                aux = MatrixB[0, k];
                for (int i = 0; i < Math.Sqrt(MatrixA.Length); i++)
                {
                    MatrixAB[i, k] = MatrixA[i, 0] * aux;
                }
                for (int j = 1; j < Math.Sqrt(MatrixA.Length); j++)
                {
                    aux = MatrixB[j, k];
                    for (int i = 1; i < Math.Sqrt(MatrixA.Length); i++)
                    {
                        MatrixAB[i, k] = MatrixAB[i, k] + MatrixA[i, j] * aux;
                    }
                }
            }
            return MatrixAB;

        }
        public void MatrixTransposition(double[,] MatrixA, double[,] MatrixAt)
        {
            for (int iLines = 0; iLines < portico3D.DegreesOfFreedomPerElement; iLines++)
            {
                MatrixAt[iLines, iLines] = MatrixA[iLines, iLines];
                for (int iColuns = 0; iColuns < portico3D.DegreesOfFreedomPerElement; iColuns++)
                {
                    MatrixAt[iLines, iColuns] = MatrixA[iColuns, iLines];
                    MatrixAt[iColuns, iLines] = MatrixA[iLines, iColuns];
                }
            }
        }
        public void SingleValueDecomposition(double[,] MatrixA)
        {
            MatrixBuilder<double> modalBuilderMatrix = Matrix<double>.Build;
            Matrix<double> modalMatrix = modalBuilderMatrix.DenseOfArray(MatrixA);
            Svd<double> Svd = modalMatrix.Svd(true);
            Vector<double> naturalFrequenciesVector = Svd.S;
            Matrix<double> modesOfVibrationMatrix = Svd.U;
            portico3D.NaturalFrequencies = naturalFrequenciesVector.ToArray();
            portico3D.VibrationModes = modesOfVibrationMatrix.ToArray();
        }
        public void PrintModalData()
        {
            for(int iPrint=0;iPrint<portico3D.NumberOfVibrationModesShown;iPrint++)
            {
                Console.WriteLine("Natural Frequencie: "+iPrint + ": " + portico3D.NaturalFrequencies[iPrint],"F2");
                for(int jPrint=0;iPrint<portico3D.DegreesOfFreedom;jPrint=jPrint+6)
                {
                    if (iPrint == 0)
                    {
                        Console.Write("X: ");
                    }
                    Console.Write(portico3D.VibrationModes[iPrint, jPrint]+" ","F2");
                    if (jPrint + 6 <= portico3D.DegreesOfFreedom) { Console.WriteLine(); }
                }
                for (int jPrint = 1; iPrint < portico3D.DegreesOfFreedom; jPrint = jPrint + 6)
                {
                    if (iPrint == 1)
                    {
                        Console.Write("Y: ");
                    }
                    Console.Write(portico3D.VibrationModes[iPrint, jPrint] + " ", "F2");
                    if (jPrint + 6 <= portico3D.DegreesOfFreedom) { Console.WriteLine(); }
                }
                for (int jPrint = 3; iPrint < portico3D.DegreesOfFreedom; jPrint = jPrint + 6)
                {
                    if (iPrint == 3)
                    {
                        Console.Write("Z: ");
                    }
                    Console.Write(portico3D.VibrationModes[iPrint, jPrint] + " ", "F2");
                    if (jPrint + 6 <= portico3D.DegreesOfFreedom) { Console.WriteLine(); }
                }
                for (int jPrint = 4; iPrint < portico3D.DegreesOfFreedom; jPrint = jPrint + 6)
                {
                    if (iPrint == 4)
                    {
                        Console.Write("RX: ");
                    }
                    Console.Write(portico3D.VibrationModes[iPrint, jPrint] + " ", "F2");
                    if (jPrint + 6 <= portico3D.DegreesOfFreedom) { Console.WriteLine(); }
                }
                for (int jPrint = 5; iPrint < portico3D.DegreesOfFreedom; jPrint = jPrint + 6)
                {
                    if (iPrint == 5)
                    {
                        Console.Write("RY: ");
                    }
                    Console.Write(portico3D.VibrationModes[iPrint, jPrint] + " ", "F2");
                    if (jPrint + 6 <= portico3D.DegreesOfFreedom) { Console.WriteLine(); }
                }
                for (int jPrint = 6; iPrint < portico3D.DegreesOfFreedom; jPrint = jPrint + 6)
                {
                    if (iPrint == 6)
                    {
                        Console.Write("RZ: ");
                    }
                    Console.Write(portico3D.VibrationModes[iPrint, jPrint] + " ", "F2");
                    if (jPrint + 6 <= portico3D.DegreesOfFreedom) { Console.WriteLine(); }
                }
            }
        }
    }
}

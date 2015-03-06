/**
    This file is part of the evaluation of the 3D human shape model as described in the paper:

    Leonid Pishchulin, Stefanie Wuhrer, Thomas Helten, Christian Theobalt and Bernt Schiele
    Building Statistical Shape Spaces for 3D Human Modeling
    ArXiv, March 2015

    Please cite the paper if you are using this code in your work.
    
    Contributors: Stefanie Wuhrer, Leonid Pishchulin.

    The code may be used free of charge for non-commercial and
    educational purposes, the only requirement is that this text is
    preserved within the derivative work. For any other purpose you
    must contact the authors for permission. This code may not be
    redistributed without permission from the authors.
*/

#include "GeneralizedProcrustes.h"


GeneralizedProcrustes::GeneralizedProcrustes()
{
        partialAverages = NULL;

        shapes.clear();
}

GeneralizedProcrustes::GeneralizedProcrustes(vector < Shape >& shapes, int numberOfShapes, int minNumNodes)
{
        int i;

        this->numberOfShapes = numberOfShapes;
        this->minNumNodes = minNumNodes;
        partialAverages = NULL;

        this->shapes.clear();
        
        for(i = 0; i < numberOfShapes; i++)
        //{  
            this->shapes.push_back(shapes[i]); 
        
       // printf("Value",shapes[i].vertexCoordinates[0][0]);
        //}
        
        
}

GeneralizedProcrustes::~GeneralizedProcrustes()
{
        if(partialAverages != NULL)
        {
                for(int i = 0; i < numberOfShapes; i++) delete [] partialAverages[i];
                delete [] partialAverages;
                partialAverages = NULL;
        }

        shapes.clear();
}

void GeneralizedProcrustes::computeGPAlignment(vector < Shape >& outShapes)
{
        int i, counter;
        double oldEnergy, newEnergy;
        if(shapes.size() == 0 || numberOfShapes < 1) return;

        centerAllShapes();

        
        printf("Performing alignment\n");
        if(partialAverages == NULL)
        {
                partialAverages = new double * [numberOfShapes];
                for(i = 0; i < numberOfShapes; i++)
                        partialAverages[i] = new double[3*minNumNodes];
                      
        }
        
        oldEnergy = newEnergy = -1;
        counter = 0;
        do
        {
                computePartialAverages();
               // printf("Antes %f \n",shapes[0].vertexCoordinates[0][0]);
                for(i = 0; i < numberOfShapes; i++) 
                {      computeRotationAlignment(shapes[i] , partialAverages[i]);
                       outShapes[i]=shapes[i];
                }
                
                
                oldEnergy = newEnergy;
                newEnergy = evaluateEnergy();
                counter++;
        }while((oldEnergy == -1) || (fabs(newEnergy-oldEnergy) > ROTATE_THRESHOLD) && (counter < MAX_WHILE));
}


void GeneralizedProcrustes::centerAllShapes()
{
        int i, j, numberOfNodes;
        double sumX, sumY, sumZ;

        for(i = 0; i < numberOfShapes; i++)
        {
                numberOfNodes = shapes[i].numberVertices;
                sumX = 0.0; sumY = 0.0; sumZ = 0.0;
                for(j = 0; j < numberOfNodes; j++)
                {
                        sumX += shapes[i].vertexCoordinates[j][0];
                        sumY += shapes[i].vertexCoordinates[j][1];
                        sumZ += shapes[i].vertexCoordinates[j][2];
                }
                sumX = sumX/(double)numberOfNodes;
                sumY = sumY/(double)numberOfNodes;
                sumZ = sumZ/(double)numberOfNodes;

                for(j = 0; j < numberOfNodes; j++) 
                {
                        shapes[i].vertexCoordinates[j][0] -= sumX;
                        shapes[i].vertexCoordinates[j][1] -= sumY;
                        shapes[i].vertexCoordinates[j][2] -= sumZ;
                }
        }
}

void GeneralizedProcrustes::computeAverage(Shape & average)
{
        int i, j;
        double sumX, sumY, sumZ;
        if(shapes.size() == 0 || numberOfShapes < 1) return;

        for(i = 0; i < minNumNodes; i++)
        {
                sumX = 0.0; sumY = 0.0; sumZ = 0.0;
                for(j = 0; j < numberOfShapes; j++)
                {
                        sumX += shapes[j].vertexCoordinates[i][0];
                        sumY += shapes[j].vertexCoordinates[i][1];
                        sumZ += shapes[j].vertexCoordinates[i][2];
                }

                average.vertexCoordinates[i][0] = sumX/(double)numberOfShapes;
                average.vertexCoordinates[i][1] = sumY/(double)numberOfShapes;
                average.vertexCoordinates[i][2] = sumZ/(double)numberOfShapes;
        }
}

void GeneralizedProcrustes::computePartialAverages()
{
        int i, j;
        double sumX, sumY, sumZ;

        for(i = 0; i < minNumNodes; i++)
        {
                sumX = 0; sumY = 0; sumZ = 0;
                for(j = 0; j < numberOfShapes; j++)
                {
                        sumX += shapes[j].vertexCoordinates[i][0];
                        sumY += shapes[j].vertexCoordinates[i][1];
                        sumZ += shapes[j].vertexCoordinates[i][2];
                }
                for(j = 0; j < numberOfShapes; j++)
                {
                        partialAverages[j][i] = (sumX - shapes[j].vertexCoordinates[i][0])/((double)numberOfShapes-1.0);
                        partialAverages[j][minNumNodes + i] = (sumY - shapes[j].vertexCoordinates[i][1])/((double)numberOfShapes-1.0);
                        partialAverages[j][2*minNumNodes + i] = (sumZ - shapes[j].vertexCoordinates[i][2])/((double)numberOfShapes-1.0);
                }
        }
}



bool GeneralizedProcrustes::computeRotationAlignment(Shape & cont, double * partialAverage)
{
        long int i, numberOfNodes, dimension, lwork, info;
        double alpha, beta, determinant;
        char job, trans;
        bool returnval;

        // Compute the best rigid alignment based on minNumNodes first vertices:
        numberOfNodes = cont.numberVertices;
        dimension = 3;
        job = 'A';
        trans = 'N';
        alpha = 1.0;
        beta = 0.0;
        lwork = 2*minNumNodes;
        double * work = new double[lwork];
        double * A = new double[minNumNodes * dimension];
        double * b = new double[minNumNodes * dimension];
        double * S = new double[dimension];
        double * U = new double[dimension*dimension];
        double * VT = new double[dimension*dimension];
        double * R = new double[dimension*dimension];
        returnval = true;

        
       //printf("cont.vertexCoordinates[i][0]:%f\n",cont.vertexCoordinates[0][0]);
        //populate the arrays A and b:
        for(i = 0; i < minNumNodes; i++)
        {
                A[i] = cont.vertexCoordinates[i][0];
                A[minNumNodes + i] = cont.vertexCoordinates[i][1];
                A[2*minNumNodes + i] = cont.vertexCoordinates[i][2];
                b[i] = partialAverage[i];
                b[minNumNodes + i] = partialAverage[minNumNodes + i];
                b[2*minNumNodes + i] = partialAverage[2*minNumNodes + i];
        }
        //compute an estimate of R:
        dgels_(&trans, &minNumNodes, &dimension, &dimension, A, &minNumNodes, b, &minNumNodes, work, &lwork, &info);
        if(info != 0)
        {
                printf("Problem with estimating R %s\n", info);
                goto align_EXIT;
        }
        //make sure that R is a valid matrix (orthonormal): b contains R
        delete [] work;
        lwork = 5*dimension;
        work = new double[lwork];
        //only copy 3 by 3 submatrix to R:
        for(i = 0; i < dimension; i++)
        {
                R[i] = b[i];
                R[dimension + i] = b[minNumNodes + i];
                R[2*dimension + i] = b[2*minNumNodes + i];
        }
        dgesvd_(&job, &job, &dimension, &dimension, R, &dimension, S, U, &dimension, VT, &dimension, work, 
                &lwork, &info);
        if(info != 0)
        {
                printf("Problem with estimating R in SVD %s\n",info);
                goto align_EXIT;
        }
        //set S to I: multiply U times VT:
       dgemm_(&trans, &trans, &dimension, &dimension, &dimension, &alpha, U, &dimension, VT, &dimension, &beta, 
                R, &dimension);
        //disallow reflections:
        determinant = R[0]*R[4]*R[8]+R[1]*R[5]*R[6]+R[2]*R[3]*R[7]-R[2]*R[4]*R[6]-R[0]*R[5]*R[7]-R[1]*R[3]*R[8];
        if(determinant < 0) 
        {
                printf("Determinant is %d\n", determinant);
                returnval = false;
                goto align_EXIT;
        }
        
        //transform ALL the coordinates:
        lwork = 2*numberOfNodes;
        delete [] work;
        work = new double[lwork];
        delete [] A;
        A = new double[numberOfNodes * dimension];
        delete [] b;
        b = new double[numberOfNodes * dimension];
        //populate the array A again as it got destroyed:
        for(i = 0; i < numberOfNodes; i++)
        {
                A[i] = cont.vertexCoordinates[i][0];
                A[numberOfNodes + i] = cont.vertexCoordinates[i][1];
                A[2*numberOfNodes + i] = cont.vertexCoordinates[i][2];
        }
        //do the rigid transformation and set it to shape:
        dgemm_(&trans, &trans, &numberOfNodes, &dimension, &dimension, &alpha, A, &numberOfNodes, R, &dimension, &beta,
                b, &numberOfNodes);
        for(i = 0; i < numberOfNodes; i++)
        {
            
               // printf("cont.vertexCoordinates[i][0]:%f  b[i] %f\n", cont.vertexCoordinates[i][0] ,b[i]);
                cont.vertexCoordinates[i][0] = b[i];
                cont.vertexCoordinates[i][1] = b[numberOfNodes + i];
                cont.vertexCoordinates[i][2] = b[2*numberOfNodes + i];
        }

align_EXIT:
        delete [] work;
        delete [] A;
        delete [] b;
        delete [] S;
        delete [] U;
        delete [] VT;
        delete [] R;

        return returnval;
}

double GeneralizedProcrustes::evaluateEnergy()
{
        int i, j, k;
        double energy, sumX, sumY, sumZ;
        energy = 0.0;

        for(i = 0; i < numberOfShapes; i++)
        {
                for(j = i+1; j < numberOfShapes; j++)
                {
                        sumX = 0.0; sumY = 0.0; sumZ = 0.0;
                        for(k = 0; k < minNumNodes; k++)
                        {
                                sumX += pow(shapes[i].vertexCoordinates[k][0]-shapes[j].vertexCoordinates[k][0], 2);
                                sumY += pow(shapes[i].vertexCoordinates[k][1]-shapes[j].vertexCoordinates[k][1], 2);
                                sumZ += pow(shapes[i].vertexCoordinates[k][2]-shapes[j].vertexCoordinates[k][2], 2);
                        }
                        energy += sqrt(sumX+sumY+sumZ);
                }
        }

        return energy;
}


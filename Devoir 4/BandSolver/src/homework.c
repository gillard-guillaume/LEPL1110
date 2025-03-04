

#include"fem.h"


#ifndef NORENUMBER 

void femMeshRenumber(femMesh *theMesh, femRenumType renumType)
{
    int i;
    
    switch (renumType) {
        case FEM_NO :
            for (i = 0; i < theMesh->nodes->nNodes; i++) 
                theMesh->nodes->number[i] = i;
            break;
// 
// A modifier :-)
// debut
//
        case FEM_XNUM : 
        case FEM_YNUM : 
            for (i = 0; i < theMesh->nodes->nNodes; i++) 
                theMesh->nodes->number[i] = i;
            break;            
// 
// end
//

        default : Error("Unexpected renumbering option"); }
}

#endif
#ifndef NOBAND 

int femMeshComputeBand(femMesh *theMesh)
{
    int myBand = theMesh->nodes->nNodes;
    return(myBand);
}


#endif
#ifndef NOBANDASSEMBLE


void femBandSystemAssemble(femBandSystem* myBandSystem, double *Aloc, double *Bloc, int *map, int nLoc)
{
    double *A = myBandSystem->A;
    double *B = myBandSystem->B;
    int band = myBandSystem->band;
    int i, j, iGlob, jGlob;

    for (i = 0; i < nLoc; i++) {
        iGlob = map[i];
    
        for (j = 0; j < nLoc; j++) {
            jGlob = map[j];
    
            if (abs(iGlob - jGlob) <= band) {
                int index = (iGlob * (band + 1)) + (jGlob - iGlob + band);
                A[index] += Aloc[i*nLoc+j];
            }
        }
    }

    return;
}


#endif
#ifndef NOBANDELIMINATE


double  *femBandSystemEliminate(femBandSystem *myBand)
{
    double  *A, *B, factor;
    int     i, j, k, jend, size, band;
    A    = myBand->A;
    B    = myBand->B;
    size = myBand->size;
    band = myBand->band;
    
    // Triangularisation of the matrix
    for (i = 0; i < size; i++) {
        double pivot = A[i * (band + 1) + band]; 

        for (j = i + 1; j <= i + band && j < size; j++) {
            factor = A[j * (band + 1) + (i - j + band)] / pivot;
            for (k = i; k <= i + band && k < size; k++) {
                A[j * (band + 1) + (i - j + band)] -= factor * A[i * (band + 1) + (i - k + band)];
            }
            B[j] -= factor * B[i];
        }

    }

    // Back-substitution
    for (i = size - 1; i >= 0; i--) {
        jend = i + band;
        if (jend >= size) jend = size - 1;
        B[i] /= A[i * (band + 1) + band];
        for (j = i + 1; j <= jend; j++) {
            B[i] -= A[j * (band + 1) + (i - j + band)] * B[j] / A[i * (band + 1) + band];
        }
    }
    


    return(myBand->B);
}


#endif


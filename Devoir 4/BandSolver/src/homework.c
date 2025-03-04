

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
    double  **A, *B, factor;
    int     i, j, k, jend, size, band;
    A    = myBand->A;
    B    = myBand->B;
    size = myBand->size;
    band = myBand->band;
    
    // A completer :-)


    return(myBand->B);
}


#endif


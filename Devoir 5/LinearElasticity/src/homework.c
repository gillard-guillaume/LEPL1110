#include "fem.h"




void geoMeshGenerate() {

    femGeo* theGeometry = geoGetGeometry();

    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;
    
    int ierr;
    double r = w/4;
    int idRect = gmshModelOccAddRectangle(0.0,0.0,0.0,w,h,-1,0.0,&ierr); 
    int idDisk = gmshModelOccAddDisk(w/2.0,h/2.0,0.0,r,r,-1,NULL,0,NULL,0,&ierr); 
    int idSlit = gmshModelOccAddRectangle(w/2.0,h/2.0-r,0.0,w,2.0*r,-1,0.0,&ierr); 
    int rect[] = {2,idRect};
    int disk[] = {2,idDisk};
    int slit[] = {2,idSlit};

    gmshModelOccCut(rect,2,disk,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr); 
    gmshModelOccCut(rect,2,slit,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr); 
    gmshModelOccSynchronize(&ierr); 

    if (theGeometry->elementType == FEM_QUAD) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshOptionSetNumber("Mesh.RecombineAll",1,&ierr);
        gmshOptionSetNumber("Mesh.Algorithm",11,&ierr);  
        gmshOptionSetNumber("Mesh.SmoothRatio", 21.5, &ierr);  
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm",1.0,&ierr); 
        gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  
        gmshModelMeshGenerate(2,&ierr);  }
  
    if (theGeometry->elementType == FEM_TRIANGLE) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshModelMeshGenerate(2,&ierr);  }
 
    return;
}


double *femElasticitySolve(femProblem *theProblem)
{

    femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete    *theSpace = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theMesh = theGeometry->theElements;
    
    
    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    int iElem,iInteg,iEdge,i,j,d,map[4],mapX[4],mapY[4];
    
    int nLocal = theMesh->nLocalNode;

    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;      
    double rho = theProblem->rho;
    double g   = theProblem->g;
    double **A = theSystem->A;
    double *B  = theSystem->B;
    
    
    // Collecting the elements (triangles or quads)
    for (iElem = 0; iElem < theMesh->nElem; iElem++){

        // Collecting the global indices of the nodes of the element and stocking them temporarily in map
        for (i = 0; i < nLocal; i++) {

            map[i] = theMesh->elem[iElem*nLocal + i]; 

            // Collecting the coordinates of the nodes
            x[i] = theNodes->X[map[i]];
            y[i] = theNodes->Y[map[i]]; 

            // Duplicating the nodes for the displacement in x and y
            mapX[i] = 2 * map[i];
            mapY[i] = 2 * map[i] + 1; 
        }
        
        // Loop over the integration points

        for (iInteg = 0; iInteg < theRule->n; iInteg++){
            // Getting the coordinates and the weights of the integration points
            double xsi = theRule->xsi[iInteg];
            double eta = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];

            // Computing the shape functions and their derivatives
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);

            // Computing Jacobian matrix
            double dxdxsi = 0.0, dydxsi = 0.0, dxdeta = 0.0, dydeta = 0.0;
            for (i = 0; i < theSpace->n; i++) {
                dxdxsi += dphidxsi[i] * x[i];
                dydxsi += dphidxsi[i] * y[i];
                dxdeta += dphideta[i] * x[i];
                dydeta += dphideta[i] * y[i]; 
            }
            double jacobian = fabs(dxdxsi * dydeta - dydxsi * dxdeta);

            // Transformation of the derivatives of the shape functions to global coordinates
            for (i = 0; i < theSpace->n; i++) {
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jacobian;
                dphidy[i] = (-dphidxsi[i] * dxdeta + dphideta[i] * dxdxsi) / jacobian; 
            }

            // Computing weighted Jacobian
            double weightJacobian = weight * jacobian;

            // Assembling the rigidity matrix A
            for (i = 0; i < theSpace->n; i++){
                for(j = 0; j < theSpace->n; j++){
                    A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + dphidy[i] * c * dphidy[j]) * weightJacobian;
                    A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + dphidy[i] * c * dphidx[j]) * weightJacobian;
                    A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + dphidx[i] * c * dphidy[j]) * weightJacobian;
                    A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + dphidx[i] * c * dphidx[j]) * weightJacobian;
                }
                // Assembling the load vector B
                B[mapY[i]] -= g * rho * phi[i] * weightJacobian;
            }
        }

    }
        
                
                
  
    int *theConstrainedNodes = theProblem->constrainedNodes;     
    for (int i=0; i < theSystem->size; i++) {
        if (theConstrainedNodes[i] != -1) {
            double value = theProblem->conditions[theConstrainedNodes[i]]->value;
            femFullSystemConstrain(theSystem,i,value); }}
                            
    return femFullSystemEliminate(theSystem);
}
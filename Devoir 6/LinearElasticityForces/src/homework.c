#include "fem.h"

void femElasticityAssembleElements(femProblem *theProblem){
    femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete    *theSpace = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theMesh = theGeometry->theElements;
    femMesh        *theEdges = theGeometry->theEdges;
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
    
    
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j=0; j < nLocal; j++) {
            map[j]  = theMesh->elem[iElem*nLocal+j];
            mapX[j] = 2*map[j];
            mapY[j] = 2*map[j] + 1;
            x[j]    = theNodes->X[map[j]];
            y[j]    = theNodes->Y[map[j]];} 
        
        for (iInteg=0; iInteg < theRule->n; iInteg++) {    
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];  
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            
            double dxdxsi = 0.0;
            double dxdeta = 0.0;
            double dydxsi = 0.0; 
            double dydeta = 0.0;
            for (i = 0; i < theSpace->n; i++) {  
                dxdxsi += x[i]*dphidxsi[i];       
                dxdeta += x[i]*dphideta[i];   
                dydxsi += y[i]*dphidxsi[i];   
                dydeta += y[i]*dphideta[i]; }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            
            for (i = 0; i < theSpace->n; i++) {    
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; }            
            for (i = 0; i < theSpace->n; i++) { 
                for(j = 0; j < theSpace->n; j++) {
                    A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + 
                                            dphidy[i] * c * dphidy[j]) * jac * weight;                                                                                            
                    A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + 
                                            dphidy[i] * c * dphidx[j]) * jac * weight;                                                                                           
                    A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + 
                                            dphidx[i] * c * dphidy[j]) * jac * weight;                                                                                            
                    A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + 
                                            dphidx[i] * c * dphidx[j]) * jac * weight; }}
             for (i = 0; i < theSpace->n; i++) {
                B[mapY[i]] -= phi[i] * g * rho * jac * weight; }}} 
}


void femElasticityAssembleNeumann(femProblem *theProblem){
    femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->ruleEdge;
    femDiscrete    *theSpace = theProblem->spaceEdge;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theEdges = theGeometry->theEdges;
    double x[2],y[2],phi[2];
    int iBnd,iElem,iInteg,iEdge,i,j,d,map[2],mapU[2];
    int nLocal = 2;
    double *B  = theSystem->B;

    for(iBnd=0; iBnd < theProblem->nBoundaryConditions; iBnd++){
        femBoundaryCondition *theCondition = theProblem->conditions[iBnd];
        femBoundaryType type = theCondition->type;
        femDomain *theDomain = theCondition->domain;
        double value = theCondition->value;

        // Checking if the boundary condition is a Neumann condition
        if (type == NEUMANN_X || type == NEUMANN_Y){
            // Looping on the domain of the boundary condition
            for (iEdge = 0; iEdge < theDomain->nElem; iEdge++){
                // Getting the global element number
                iElem = theDomain->elem[iEdge];

                for (i = 0; i < nLocal; i++){
                    map[i] = theEdges->elem[iElem*nLocal + i];
                    mapU[i] = (type == NEUMANN_X) ? 2 * map[i] : 2 * map[i] + 1;
                    x[i] = theNodes->X[map[i]];
                    y[i] = theNodes->Y[map[i]];
                }
                
                // Computing the 1D jacobian
                double dx = x[1] - x[0];
                double dy = y[1] - y[0];
                double jacobian = sqrt(dx*dx + dy*dy)/2;

                // Integrating in 1D
                for (iInteg = 0; iInteg < theRule->n; i++){
                    double xsi = theRule->xsi[iInteg];
                    double weight = theRule->weight[iInteg];

                    femDiscretePhi(theSpace, xsi, phi);

                    // Assembling Forces in B
                    for (i = 0; i < theSpace->n; i++){
                        B[mapU[i]] += phi[i] * value * jacobian * weight;
                    }
                    
                }
            }
        }

    }
}


double **A_copy = NULL;
double *B_copy  = NULL;
int size = 0;


double *femElasticitySolve(femProblem *theProblem){
    femFullSystem *theSystem = theProblem->system;
    double **A = theSystem->A;
    double *B = theSystem->B;
    double *U = theProblem->soluce;

    // Assembling the system
    femElasticityAssembleElements(theProblem);
    femElasticityAssembleNeumann(theProblem);

    size = theSystem->size;

    // Copying the system to not
    if (A_copy == NULL)
    {
        A_copy = (double **) malloc(sizeof(double *) * size);
        for (int i = 0; i < size; i++) { A_copy[i] = (double *) malloc(sizeof(double) * size); }
    }
    if (B_copy == NULL) { B_copy = (double *) malloc(sizeof(double) * size); }

    // Copy the stiffness matrix A and the load vector B
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++) { A_copy[i][j] = theSystem->A[i][j]; }
        B_copy[i] = theSystem->B[i];
    }


    // Applying Dirichlet boundary conditions
    int size = theSystem->size;
    int *contrainedNodes = theProblem->constrainedNodes;
    for (int i = 0; i < size; i++){
        if (contrainedNodes[i] != -1){
            double value = theProblem->conditions[contrainedNodes[i]]->value;
            femFullSystemConstrain(theSystem, i, value);
        }
    }

    // Solving the system
    U = femFullSystemEliminate(theSystem);
    memcpy(theProblem->soluce, theSystem->B, theSystem->size * sizeof(double));

    return theProblem->soluce;
}

double * femElasticityForces(femProblem *theProblem){        
    femFullSystem *theSystem = theProblem->system;
    double **A = theSystem->A;
    double *B = theSystem->B;
    double *U = theProblem->soluce;
    double *R = theProblem->residuals;
    int nNode = theProblem->geometry->theNodes->nNodes;

    size = theSystem->size;

    if (R == NULL) { R = (double *) malloc(sizeof(double) * size); }
    
    for (int i = 0; i < size; i++) { R[i] = 0.0; }

    // Computing the residuals R = B - A*U
    // If a node is not contrained, the equation A*U = B is satisfied and the residual is 0
    // If a node is contrained, the equation A*U = B is not satisfied and the residual give the reaction force
    for(int i = 0; i < size; i++){

        for(int j = 0; j < size; j++){
            R[i] -= A_copy[i][j] * U[j];
        }
        R[i] -= B_copy[i];

    }

    for (int i = 0; i < size; i++) { free(A_copy[i]); A_copy[i] = NULL;}
    free(A_copy); free(B_copy);
    A_copy = NULL; B_copy = NULL;

    return R;
}
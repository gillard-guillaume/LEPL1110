#include "fem.h"


double hermite_interpolation(double x, double h0, double hstar, double dstar) {
    // f(x) = a*x^3 + b*x^2 + c*x + d
    double a = 2 * (h0 - hstar) / pow(dstar, 3);
    double b = -3 * (h0 - hstar) / pow(dstar, 2);
    double c = 0;
    double d = h0;
    return a * pow(x, 3) + b * pow(x, 2) + c * x + d;
}

double geoSize(double x, double y){

    femGeo* theGeometry = geoGetGeometry();
    
    double h = theGeometry->h;
    double x0 = theGeometry->xNotch;
    double y0 = theGeometry->yNotch;
    double r0 = theGeometry->rNotch;
    double h0 = theGeometry->hNotch;
    double d0 = theGeometry->dNotch;
  
    
    double x1 = theGeometry->xHole;
    double y1 = theGeometry->yHole;
    double r1 = theGeometry->rHole;
    double h1 = theGeometry->hHole;
    double d1 = theGeometry->dHole;

    // Computing the distance to the notch and the hole
    double dNotch = sqrt(pow(x - x0, 2) + pow(y - y0, 2)) - r0;
    double dHole  = sqrt(pow(x - x1, 2) + pow(y - y1, 2)) - r1;

    // Notch interpolation 
    if (dNotch <= d0) {
        h = hermite_interpolation(dNotch, h0, h, d0);
    }

    // Hole interpolation
    if (dHole <= d1) {
        h = hermite_interpolation(dHole, h1, h, d1);
    }
        
    
     
    return h;
    
//   
// Your contribution ends here :-)
//

}


#define ___ 0

void geoMeshGenerate() {

    femGeo* theGeometry = geoGetGeometry();

    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;
     
    double x0 = theGeometry->xNotch;
    double y0 = theGeometry->yNotch;
    double r0 = theGeometry->rNotch;
    
    
    double x1 = theGeometry->xHole;
    double y1 = theGeometry->yHole;
    double r1 = theGeometry->rHole;
 
//
//  -1- Construction de la géométrie avec OpenCascade
//      On crée le rectangle
//      On crée les deux cercles
//      On soustrait les cercles du rectangle :-)
//
 
    int ierr;
    int idPlate = gmshModelOccAddRectangle(x0, y0, 0, w, h, -1, 0, &ierr);   
    ErrorGmsh(ierr);
    int idNotch = gmshModelOccAddDisk(x0, y0, 0, r0, r0, -1, NULL, 0, NULL, 0, &ierr); 
    ErrorGmsh(ierr);
    int idHole  = gmshModelOccAddDisk(x1, y1, 0, r1, r1, -1, NULL, 0, NULL, 0, &ierr);    
    ErrorGmsh(ierr);
    
    int plate[] = {2, idPlate};
    int notch[] = {2, idNotch};
    int hole[]  = {2, idHole};
    gmshModelOccCut(plate, 2, notch, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr); 
    ErrorGmsh(ierr);
    gmshModelOccCut(plate, 2, hole, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr); 
    ErrorGmsh(ierr);
 
//
//  -2- Définition de la fonction callback pour la taille de référence
//      Synchronisation de OpenCascade avec gmsh
//      Génération du maillage (avec l'option Mesh.SaveAll :-)
                  
   
    geoSetSizeCallback(geoSize);
                                  
    gmshModelOccSynchronize(&ierr);       
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr);  
       
//
//  Generation de quads :-)
//
//    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
//    gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
//    gmshOptionSetNumber("Mesh.Algorithm", 8, &ierr);  chk(ierr);
//    gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr);  chk(ierr);
//    gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  chk(ierr);
//    gmshModelMeshGenerate(2, &ierr);  
   
 
//
//  Plot of Fltk
//
//   gmshFltkInitialize(&ierr);
//   gmshFltkRun(&ierr);  chk(ierr);
//
    
}
#include <stdio.h>
#include <math.h>
#include "glfem.h"


// Guillaume Gillard
// Igor Gregoire


double integrate(double x[3], double y[3], double (*f) (double, double))
{
    double I = 0;
    double xLoc[3];
    double yLoc[3];

    double xi[3] = {1/6, 1/6, 2/3};
    double eta[3] = {1/6, 2/3, 1/6};
    double w[3] = {1/6, 1/6, 1/6};


    for (int i = 0; i < 3; i++)
    // applying variable change on xi et eta and then the integration formula
    {
        xLoc[i] = x[0] * (1 - xi[i] - eta[i]) + x[1] * xi[i] + x[2] * eta[i];
        yLoc[i] = y[0] * (1 - xi[i] - eta[i]) + y[1] * xi[i] + y[2] * eta[i];
        I += w[i] * f(xLoc[i], yLoc[i]);
    }




    glfemSetColor(GLFEM_BLACK); glfemDrawElement(x,y,3);
    glfemSetColor(GLFEM_BLUE);  glfemDrawNodes(x,y,3);
    glfemSetColor(GLFEM_RED);   glfemDrawNodes(xLoc,yLoc,3);
    


    return I;
}

double integrateRecursive(double x[3], double y[3], double (*f)(double,double), int n)
{

    double I = 0;

    // decomposing into subtriangles by placing points on the middle of corners
    if (n > 0)
    {
        double newpoints_x[3] = { (x[0] + x[1]) / 2, (x[1] + x[2]) / 2, (x[2] + x[0]) / 2 };
        double newpoints_y[3] = { (y[0] + y[1]) / 2, (y[1] + y[2]) / 2, (y[2] + y[0]) / 2 };

        triangle1_x = {x[0], newpoints_x[0], newpoints_x[2]};
        triangle1_y = {y[0], newpoints_y[0], newpoints_y[2]};

        triangle2_x = {newpoints_x[0], x[1], newpoints_x[1]};
        triangle2_y = {newpoints_y[0], y[1], newpoints_y[1]};

        triangle3_x = {newpoints_x[2], newpoints_x[1], x[2]};
        triangle3_y = {newpoints_y[2], newpoints_y[1], y[2]};

        triangle4_x = {newpoints_x[0], newpoints_x[1], newpoints_x[2]};
        triangle4_y = {newpoints_y[0], newpoints_y[1], newpoints_y[2]};

        I += integrateRecursive(triangle1_x, triangle1_y, f, n-1);
        I += integrateRecursive(triangle2_x, triangle2_y, f, n-1);
        I += integrateRecursive(triangle3_x, triangle3_y, f, n-1); 
        I += integrateRecursive(triangle4_x, triangle4_y, f, n-1);
    }
    
    else if (n == 0)
    {
        I += integrate(x,y,f);
    }

    return I;
}

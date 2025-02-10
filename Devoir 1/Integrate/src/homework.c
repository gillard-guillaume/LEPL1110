#include <stdio.h>
#include <math.h>
#include "glfem.h"


// Guillaume Gillard
// Igor Gregoire


double integrate(double x[3], double y[3], double (*f) (double, double))
{
    double I = 0.0;
    double xLoc[3];
    double yLoc[3];

    double xi[3] = {1.0/6, 1.0/6, 2.0/3};
    double eta[3] = {1.0/6, 2.0/3, 1.0/6};
    double w[3] = {1.0/6, 1.0/6, 1.0/6};

    double J = fabs((x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0]));

    for (int i = 0; i < 3; i++)
    // applying variable change on xi et eta
    {
        xLoc[i] = x[0] * (1 - xi[i] - eta[i]) + x[1] * xi[i] + x[2] * eta[i];
        yLoc[i] = y[0] * (1 - xi[i] - eta[i]) + y[1] * xi[i] + y[2] * eta[i];
    }

    for (int i = 0; i < 3; i++)
    {
        I += J * w[i] * f(xLoc[i], yLoc[i]);
    }


    glfemSetColor(GLFEM_BLACK); glfemDrawElement(x,y,3);
    glfemSetColor(GLFEM_BLUE);  glfemDrawNodes(x,y,3);
    glfemSetColor(GLFEM_RED);   glfemDrawNodes(xLoc,yLoc,3);
    
    return I;
}

double integrateRecursive(double x[3], double y[3], double (*f)(double,double), int n)
{

    double I = 0.0;

    if (n == 1) return integrate(x, y, f);

    // decomposing into subtriangles by placing points on the middle of corners

    double newpoints_x[3] = { (x[0] + x[1]) / 2, (x[1] + x[2]) / 2, (x[2] + x[0]) / 2 };
    double newpoints_y[3] = { (y[0] + y[1]) / 2, (y[1] + y[2]) / 2, (y[2] + y[0]) / 2 };

    double triangle1_x[3] = {x[0], newpoints_x[0], newpoints_x[2]};
    double triangle1_y[3] = {y[0], newpoints_y[0], newpoints_y[2]};

    double triangle2_x[3] = {newpoints_x[0], x[1], newpoints_x[1]};
    double triangle2_y[3] = {newpoints_y[0], y[1], newpoints_y[1]};

    double triangle3_x[3] = {newpoints_x[2], newpoints_x[1], x[2]};
    double triangle3_y[3] = {newpoints_y[2], newpoints_y[1], y[2]};

    double triangle4_x[3] = {newpoints_x[0], newpoints_x[1], newpoints_x[2]};
    double triangle4_y[3] = {newpoints_y[0], newpoints_y[1], newpoints_y[2]};

    I += integrateRecursive(triangle1_x, triangle1_y, f, n-1);
    I += integrateRecursive(triangle2_x, triangle2_y, f, n-1);
    I += integrateRecursive(triangle3_x, triangle3_y, f, n-1); 
    I += integrateRecursive(triangle4_x, triangle4_y, f, n-1);

    return I;
}

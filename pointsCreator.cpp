#include <stdio.h>
#include <iostream>
#include <fstream>
#include <ctime>

using namespace std;

int ALL_POINTS = 500;

double sqr(double x)
{
    return x*x;
}

bool isInsideEllipse(double x, double y, double a, double b, double x0, double y0)
{
    if (sqr(x-x0)/sqr(a)+sqr(y-y0)/sqr(b) <= 1)
        return true;
    return false;
}

int LEFT = -200;
int RIGHT = 300;

int main()
{   
    srand (time (0));
    ofstream pointsfile;
    pointsfile.open ("input.txt");

    int count = 0;
    double x, y;
    //200*((double) rand() / RAND_MAX)-100

    while (count < ALL_POINTS)
    {
        x = (RIGHT-LEFT)*((double) rand() / RAND_MAX)-RIGHT;
        y = (RIGHT-LEFT)*((double) rand() / RAND_MAX)-RIGHT;

        if (isInsideEllipse(x, y, 13, 77, -67, -13))
        {
            pointsfile << x << " " << y << " 0\n";
            cout << x << " " << y << endl;
            count++;
        }
        else if (isInsideEllipse(x, y, 32, 40, 40, -50))
        {
            pointsfile << x << " " << y << " 1\n";
            cout << x << " " << y << endl;
            count++;
        }
        else if (isInsideEllipse(x, y, 34, 40, 40, 50))
        {
            pointsfile << x << " " << y << " 2\n";
            cout << x << " " << y << endl;
            count++;
        }
    }
    
    return 0;
};
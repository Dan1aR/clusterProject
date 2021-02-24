#include <stdio.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <math.h>
#include <array>

//set pointsize 3
//plot "< awk '{if($3 == \"0\") print}' input.txt" u 1:2 t "red" w p pt 2, "< awk '{if($3 == \"1\") print}' input.txt" u 1:2 t "green" w p pt 2, "< awk '{if($3 == \"2\") print}' input.txt" u 1:2 t "blue" w p pt 2

//plot "< awk '{if($3 == \"0\") print}' output.txt" u 1:2 t "red" w p pt 2, "< awk '{if($3 == \"1\") print}' output.txt" u 1:2 t "green" w p pt 2, "< awk '{if($3 == \"2\") print}' output.txt" u 1:2 t "blue" w p pt 2, "< awk '{if($3 == \"3\") print}' output.txt" u 1:2 t "red" w p pt 2

using namespace std;

int K = 3;
int INTERVAL = 50;
string FILENAME = "input.txt";

//--------------------Functions------------------------------------------
double sqr(double x)
{
    return x * x;
}
//-----------------------------------------------------------------------


class Point
{
    public:
        double x;
        double y;
        int nearestCentroidNum;

        void Init(double px, double py)
        {
            x = px;
            y = py;
        }
};

class SetPoints
{
    public:
        vector<Point> points;

        Point mostLeft;
        Point mostRight;
        Point mostBottom;
        Point mostHigh;

        void defineSpecPoint(double x, double y)
        {
            if (x < mostLeft.x)
            {
                mostLeft.Init(x, y);
            }
            if (x > mostRight.x)
            {
                mostRight.Init(x, y);
            }
            if (y < mostBottom.y)
            {
                mostBottom.Init(x, y);
            }
            if (y > mostHigh.y)
            {
                mostHigh.Init(x, y);
            }
        }

        void pointsInit(string fileName)
        {
            double x;
            double y;
            int clasterNumber;

            ifstream pointsfile;
            pointsfile.open(fileName);

            while (pointsfile >> x >> y >> clasterNumber)
            {
                Point point;
                point.Init(x, y);
                points.push_back(point);
                
                defineSpecPoint(x, y);
            }

            pointsfile.close();
        }
};

class Clusters
{
    private:
        double Rmin(SetPoints U, SetPoints V)
        {
            double minR = 1e10;
            for (auto& u : U.points)
            {
                for (auto& v : V.points)
                {
                    double R = sqrt(sqr(u.x - v.x) + sqr(u.y - v.y));
                    if (R < minR)
                    {
                        minR = R;
                    }
                }
            }
            return minR;
        }
    public:
        vector<SetPoints> clusters;

        void pointsToClusters(SetPoints *field)
        {
            for (auto& point : field->points)
            {
                SetPoints newCluster;
                newCluster.points.push_back(point);
                clusters.push_back(newCluster);
            }
        }

        array<double, 3> getNearestClustersIdxs()
        {
            double minR = 1e10;
            int minI, minJ;

            for (int i = 0; i < clusters.size(); ++i)
            {
                for (int j = 0; j < clusters.size(); ++j)
                {
                    double R = Rmin(clusters[i], clusters[j]);
                    if ((i != j) && (R < minR))
                    {
                        minR = R;
                        minI = i;
                        minJ = j;
                    }
                }
            }
            array<double, 3> IJR = { (double)minI, (double)minJ, minR};
            cout << IJR[0] << " " << IJR[1] << " " << IJR[2] << endl;
            return IJR;
        }

        void clustersUnion(int i, int j)
        {
            for (auto& point : clusters[ j ].points)
            {
                clusters[ i ].points.push_back(point);
            }
            clusters.erase(clusters.begin() + j );
        }
};

class Algorithms
{
    private:
        SetPoints workGridField;
        SetPoints workField;
        vector<Point> centroids;

        void printResult(SetPoints *pworkField)
        {
            ofstream pointsfile;
            pointsfile.open("output.txt");
            //int i = 0;
            for (auto& point : pworkField->points)
            {
                
                if ( (point.x != 0) && (point.y != 0) )
                {
                    //cout << point.x << " " << point.y << " " << point.nearestCentroidNum << "\n";
                    pointsfile << point.x << " " << point.y << " " << point.nearestCentroidNum << "\n";
                    //i++;
                }
            }
            //cout << i << endl;
        }
        
        //-----------------------------------------------------------------------

        SetPoints redefineField(int interval)
        {
            SetPoints gridField;

            double deltaX = (workField.mostRight.x - workField.mostLeft.x) / interval;
            double deltaY = (workField.mostHigh.y - workField.mostBottom.y) / interval;

            double y0 = workField.mostBottom.y;
            double y1 = y0 + deltaY;
            
            vector<Point> newPoints;

            for (int i = 0; i < interval; ++ i)
            {
                double x0 = workField.mostLeft.x;
                double x1 = x0 + deltaX;

                for (int j = 0; j < interval; ++j)
                {
                    double cx = 0;
                    double cy = 0;
                    double cn = 0;

                    for (auto& point : workField.points)
                    {
                        if ((point.x > x0) && (point.x < x1) && (point.y > y0) && (point.y < y1))
                        {
                            cx += point.x;
                            cy += point.y;
                            cn++;
                        }
                    }

                    Point newPoint;
                    newPoint.Init(0, 0);
                    newPoint.nearestCentroidNum = 0;
                    if (cn > 0)
                    {
                        newPoint.x = cx/cn; 
                        newPoint.y = cy/cn;   
                    }
                    newPoints.push_back(newPoint);
                    x0 = x1;
                    x1 += deltaX;
                }
                
                y0 = y1;
                y1 += deltaY;
            }

            //cout << newPoints.size() << endl;

            gridField.points = newPoints;
            for (auto& point : gridField.points)
            {
                gridField.defineSpecPoint(point.x , point.y);
            }

            return gridField;
        }
        
        void recConnectCluster(SetPoints *field, int interval, int myi, int myj, int mycn)
        {
            if ((myi >= 0) && (myj >= 0) && (myi < interval) && (myj < interval))
            {
                if ( (field->points[interval*(myi)+ myj].x != 0) && (field->points[interval*(myi)+ myj].nearestCentroidNum == 0) )
                {
                    //cout << mycn << ") " << myi << " " << myj << endl;
                    field->points[interval*(myi)+ myj].nearestCentroidNum = mycn;

                    recConnectCluster(field, interval, myi+1, myj, mycn);
                    recConnectCluster(field, interval, myi-1, myj, mycn);
                    recConnectCluster(field, interval, myi, myj+1, mycn);
                    recConnectCluster(field, interval, myi, myj-1, mycn);

                    recConnectCluster(field, interval, myi+1, myj+1, mycn);
                    recConnectCluster(field, interval, myi+1, myj-1, mycn);
                    recConnectCluster(field, interval, myi-1, myj+1, mycn);
                    recConnectCluster(field, interval, myi-1, myj-1, mycn);
                    
                }
            }
        }

        void connectComponents(SetPoints *field, int interval)
        {
            int cn = 1;
            //int grcn = 0;
            for (int i = 0; i < interval; ++i)
            {
                for (int j = 0; j < interval; ++j)
                {
                    if (field->points[interval*(i)+ j].x != 0)
                    {
                        if (field->points[interval*(i)+ j].nearestCentroidNum == 0)
                        {
                            
                            field->points[interval*(i)+ j].nearestCentroidNum = cn;
                            
                            recConnectCluster(field, interval, i+1, j, cn);
                            recConnectCluster(field, interval, i-1, j, cn);
                            recConnectCluster(field, interval, i, j+1, cn);
                            recConnectCluster(field, interval, i, j-1, cn);

                            recConnectCluster(field, interval, i+1, j+1, cn);
                            recConnectCluster(field, interval, i+1, j-1, cn);
                            recConnectCluster(field, interval, i-1, j+1, cn);
                            recConnectCluster(field, interval, i-1, j-1, cn);

                            cn++;
                        }
                        //grcn++;
                    }
                }
            }
            //cout << grcn << endl;
            
        }

        void gridToWorkField(int interval)
        {
            double deltaX = (workField.mostRight.x - workField.mostLeft.x) / interval;
            double deltaY = (workField.mostHigh.y - workField.mostBottom.y) / interval;

            double y0 = workField.mostBottom.y;
            double y1 = y0 + deltaY;

            for (int i = 0; i < interval; ++ i)
            {
                double x0 = workField.mostLeft.x;
                double x1 = x0 + deltaX;

                for (int j = 0; j < interval; ++j)
                {
                    for (auto& point : workField.points)
                    {
                        if ((point.x > x0) && (point.x < x1) && (point.y > y0) && (point.y < y1))
                        {
                            point.nearestCentroidNum = workGridField.points[interval*i + j].nearestCentroidNum;
                        }
                    }

                    x0 = x1;
                    x1 += deltaX;
                }
                
                y0 = y1;
                y1 += deltaY;
            }


        }
        //--------------------------------------------------------------------------------------------------------------

        void defineKCentroids(int k)
        {
            srand (time (0));
            for (int i = 0; i <= k; ++i)
            {
                double x = (workField.mostRight.x - workField.mostLeft.x)*((double) rand() / RAND_MAX) - workField.mostRight.x;
                double y = (workField.mostHigh.y - workField.mostBottom.y)*((double) rand() / RAND_MAX) - workField.mostHigh.y;
                Point centroid;
                centroid.Init(x, y);
                centroids.push_back(centroid);
                cout << centroids[i].x << " " << centroids[i].y << endl;
            }
        }

        int findNearestCentroid(int k, double x, double y)
        {
            double lMin = 1e10;
            double l;
            int nearestCentoidNum;


            for (int i = 1; i <= k; ++i)
            {
                l = sqrt(sqr(x - centroids[i].x) + sqr(y - centroids[i].y));
                //cout << "l" << i << " : " << l << endl;
                if (l < lMin)
                {
                    lMin = l;
                    nearestCentoidNum = i;
                }
            }

            return nearestCentoidNum;
        }

        Point resetCenroid(int k, int centroidNum)
        {
            //Some bad code here
            double cx = 0;
            double cy = 0;
            int cn = 0;

            for (auto& point : workField.points)
            {
                if (point.nearestCentroidNum == centroidNum)
                {
                    cx += point.x;
                    cy += point.y;
                    cn++;
                }
            }
            Point newCentroid;
            newCentroid.Init(cx/cn, cy/cn);
            return newCentroid;
        }
        //-----------------------------------------------------------------------------------

        void clastersToWorkField(Clusters *cls)
        {
            int i = 0;
            int clsNum = 1;
            for (auto& sp : cls->clusters)
            {
                for (auto& point : sp.points)
                {
                    point.nearestCentroidNum = clsNum;
                    workField.points[i] = point;
                    i++;
                }
                clsNum++;
            }
        }
        
    
    public:
        void Init(SetPoints pfield)
        {
            workField = pfield;
        }

        void kmeans(int k)
        {
            defineKCentroids(k);
            //REPLACE Iterations by Delta
            for (int i = 0; i < 100; ++i)
            {
                for (auto& point : workField.points)
                {
                    point.nearestCentroidNum = findNearestCentroid(k, point.x, point.y);
                }
                
                for (int j = 1; j <= k; ++j)
                {
                    centroids[j] = resetCenroid(k, j);
                }
                
            }
            printResult(&workField);
        }

        void wavecluster(int interval)
        {
            workGridField = redefineField(interval);
            
            connectComponents(&workGridField, interval);   

            gridToWorkField(interval); 

            printResult(&workField);
        }

        void hierarchy()
        {
            Clusters allClusters;
            allClusters.pointsToClusters(&workField);
            array<double, 3> ijr;

            double curDelt = 1;
            double prevDelt = 1;
            while (allClusters.clusters.size() > 1)
            {
                ijr = allClusters.getNearestClustersIdxs();
                curDelt = ijr[2];
                if (curDelt / prevDelt < 2)
                {
                    allClusters.clustersUnion((int)ijr[0], (int)ijr[1]);
                    prevDelt = curDelt;
                }
                else
                {
                    break;
                }


            }

            //!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            int i = 1;
            for (auto& cl : allClusters.clusters)
            {
                cout << i << ") " << cl.points.size() << endl;
                i++;
            }
            //!!!!!!!!!!!!!!!!!!!!!!!!
            
            clastersToWorkField(&allClusters);
            printResult(&workField);
            
        }
};


int main()
{   
    srand(time(NULL));

    SetPoints field;
    Algorithms alg;

    field.pointsInit(FILENAME);
    alg.Init(field);


    alg.kmeans(K);
    //alg.wavecluster(INTERVAL);
    //alg.hierarchy();

    return 0;
}
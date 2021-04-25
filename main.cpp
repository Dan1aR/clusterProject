//gnuplot visulization commands (Prefer to use analytics.ipynb)
//plot "< awk '{if($3 == \"0\") print}' input.txt" u 1:2 t "red" w p pt 2, "< awk '{if($3 == \"1\") print}' input.txt" u 1:2 t "green" w p pt 2, "< awk '{if($3 == \"2\") print}' input.txt" u 1:2 t "blue" w p pt 2
//plot "< awk '{if($3 == \"1\") print}' output.txt" u 1:2 t "green" w p pt 2, "< awk '{if($3 == \"2\") print}' output.txt" u 1:2 t "blue" w p pt 2, "< awk '{if($3 == \"3\") print}' output.txt" u 1:2 t "purpl" w p pt 2

#include <stdio.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <math.h>
#include <array>
#include <algorithm>
#include <iomanip>
#include <random>
#include <map>

using namespace std;

// GLOBAL VARIABLES
int K = 3;
int INTERVAL = 50; //Change to Treshhold!!!
double CDdivPD = 2.3;
float MAX = 1000;
string FILENAME_IN = "input.txt";
string FILENAME_OUT = "output.txt";


// Classes Prototypes;
class Point;
class Oval;
class SetPoints;
class KMeans;
class WaveCluster;
class HierarchyNew;
class Search;
class Interface;


// Usefull Functions
double sqr(double x){ return x * x; }
double degrees_to_radians(double degrees){ return degrees * M_PI / 180; }


// Classes
class Point
{
    public:
        double x;
        double y;
        unsigned short cluster_mark;

        Point(double px, double py, unsigned short pcm = 0)
        {
            x = px;
            y = py;
            cluster_mark = pcm;
        }

        void change_coordinates(double px, double py)
        {   
            x = px;
            y = py;
        }
};

class Oval
{
    private:
        // (x-x0)^2/a^2 + (y-y0)^2/b^2 <= r^2
        float x0;
        float y0;
        float a;
        float b;
        // alpha - смещение оси
        float alf;
        // number of points - количество точек в овале
        int number_of_points;

        bool is_inside(float x, float y)
        {
            if ( sqr((x-x0)*cos(alf) - (y-y0)*sin(alf))/sqr(a) + sqr((y-y0)*cos(alf) + (x-x0)*sin(alf))/sqr(b) <= 1 )
            {
                return true;
            }
            return false;
        }

    public:
        vector<Point> points;

        Oval(float px0, float py0, float pa, float pb, float palf, int pnop)
        {
            x0 = px0;
            y0 = py0;
            a = pa;
            b = pb;
            alf = palf;
            number_of_points = pnop;
        }

        void fill_oval()
        {
            int count = 0;
            double x, y;

            float RIGHT = max(x0+max(abs(a), abs(b)) , y0+max(abs(a), abs(b)));
            float LEFT = min(x0+min(-abs(a), -abs(b)) , y0+min(-abs(a), -abs(b)));

            while (count <= number_of_points)
            {
                x = (RIGHT-LEFT)*((double) rand() / RAND_MAX)+LEFT;
                y = (RIGHT-LEFT)*((double) rand() / RAND_MAX)+LEFT;

                if (is_inside(x, y))
                {
                    Point point(x, y, 0);
                    points.push_back(point);
                    count++;
                }
                
            }
                
        }

        void fill_oval(string T)
        {
            int count = 0;
            double x, y;

            float RIGHT = max(x0+max(abs(a), abs(b)) , y0+max(abs(a), abs(b)));
            float LEFT = min(x0+min(-abs(a), -abs(b)) , y0+min(-abs(a), -abs(b)));

            random_device rd{};
            mt19937 gen{rd()};
            // values near the mean are the most likely
            // standard deviation affects the   dispersion of generated values from the mean
            normal_distribution<> d{0.5, 0.25};

            while (count < number_of_points)
            {
                x = (RIGHT-LEFT)*(d(gen))+LEFT;
                y = (RIGHT-LEFT)*(d(gen))+LEFT;

                if (is_inside(x, y))
                {
                    Point point(x, y, 0);
                    points.push_back(point);
                    count++;
                }
                
            }
        }

        void move(float d_x0, float d_y0, float d_alf, float d_radius_vector_alf)
        {
            x0 += d_x0;
            y0 += d_y0;

            x0 = cos(d_radius_vector_alf)*x0 + sin(d_radius_vector_alf)*y0;
            y0 = -sin(d_radius_vector_alf)*x0 + cos(d_radius_vector_alf)*y0;

            alf += d_alf;

            for (auto& point : points)
            {
                point.x += d_x0;
                point.y += d_y0;

                point.x = cos(d_radius_vector_alf)*point.x + sin(d_radius_vector_alf)*point.y;
                point.y = -sin(d_radius_vector_alf)*point.x + cos(d_radius_vector_alf)*point.y;
                
                float m_x = point.x - x0;
                float m_y = point.y - y0;
                m_x = cos(d_alf)*m_x + sin(d_alf)*m_y;
                m_y = -sin(d_alf)*m_x + cos(d_alf)*m_y;
                point.x = m_x + x0;
                point.y = m_y + y0;
            }
        }
};

class SetPoints
{
    public:
        vector<Point> points;

        Point most_left = {MAX, 0};
        Point most_right = {-MAX, 0};
        Point most_bottom = {0, MAX};
        Point most_high = {0, -MAX};

        void define_spec_point(double x, double y)
        {
            if (x < most_left.x)
            {
                most_left.change_coordinates(x, y);
            }
            if (x > most_right.x)
            {
                most_right.change_coordinates(x, y);
            }
            if (y < most_bottom.y)
            {
                most_bottom.change_coordinates(x, y);
            }
            if (y > most_high.y)
            {
                most_high.change_coordinates(x, y);
            }

            if (most_high.y > MAX)
            {
                MAX = most_high.y + 1;
            }
            if (most_right.x > MAX)
            {
                MAX = most_right.x + 1;
            }
        }

        void init_points_from_file(string fileName)
        {
            double x;
            double y;
            int clasterNumber;

            ifstream pointsfile;
            pointsfile.open(fileName);

            while (pointsfile >> x >> y >> clasterNumber)
            {
                Point point(x, y, 0);
                points.push_back(point);
                define_spec_point(x, y);
            }

            pointsfile.close();
        }

        void add_oval(Oval *oval)
        {
            for (auto& oval_point : oval->points)
            {
                points.push_back(oval_point);
                define_spec_point(oval_point.x, oval_point.y);
            }
        }
};

//Old Ones Start
class Clusters
{
    private:
        double Rmin(SetPoints& U, SetPoints& V)
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

        void points_to_clusters(SetPoints *field)
        {
            for (auto& point : field->points)
            {
                SetPoints newCluster;
                newCluster.points.push_back(point);
                clusters.push_back(newCluster);
            }
        }

        void create_matrix_distance(vector< vector<float> > *mtx)
        {

            for (int i = 0; i < clusters.size(); ++i)
            {   
                vector<float> distances = {};
                for (int j = 0; j < clusters.size(); ++j)
                {
                    double R = Rmin(clusters[i], clusters[j]);
                    if (i != j)
                    {
                        distances.push_back(R);
                    }
                    if (i == j)
                    {
                        distances.push_back(MAX);
                    }
                }
                mtx->push_back(distances);
            }
        }

        array<double, 3> get_nearest_clusters_idxs()
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

        array<double, 3> get_nearest_clusters_idxs(vector< vector<float> > *mtx)
        {
            vector<float> mins;
            vector<int> minsJ;
            
            for (int i = 0; i < mtx->size(); ++i)
            {
                int minJ = distance(  (*mtx)[i].begin(), min_element((*mtx)[i].begin(), (*mtx)[i].end())  );
                //cout << minJ << endl;
                mins.push_back( (*mtx)[i][minJ] );
                minsJ.push_back(minJ);
            }

            int minIJ = distance(mins.begin(), min_element(mins.begin(), mins.end()) );
            //cout << "!!!!!!!!!!!!!!!!!! :: " << mins[0] << endl;
            array<double, 3> ijr = { (double)minIJ, (double)minsJ[minIJ], mins[minIJ]};
            return ijr;
        }

        void clusters_union(int i, int j)
        {
            for (auto& point : clusters[ j ].points)
            {
                clusters[ i ].points.push_back(point);
            }
            clusters.erase(clusters.begin() + j );
        }

        void clusters_union(int i, int j, vector< vector<float> > *mtx)
        {
            //cout << "FLAG " << i << " " << j << endl;
            for (int l = 0; l < clusters.size(); ++l)
            {
                (*mtx)[j][l] = MAX;
                (*mtx)[l][j] = MAX;
            }

            for (auto& point : clusters[ j ].points)
            {
                clusters[ i ].points.push_back(point);
            }
            clusters[j].points = {};

            
            for (int l = 0; l < clusters.size(); ++l)
            {   
                if ( (i != j) && ( (*mtx)[i][l] != MAX))
                {   
                    double R = Rmin(clusters[i], clusters[l]);
                    (*mtx)[i][l] = R;
                    (*mtx)[l][i] = R;
                }
                else
                {
                    (*mtx)[i][l] = MAX;
                }
            }
        }
};
class Hierarchy
{
    private:
        SetPoints *workField;

        void clastersToWorkField(Clusters *cls)
        {
            int i = 0;
            int clsNum = 1;
            for (auto& sp : cls->clusters)
            {
                if (sp.points.size() > 0)
                {
                    for (auto& point : sp.points)
                    {
                        point.cluster_mark = clsNum;
                        workField->points[i] = point;
                        i++;
                    }
                    clsNum++;
                }
            }
        }

    public:
        Hierarchy(SetPoints *p_work_field)
        {
            workField = p_work_field;
        }

        void hierarchy()
        {
            Clusters allClusters;
            allClusters.points_to_clusters(workField);
            vector< vector<float> > matrixDistnce;
            array<double, 3> ijr;

            allClusters.create_matrix_distance(&matrixDistnce);

            double curDelt = 1;
            double prevDelt = 1;

            while (true)
            {
                ijr = allClusters.get_nearest_clusters_idxs(&matrixDistnce);
                prevDelt = curDelt;
                curDelt = ijr[2];

                cout << curDelt << " " << prevDelt << " :: " << curDelt/prevDelt << endl;
                if (curDelt/prevDelt >= CDdivPD)
                {
                    break;
                }

                allClusters.clusters_union((int)ijr[0], (int)ijr[1], &matrixDistnce);
            }
            
            clastersToWorkField(&allClusters);
        }
};
//Old Ones End

// Cluster-Algorithms Goes Here
class KMeans
{
    private:
        int k;
        SetPoints *workField;
        vector<Point> centroids;
        void define_k_centroids(int k)
        {
            for (int i = 0; i < k; ++i)
            {
                double x = (workField->most_right.x - workField->most_left.x)*((double) rand() / RAND_MAX) - workField->most_right.x;
                double y = (workField->most_high.y - workField->most_bottom.y)*((double) rand() / RAND_MAX) - workField->most_high.y;
                Point centroid(x, y);
                centroids.push_back(centroid);
                //cout << centroids[i].x << " " << centroids[i].y << endl;
            }
        }

        int find_nearest_centroid(int k, double x, double y)
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

        Point reset_cenroid(int k, int centroid_num)
        {
            //Some bad code here
            double cx = 0;
            double cy = 0;
            int cn = 0;

            for (auto& point : workField->points)
            {
                if (point.cluster_mark == centroid_num)
                {
                    cx += point.x;
                    cy += point.y;
                    cn++;
                }
            }
            Point newCentroid(cx/cn, cy/cn);
            return newCentroid;
        }

    public:
        KMeans(int p_k, SetPoints *p_work_field)
        {
            k = p_k;
            workField = p_work_field;
        }

        void kmeans()
        {
            cout << "K is :: " << k << endl;
            define_k_centroids(k);
            //REPLACE Iterations by Delta
            for (int i = 0; i < 100; ++i)
            {
                for (auto& point : workField->points)
                {
                    point.cluster_mark = find_nearest_centroid(k, point.x, point.y);
                }
                
                for (int j = 1; j <= k; ++j)
                {
                    centroids[j] = reset_cenroid(k, j);
                }
                
            }
        }

};

class WaveCluster
{
    private:
        SetPoints workGridField;
        SetPoints *workField;
        int interval;

        SetPoints redefineField(int interval)
        {
            SetPoints gridField;

            double deltaX = (workField->most_right.x - workField->most_left.x) / interval;
            double deltaY = (workField->most_high.y - workField->most_bottom.y) / interval;

            double y0 = workField->most_bottom.y;
            double y1 = y0 + deltaY;
            
            vector<Point> newPoints;

            for (int i = 0; i < interval; ++ i)
            {
                double x0 = workField->most_left.x;
                double x1 = x0 + deltaX;

                for (int j = 0; j < interval; ++j)
                {
                    double cx = 0;
                    double cy = 0;
                    double cn = 0;

                    for (auto& point : workField->points)
                    {
                        if ((point.x > x0) && (point.x < x1) && (point.y > y0) && (point.y < y1))
                        {
                            cx += point.x;
                            cy += point.y;
                            cn++;
                        }
                    }

                    Point newPoint(0, 0, 0);
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
                gridField.define_spec_point(point.x , point.y);
            }

            return gridField;
        } 
        
        void recConnectCluster(SetPoints *field, int interval, int myi, int myj, int mycn)
        {
            if ((myi >= 0) && (myj >= 0) && (myi < interval) && (myj < interval))
            {
                if ( (field->points[interval*(myi)+ myj].x != 0) && (field->points[interval*(myi)+ myj].cluster_mark == 0) )
                {
                    //cout << mycn << ") " << myi << " " << myj << endl;
                    field->points[interval*(myi)+ myj].cluster_mark = mycn;

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
                        if (field->points[interval*(i)+ j].cluster_mark == 0)
                        {
                            
                            field->points[interval*(i)+ j].cluster_mark = cn;
                            
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
            double deltaX = (workField->most_right.x - workField->most_left.x) / interval;
            double deltaY = (workField->most_high.y - workField->most_bottom.y) / interval;

            double y0 = workField->most_bottom.y;
            double y1 = y0 + deltaY;

            for (int i = 0; i < interval; ++ i)
            {
                double x0 = workField->most_left.x;
                double x1 = x0 + deltaX;

                for (int j = 0; j < interval; ++j)
                {
                    for (auto& point : workField->points)
                    {
                        if ((point.x > x0) && (point.x < x1) && (point.y > y0) && (point.y < y1))
                        {
                            point.cluster_mark = workGridField.points[interval*i + j].cluster_mark;
                        }
                    }

                    x0 = x1;
                    x1 += deltaX;
                }
                
                y0 = y1;
                y1 += deltaY;
            }
        }

    public:
        WaveCluster(int p_interval, SetPoints *p_work_field)
        {
            interval = p_interval;
            workField = p_work_field;
        }

        void wavecluster()
        {  
            workGridField = redefineField(interval);
            
            connectComponents(&workGridField, interval);   

            gridToWorkField(interval);
        }
};

class HierarchyNew
{
    private:
        SetPoints *work_field;
        unsigned short N;

        vector <double> distances;
        vector < array<int, 2> > pairs;
        map <int, int> num_points_in_cluster;
        map <int, int> change_colors;

        double R(Point &p1, Point &p2){ return sqrt(sqr(p1.x - p2.x) + sqr(p1.y - p2.y)); }
        
        int num_non_zeros_clusters(map<int, int>& M)
        {
            int n = 0;
            for (auto it : M)
            {
                if (it.second > 0)
                {
                    ++n;
                }
            }
            return n;
        }

        double partition (vector<double> *arr, int low, int high, vector< array<int, 2> > *sub_arr)
        {
            double pivot = arr->at(high);
            int i = (low - 1);
        
            for (int j = low; j <= high- 1; j++)
            {
                if (arr->at(j) <= pivot)
                {
                    i++;
                    swap(arr->at(i), arr->at(j));
                    swap(sub_arr->at(i), sub_arr->at(j));
                }
            }
            swap(arr->at(i + 1), arr->at(high));
            swap(sub_arr->at(i + 1), sub_arr->at(high));
            return (i + 1);
        }

        void quickSort(vector<double> *arr, int low, int high, vector< array<int, 2> > *sub_arr)
        {
            if (low < high)
            {
                int pi = partition(arr, low, high, sub_arr);
                quickSort(arr, low, pi - 1, sub_arr);
                quickSort(arr, pi + 1, high, sub_arr);
            }
        }

        void create_distances_vector()
        {   
            for (int i = 0; i < work_field->points.size(); ++i)
            {
                for (int j = i+1; j < work_field->points.size(); ++j)
                {
                    distances.push_back(R(work_field->points[i], work_field->points[j]));
                    pairs.push_back({i, j});
                }
            }

        }

        void find_real_cluster(unsigned short *_cm)
        {
            while (num_points_in_cluster[ *_cm ] == 0) 
            { 
                *_cm = change_colors[*_cm]; 
            }
        }

        void color_pair(int _i, int _j, int *_cm)
        {
            if ( (work_field->points[_i].cluster_mark == 0) && (work_field->points[_j].cluster_mark == 0) )
            {
                work_field->points[_i].cluster_mark = *_cm;
                work_field->points[_j].cluster_mark = *_cm;
                num_points_in_cluster[*_cm] = 2; 
                *_cm += 1;
            }          
            else if ( work_field->points[_i].cluster_mark * work_field->points[_j].cluster_mark == 0 )
            {
                unsigned short new_cm = max(work_field->points[_i].cluster_mark, work_field->points[_j].cluster_mark);

                find_real_cluster(&new_cm);

                work_field->points[_j].cluster_mark = new_cm;
                work_field->points[_i].cluster_mark = new_cm;

                num_points_in_cluster[new_cm] += 2;
            }  
            else if ( (work_field->points[_i].cluster_mark != 0) && (work_field->points[_j].cluster_mark != 0) )
            {
                find_real_cluster(&work_field->points[_i].cluster_mark);
                find_real_cluster(&work_field->points[_j].cluster_mark);

                if (num_points_in_cluster[work_field->points[_i].cluster_mark] < num_points_in_cluster[work_field->points[_j].cluster_mark])
                {
                    change_colors[ work_field->points[_i].cluster_mark ] = work_field->points[_j].cluster_mark;
                    num_points_in_cluster[work_field->points[_i].cluster_mark] = 0;
                    num_points_in_cluster[work_field->points[_j].cluster_mark] += 1;
                }
                else
                {
                    change_colors[ work_field->points[_j].cluster_mark ] = work_field->points[_i].cluster_mark;
                    num_points_in_cluster[work_field->points[_j].cluster_mark] = 0;
                    num_points_in_cluster[work_field->points[_i].cluster_mark] += 1;
                }

            }
        }

        void color_last_one()
        {
            for (auto& point : work_field->points)
            {
                find_real_cluster(&point.cluster_mark);
            }
        }

    public:
        HierarchyNew(SetPoints *p_work_field, unsigned short p_N)
        {
            work_field = p_work_field;
            N = p_N;
        }

        void hierarchy()
        {
            create_distances_vector();
            //cout << work_field->points.size() << " : " << distances.size() << endl;
            quickSort(&distances, 0, distances.size()-1, &pairs);

            cout << "Sorted..." << endl;
            
            int cluster_mark = 1;

            for (int i = 0; i < distances.size(); ++i)
            {
                //cout << i << " ";
                color_pair(pairs[i][0], pairs[i][1], &cluster_mark);

                if ((i > N*N*N) && (num_non_zeros_clusters(num_points_in_cluster) <= N))
                {
                    break;
                }
            }
            color_last_one();
        }
};

class Search
{
    private:
        void preprocess_field_to_save(SetPoints *pwork_field)
        {
            vector<int> cluster_marks;
            for (auto& point : pwork_field->points)
            {
                if ( find(cluster_marks.begin(), cluster_marks.end(), point.cluster_mark) == cluster_marks.end() )
                {
                    cluster_marks.push_back( point.cluster_mark );
                }
            }
            for (auto& point : pwork_field->points)
            {
                point.cluster_mark = find(cluster_marks.begin(), cluster_marks.end(), point.cluster_mark) - cluster_marks.begin() + 1;
            }
        }

        void save_field_to_file(SetPoints *pwork_field, string file_name)
        {
            ofstream pointsfile;
            pointsfile.open(file_name);
            //int i = 0;
            preprocess_field_to_save(pwork_field);
            for (auto& point : pwork_field->points)
            {
                if ( (point.x != 0) && (point.y != 0) )
                {
                    //cout << point.x << " " << point.y << " " << point.nearestCentroidNum << "\n";
                    pointsfile << point.x << " " << point.y << " " << point.cluster_mark << "\n";
                    //i++;
                }
            }
            //cout << i << endl;
        }

        

    public:
        SetPoints *workField;

        Search(SetPoints *pfield)
        {
            workField = pfield;
        }

        void kmeans(int k)
        {
            KMeans kmeans_alg(k, workField);
            kmeans_alg.kmeans();
            save_field_to_file(workField, FILENAME_OUT);
        }

        void wavecluster(int interval)
        {
            WaveCluster wavecluster_alg(interval, workField);
            wavecluster_alg.wavecluster();
            save_field_to_file(workField, FILENAME_OUT);
        }

        void hierarchy()
        {
            Hierarchy hierarchy_alg(workField);
            hierarchy_alg.hierarchy();
            save_field_to_file(workField, FILENAME_OUT);
        }

        void hierarchyNew(int k)
        {
            HierarchyNew hierarchyNew_alg(workField, k);
            hierarchyNew_alg.hierarchy();
            save_field_to_file(workField, FILENAME_OUT);
        }
};

class Interface
{
    private:
        SetPoints *work_field;
        Search *alg;
    
        void place_ovals()
        {
            cout << "1 :: Get clouds positions from file" << endl;
            cout << "2 :: Get clouds positions one by one" << endl;
            int command = 1;
             
            cout << ">>> "; cin >> command;
            if (command == 1)
            {
                string ovals_file_name;
                cout << "File name with clouds:: ";
                 
                cout << ">>> "; cin >> ovals_file_name;

                ifstream pointsfile;
                pointsfile.open(ovals_file_name);

                float p_x0;
                float p_y0;
                float p_a;
                float p_b;
                float p_alf;
                int p_nop;

                while (pointsfile >> p_x0 >> p_y0 >> p_a >> p_b >> p_alf >> p_nop)
                {
                    p_alf = degrees_to_radians(p_alf);
                    Oval oval(p_x0, p_y0, p_a, p_b, p_alf, p_nop);
                    oval.fill_oval("N");
                    work_field->add_oval(&oval);
                }
                pointsfile.close();

            }
            if (command == 2)
            {
                while (command != 0)
                {
                    cout << "1 :: Add cloud" << endl;
                    cout << "0 :: end" << endl;
                     
                    cout << ">>> "; cin >> command;

                    if (command == 1)
                    {
                        cout << "Enter params in format :: x0, y0, a, b, alf, number_of_points" << endl;
                        float p_x0, p_y0, p_a, p_b, p_alf, p_nop;
                        cout << ">>> "; cin >> p_x0 >> p_y0 >> p_a >> p_b >> p_alf >> p_nop;
                        p_alf = degrees_to_radians(p_alf);
                        Oval oval(p_x0, p_y0, p_a, p_b, p_alf, p_nop);
                        oval.fill_oval("N");

                        int cloud_move_command = 1;
                        while (cloud_move_command)
                        {
                            cout << "Want to move cloud?" << endl;
                            cout << "1 :: YES" << endl;
                            cout << "0 :: NO" << endl;
                            cout << ">>> "; cin >> cloud_move_command;
                            
                            if (cloud_move_command)
                            {
                                cout << "Enter params in format :: d_x0, d_y0, d_alf, d_radius_vector_alf" << endl;
                                float d_x0, d_y0, d_alf, d_radius_vector_alf;
                                cout << ">>> "; cin >> d_x0 >> d_y0 >> d_alf >> d_radius_vector_alf;
                                d_alf = degrees_to_radians(d_alf);
                                d_radius_vector_alf = degrees_to_radians(d_radius_vector_alf);
                                oval.move(d_x0, d_y0, d_alf, d_radius_vector_alf);
                            }
                        }

                        work_field->add_oval(&oval);
                    } 
                }
            }
        }

    public:
        void interact()
        {   
            cout << "Heyyyyyy, you are in a clustering programm!" << endl;

            int command = 1;
            while (command != 0)
            {
                cout << "1 :: Define your points by ellipse-cloud" << endl;
                cout << "2 :: start K-mean (Needs number of Clusters)" << endl;
                cout << "3 :: start WaveCluster (Needs number of Grid-Intervals)" << endl;
                cout << "4 :: start Hierarchy (Doesn't need clusters number)" << endl;
                cout << "5 :: start HierarchyNew (Needs number of Clusters)" << endl; 
                cout << "0 :: Exit" << endl;

                cout << ">>> "; cin >> command;
                if (command == 1)
                {
                    place_ovals();
                    /*alg->save_field_to_file(work_field, FILENAME_OUT);*/
                    cout << "All clouds are defined!" << endl;
                }
                if (command == 2)
                {
                    cout << "Define K::" << endl;
                    cout << ">>> "; cin >> K;
                    alg->kmeans(K);
                    cout << "K-means results saved to :: " << FILENAME_OUT << endl;
                }
                if (command == 3)
                {
                    cout << "Define Grid Interval::" << endl;
                    cout << ">>> "; cin >> INTERVAL;
                    alg->wavecluster(INTERVAL);
                    cout << "WaveCluster results saved to :: " << FILENAME_OUT << endl;
                }
                if (command == 4)
                {
                    alg->hierarchy();
                    cout << "Hierarchy results saved to :: " << FILENAME_OUT << endl;
                }
                if (command == 5)
                {
                    cout << "Define K::" << endl;
                    cout << ">>> "; cin >> K;
                    alg->hierarchyNew(K);
                    cout << "HierarchyNew results saved to :: " << FILENAME_OUT << endl;
                }
            } 

        }

        Interface(SetPoints *p_field, Search *p_alg)
        {
            work_field = p_field;
            alg = p_alg;
        }

};


int main(int argc, char *argv[])
{   
    srand(time(NULL));

    SetPoints field;
    Search cluster_search_algorithms(&field);
    Interface user_interface(&field, &cluster_search_algorithms);

    user_interface.interact();

    return 0;
}
/kmeans clustering program


#include <iostream>
#include <cmath>
#include <time.h>
#include <vector>


using namespace std;
void  calcCentriod(int x, int y, vector<int> Coor, double size){
  int xMean, yMean, xory;
  xory = 0;
  cout << "calc centroid" << endl;
  for (auto i = Coor.begin(); i != Coor.end(); ++i){
    if(xory%2 == 0)
      xMean += Coor[*i];
    if(xory%2 != 0)
      yMean += Coor[*i];
    xory++;
  }
  x = xMean / (size/2);
  y = yMean / (size/2);
}

double distanceCalc(int x1, int x2, int y1, int y2){
  return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2) * 1.0);
}


int main(){

  srand(time(NULL));

  int NumberOfIterations = 10;
  int distanceThresh = rand() % 15 + 1;
  int kClusters = rand() % 10 + 1;
  cout << "You have "  << kClusters -2 << " clusters." << endl;
  int kPoints = rand()% 50 + 1;
  cout << "You have "  << kPoints << " points." << endl;

  int xPoints[kPoints], yPoints[kPoints];

  int xCentriods[kClusters], yCentriods[kClusters];

  cout << "pushing intial clusters" << endl;
  vector< vector<int> > clusters(kClusters+1);
  for(int i = 0; i < kClusters + 1; ++i){
    clusters[i].push_back(i);
  }

  cout << "pushing centriods" << endl;
  for(int i = 0; i < kClusters - 1; ++i){
    xCentriods[i] = rand() % 101;
    yCentriods[i] = rand() % 101;
  }

  cout << "pushing points" << endl;
  for(int i = 0; i < kPoints; ++i){
    xPoints[i] = rand() % 101;
     std::cout << xPoints[i] << std::endl;
    yPoints[i] = rand() % 101;
    std::cout << yPoints[i] << std::endl;
  }

  for(int iter = 1; iter < NumberOfIterations; ++iter){

  //  double lowestDistance = 100.0;
  int closestCluster = 0;
  for(int i  = 0; i < kPoints; ++i){
    double lowestDistance = 100.0;
    for(int j = 0; j < kClusters; ++j){
      double currDistance = distanceCalc(xPoints[i], xCentriods[j], yPoints[i], yCentriods[j]);
      lowestDistance = min(currDistance, lowestDistance);
      if(currDistance == lowestDistance)
        closestCluster = j;
    }
    clusters[closestCluster].push_back(xPoints[i]);
      //std::cout << clusters[closestCluster] << std::endl;
    clusters[closestCluster].push_back(yPoints[i]);
    //std::cout << *clusters[closestCluster] << std::endl;

  //vector<int> temp;

  for(int u = 0; u < kClusters - 2; ++u){
    vector<int>::iterator IT;
    vector<int> temp(clusters[u].size() - 1, 0);
    for( IT = clusters[u].begin()+1 ; IT != clusters[u].end(); IT++){
      temp.push_back(*IT);
      std::cout << *IT << ", " << std::endl;
 for(int u = 0; u < kClusters - 2; ++u){
    vector<int>::iterator IT;
    vector<int> temp(clusters[u].size() - 1, 0);
    for( IT = clusters[u].begin()+1 ; IT != clusters[u].end(); IT++){
      temp.push_back(*IT);
      std::cout << *IT << ", " << std::endl;
    }
    calcCentriod(xCentriods[u], yCentriods[u], temp, ((clusters[u].size()-1)));
    std::cout << "New Centriod " << u << " : " << xCentriods[u] << ", " << yCentriods[u] << endl;
  }
}
  }
}


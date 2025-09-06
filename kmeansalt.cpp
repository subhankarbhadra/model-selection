#include <RcppArmadillo.h>

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
// Function to assign each point to the nearest cluster
List assignPointsToClusters(arma::mat points, arma::cube Projections) {
  arma::uword nPoints = points.n_rows;
  arma::uword nClusters = Projections.n_slices;
  arma::ivec clusterAssignments(nPoints);
  arma::vec dist(nPoints);
  
  for (arma::uword i = 0; i < nPoints; i++) {
    double minDistance = INFINITY;
    int minIndex = -1;
    for (arma::uword j = 0; j < nClusters; j++) {
      double distance = arma::norm(points.row(i).t() - Projections.slice(j) * points.row(i).t(), 2);
      if (distance < minDistance) {
        minDistance = distance;
        minIndex = j;
      }
    }
    dist(i) = pow(minDistance, 2);
    clusterAssignments(i) = minIndex;
  }
  return List::create(Named("dist") = dist, Named("clusters") = clusterAssignments);
}

// Function to update the projection matrix for each cluster
// [[Rcpp::export]]
arma::cube updateProjections(arma::mat points, arma::ivec clusterAssignments, int nClusters, int d) {
  int nDimensions = points.n_cols;
  arma::cube Projections(nDimensions, nDimensions, nClusters);
  NumericVector nPointsInCluster(nClusters);
  for (arma::uword i = 0; i < points.n_rows; i++) {
    int clusterIndex = clusterAssignments(i);
    nPointsInCluster[clusterIndex]++;
  }
  for (int i = 0; i < nClusters; i++) {
    if (nPointsInCluster[i] > 0) {
      arma::mat X = points.rows(arma::find(clusterAssignments == i));
      arma::mat U, V;
      arma::vec s;
      arma::svd_econ(U, s, V, X.t(), "left");
      arma::mat Ud = U.cols(0, d-1);
      arma::mat proj = Ud * trans(Ud);
      Projections.slice(i) = proj;
    }
  }
  return Projections;
}

// [[Rcpp::export]]
List kmeansAlt(arma::mat points, int nClusters, int d, int maxIterations, int nStart, double tol=0.002) {
  int nPoints = points.n_rows;
  int nDimensions = points.n_cols;
  arma::ivec clusterAssignmentsFinal(nPoints);
  double mintotaldist = INFINITY;
  int iterationFinal = -1;
  
  for(int l = 0; l < nStart; l++) {
    // Initialize cluster assignments randomly
    arma::ivec clusterAssignments = arma::randi(nPoints, arma::distr_param(0, nClusters-1));
    arma::cube Projections(nDimensions, nDimensions, nClusters);
    arma::vec dist(nPoints);
    
    double totaldist = INFINITY;
    double oldtotaldist = INFINITY;
    double gap = INFINITY;
    int iteration = 0;
    
    while(gap > tol && iteration < maxIterations) {
      oldtotaldist = totaldist;
      // Update the centroids for each cluster
      Projections = updateProjections(points, clusterAssignments, nClusters, d);
      
      // Assign each point to the nearest cluster
      List out = assignPointsToClusters(points, Projections);
      clusterAssignments = as<arma::ivec>(out["clusters"]);
      dist = as<arma::vec>(out["dist"]);
      totaldist = sum(dist);
      gap = std::abs(oldtotaldist - totaldist);
      iteration++;
    }
    if(mintotaldist > totaldist) {
      mintotaldist = totaldist;
      clusterAssignmentsFinal = clusterAssignments;
      iterationFinal = iteration;
    }
  }
  
  return List::create(Named("cluster") = clusterAssignmentsFinal + 1, Named("error") = mintotaldist, Named("iter") = iterationFinal);
}

// [[Rcpp::export]]
List kmeansAltV2(arma::mat points, int nClusters, int d, int maxIterations, int nStart, double tol=0.01, int batch_size=50) {
  int nPoints = points.n_rows;
  int nDimensions = points.n_cols;
  arma::ivec clusterAssignmentsFinal(nPoints);
  double mintotaldist = INFINITY;
  int iterationFinal = -1;
  
  for(int l = 0; l < nStart; l++) {
    // Initialize cluster assignments randomly
    arma::ivec clusterAssignments = arma::randi(nPoints, arma::distr_param(0, nClusters-1));
    arma::cube Projections(nDimensions, nDimensions, nClusters);
    arma::vec dist(nPoints);
    arma::uvec randvec(batch_size);
    
    double totaldist = INFINITY;
    //double oldtotaldist = INFINITY;
    //double gap = INFINITY;
    arma::ivec oldclusterAssignments(nPoints);
    arma::uword numMismatches = 1;
    double prop = 1;
    int iteration = 0;
    
    while(prop > tol && iteration < maxIterations) {
      oldclusterAssignments = clusterAssignments;
      //oldtotaldist = totaldist;
      // Update the centroids for each cluster
      randvec = arma::randperm(nPoints, batch_size);
      Projections = updateProjections(points.rows(randvec), clusterAssignments(randvec), nClusters, d);
      
      // Assign each point to the nearest cluster
      List out = assignPointsToClusters(points, Projections);
      clusterAssignments = as<arma::ivec>(out["clusters"]);
      
      numMismatches = arma::accu(clusterAssignments != oldclusterAssignments);
      
      // Calculate the proportion of mismatches
      prop = static_cast<double>(numMismatches) / clusterAssignments.size();
      
      dist = as<arma::vec>(out["dist"]);
      totaldist = sum(dist);
      // gap = abs(oldtotaldist - totaldist);
      iteration++;
    }
    if(mintotaldist > totaldist) {
      mintotaldist = totaldist;
      clusterAssignmentsFinal = clusterAssignments;
      iterationFinal = iteration;
    }
  }
  
  return List::create(Named("cluster") = clusterAssignmentsFinal + 1, Named("error") = mintotaldist, Named("iter") = iterationFinal);
}


#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cctype>
#include <cstring>
#include <Rcpp.h>
#include <vector>
#include <string>
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::Range]]
#include <progress.hpp>
#include <progress_bar.hpp>

using namespace std;
// All variables that do not depend on R data
#define sep ' '
int Z;                         //number of members of most recent spheres
double reward;                 //reward value
//number of samples


char buf[5000];                //input buffer
char nm[1000];                 //name complilation buffer
double minn,MAXX;              //min,max data distance

// Stuff for error reporting




// All functions required


double MyD(int i,int j, int& features, Rcpp::NumericMatrix& data){//Euclidian distance between lines of data

  double delta,accu;  //coordinate difference and acumulator
  int k;              //loop index

  accu=0.0;  //zero the accumulator

  for(k=0;k<features;k++){//loop over coordinates
    delta=data(i, k)-data(j, k);  //compute difference
    accu+=(delta*delta);          //sum to accumulator
  }

  return(accu);  //return distance

}

void findminMAX(int&lines, int& features, Rcpp::NumericMatrix& distMat){//find the distance parameters

  int i,j;     //loop index variables
  double cp;   //comparison variables

  minn=MAXX=distMat(0,1);  //initialize distances
  for(i=0;i<lines-1;i++)for(j=i+1;j<lines;j++){
    cp=distMat(i, j); //get next number for comparison
    if(minn>cp)minn=cp;  //new minimum?
    if(MAXX<cp)MAXX=cp;  //new maximum?
  }

}

void addTo(int addVect[], int*** arrayPnt, int locationX, int locationY, int nSamples){
  for(int i = 0; i < nSamples; i++){
    arrayPnt[locationX][locationY][i] = int (addVect[i]);
  }
}

void memberate(int v,double r, Rcpp::NumericMatrix& orderMat, Rcpp::NumericMatrix& distMat, int* member, int*** mArray ){//find membership

  int size = distMat.ncol();
  int ls = 0;
  int rs = size - 1;
  int check;

  while (rs - ls > 1){
    check = floor((rs - ls) / 2) + ls;
    int dist_index = orderMat(v, check);
    if (distMat(v, dist_index) < r){
      ls = check;}
      else
       { rs = check;}
    }
    double Z = ls+1;

  reward=1.0/(Z); //compute the reward
  int firstEntry = mArray[v][ls][0];
 if(firstEntry!=-1){
   for(int i=0;i<size;i++){
     member[i] = mArray[v][ls][i];
   }}
 else{
   mArray[v][ls][0]=0;
   for(int i=0;i<ls+1;i++){
     int j=orderMat(v,i);
     mArray[v][ls][j]=1;
   }
   for(int i=0;i<size;i++){
     member[i] = mArray[v][ls][i];
   }
 }
}

/*

 void memberate(int v,double r, int& lines, int& features, Rcpp::NumericMatrix& distMat, int* member){//find membership

  int i;  //loop index

  Z=0;  //number of members
  //find sphere members
  for(i=0;i<lines;i++){
    double dist=distMat(v,i);
   {if(dist<=r){member[i]=1;Z++;} else member[i]=0;}}

  reward=1.0/((double)Z); //compute the reward
}
 */



void Flatten(int&lines, Rcpp::NumericMatrix& A){//zero the assoicator

  int i,j;  //loop index

  //zero out the assiciator
  for(i=0;i<lines;i++)for(j=0;j<lines;j++)A(i, j)=0.0;

}

void Associate(int& lines, int* member, Rcpp::NumericMatrix& A){//update the associator with the reward

  int i,j;  //loop index

  for(i=0;i<lines;i++)if(member[i]==1){//inside the sphere
    for(j=0;j<lines;j++)if((member[j]==1)&&(i!=j)){//active pair
      A(i, j)+=reward;  //add in associator
    }
  }
}

double RadiiDistU(double MAXX, double sizeB,double minn){ //Generates uniform random radius between minn and MAXX
  double randNum=unif_rand();
  return(randNum*MAXX*sizeB+minn);
}

double RadiiDistQ(double MedianRadius){ //quadratically favors shorter radii
  double u=0;
  do{u=unif_rand();}while(u<0.0001); //get a random number that will not cause instability
  return(MedianRadius*u/(1-u));
}

//Warning -- this function uses a Median Radius about 75% of the real median radius to work right
double RadiiDistE(double MedianRadius){ //exponenitally favors shorter radii

  double u=0;
  do{u=unif_rand();}while(u<0.0001); //get a random number that will not cause instability

  return(MedianRadius*(-log(1-u)));   //recall that the log() in math.h IS natural log
}


//Generate a distance matrix

Rcpp::NumericMatrix gen_dist(Rcpp::NumericMatrix data){
	Rcpp::NumericMatrix returnMat(data.nrow(), data.nrow());
	int features = data.ncol();
	 for (int i = 0; i < returnMat.nrow(); i++) {
      for (int j = 0; j < i; j++) {
		  double dist = MyD(i, j, features, data);
		  returnMat(i,j) = dist;
		  returnMat(j,i) = dist;
	  }
	 }
	return returnMat;
}


// Generate matrix with samples in ascending distance order

Rcpp::NumericMatrix orderMatrix(Rcpp::NumericMatrix data_matrix){
  Rcpp::NumericMatrix distanceMat(data_matrix);
  int nl = distanceMat.nrow();
  Rcpp::NumericMatrix returnMat(nl);

  for (int k=0;k<nl;k++){ // for each possible center
    vector<int> rowOrder(nl, -1 );
    //smallest distance is always itself
    rowOrder.insert(rowOrder.begin(),k);
    for (int j=0;j<nl;j++){ //for each possible point has a distance with
      if (j==k){continue;} // k already in vector
      for (int i=0; i<nl; i++){//go through the vector from start to finish
        if (rowOrder[i] == -1){rowOrder.insert(rowOrder.begin()+i, j); break;}
        // if the distance between k and j is larger than k and rowOrder[i]
        //keep looking
        int ithPos;
        ithPos = rowOrder[i];
        if (distanceMat(k,ithPos)<distanceMat(k,j)) { continue;}
        else{ // place j into the position where dist(j) switches from being
          //smaller to larger or equal to dist(i) and stop
          rowOrder.insert(rowOrder.begin()+i, j);
          break;
        }
      }
    // populate row of NumericMatrix
    for (int l=0;l<nl;l++){
      returnMat(k,l) = rowOrder[l];
    }
  }
}
  return returnMat;
}

// Main R program
RcppExport SEXP bubble(SEXP data_matrix, SEXP samples, SEXP columns, SEXP bubbleLevel, SEXP bubbleSize, SEXP RadiiDist, SEXP RewardType, SEXP logic){//main routine
  //extract the information that is passed to C++ and put it in usable forms

  Rcpp::NumericVector error1(1);

  Rcpp::NumericMatrix data(data_matrix); //data matrix
  Rcpp::NumericVector dim1(samples);
  int lines=dim1[0]; //extract #lines from R data
  Rcpp::NumericVector dim2(columns);
  int features=dim2[0]; //extract # of features from R data
  const int nSample = lines; // used for generating array sizes
  Rcpp::NumericVector Type(RadiiDist);
  int RadiiT=Type[0]; // extract the radii distribution type

  //generate all variables that require information from R
  Rcpp::NumericMatrix A(lines, lines); //make R data varaibale
  Rcpp::NumericMatrix distMat(lines, lines);
  // initialize Associator matrix with small + value
  fill(A.begin(), A.end(), 0.000001);

  int member[lines];             //characteristic function
  int i,j;        //loop indices
  double radius;  //critical radius
  Rcpp::NumericVector dim3(bubbleLevel);
  Rcpp::NumericVector dim4(bubbleSize);

  int numB=dim3[0]; //extract Bubble level from R data
  long double sizeB=dim4[0]; // extract Bubble max size

  Rcpp::NumericVector dim5(RewardType);
  Rcpp::NumericVector logicV(logic);

  int SAM=lines*pow(10, numB);

  //generate initial membership array
  int*** memArray = new int**[nSample];
  for(int i = 0; i < nSample; i++){
    memArray[i] = new int*[nSample];
    for(int j = 0; j < nSample; j++){
      memArray[i][j] = new int[nSample];
    }
  }
    // and initialize
    for(int i = 0; i < nSample; i++){
      for(int j = 0; j < nSample; j++){
        for(int k=0; k< nSample; k++){
          memArray[i][j][k] = -1;
        }
      }

  }

  //generate distance matrix
  if (logicV[0]==0){
    distMat = gen_dist(data);
  }
  if (logicV[0]==1){
    distMat=data;
  }
  bool verbose=TRUE;

  if (logicV[1]==0){
    verbose = FALSE;
  }
  if (logicV[1]==1){
    verbose = TRUE;
  }

  findminMAX(lines, features, distMat);   //find the minimum and maximum values
  MAXX-=minn;     //normalize the maximum
  //cout<<MAXX<<endl;
  Progress p(SAM, verbose);

  //generate the sorted sample matrix

  Rcpp::NumericMatrix  oMat = orderMatrix(distMat);

  for(i=0;i<SAM;i++){//loop over samples

   //generate a random critical radius
    if(RadiiT==1){radius=RadiiDistU(MAXX, sizeB, minn);}
    if(RadiiT==2){
      double MedianRadius=((MAXX-minn)/2)+minn;
      radius=RadiiDistQ(MedianRadius);}
    if(RadiiT==3){
      double MedianRadius=(((MAXX-minn)/2)+minn)*.75;
      radius=RadiiDistE(MedianRadius);}
    if(RadiiT==4){
      double multiplier=0.49625006+0.39328344*sin(3.89098056*unif_rand()-1.93639904);
      radius=multiplier*MAXX*sizeB+minn;}
    if(RadiiT==5){
      double x=unif_rand();
      double multiplier=0.82947314*x*x*x-1.191853*x*x+1.1810129*x+0.09299366;
      radius=multiplier*MAXX*sizeB+minn;
    }


    //if(i%100==0){cout<<radius<<endl;}
    j=i%lines;          //generate core vertex
    memberate(j,radius, oMat, distMat, member, memArray);        //find the membership function
    Associate(lines, member, A);                //update the associator
    //print progress bar and check for abort
    if (Progress::check_abort() ){ return error1;}
    p.increment(); // update progress

  }
  cout<<endl;

  //clean up memory
  for (int i = 0; i < lines; i++)
  {
    for (int j = 0; j < lines; j++)
      delete[] memArray[i][j];

    delete[] memArray[i];
  }

  delete[] memArray;

  return A;     //pass the associator matrix back to R
}

void nextdiff(int &ind, int *label, int length, int other){
  int j;
  j=ind+1;

  while ((j<length)&&(label[other]==label[j])) j++;
  ind=j;
}

RcppExport SEXP ATree(SEXP associator_matrix){
  Rcpp::NumericVector error1(1);
  Rcpp::NumericMatrix Assoc(associator_matrix);
  int nl=Assoc.nrow(); //extract #lines from R data
  Rcpp::NumericMatrix Merge(nl-1 , 2); // empty merge table to be returned to R
  Rcpp::NumericVector Height(nl-1); //empty height vector to be returned to R
  Rcpp::NumericVector Order(nl); //empty order vector to be returned to R

  int mem[nl];    //create subtree membership
  int mergeRow[nl]; // contains last level of merge table to contain item


  //Make diagonal smaller to avoid it linking items to themselves

  for (int i=0; i<nl; i++) {
    Assoc(i,i)=Assoc(i,i)/10.0;
  }


  // fill subtree membership
  for (int i=0; i<nl; i++) {
    mem[i]=i;
    mergeRow[i]=0;
  }

  int a=0 , b =0;    //nodes to be joined
  int rowA=-1, rowB =-1; // record the row of a and b in the merge table if present
  int numal=0; //number of active leaves

  for(numal=0;numal<nl-1;numal++){//do nl-1 joins
    //find the join to do
    b=a;
    //nextdiff(b, mem, nl, a);
    for(int i=0;i<nl-1;i++){
      int j=i;
      while (j<nl){
        nextdiff(j, mem, nl, i);
        if (j<nl){
          if(Assoc(i,j)>Assoc(a,b)){a=i;b=j;}

        }
      }
    }

    int compare=mem[b];

    for(int i=0;i<nl;i++){

      if(mem[i]==compare) {

        mem[i]=mem[a];          //change all the b's to a's in mem
      }
    }

    // Print mem results to try and de-bug

    //  for(int i=0; i<nl; i++) { cout << mem[i] << " ";    }

    // cout << endl;

    //determine if a and be are already in the merge table and if so what row they are in

    for(int i=0; i<nl; i++) {
      if (mergeRow[a]>0) rowA=mergeRow[a];
      else rowA=-1;
    }

    for(int i=0; i<nl; i++) {
      if (mergeRow[b]>0) rowB=mergeRow[b];
      else rowB=-1;
    }

    //add the appropriate values to the merge table, a or b already present include +row else -a and -b
    if (rowA<0) Merge(numal,0)=-(a+1);
    else Merge(numal,0)= (rowA);
    if (rowB<0) Merge(numal,1)=-(b+1);
    else Merge(numal,1)= (rowB);

    //rearange pairs so that interior node always on the right

    if (Merge(numal,0)>Merge(numal,1)){
      int placehold=Merge(numal,1);
      Merge(numal,1)=Merge(numal,0);
      Merge(numal,0)=placehold;
    }

    for(int i=0;i<nl;i++){  //update MergeRow
      if (mem[i]==mem[a]) mergeRow[i]=numal+1;
    }
    Height[numal]= 1/Assoc(a,b);// get the hight of the interior nodes from the Association Matrix
  }

  vector<int> C_Order(nl);
  for (int i=0; i<2; i++) C_Order[i]= Merge(nl-2,i);  //grab the last line of the merge table -1 for index -1 because n-1 # of joins
  for (int i=0; i<nl; i++){//go through the vector from start to finish
    if (C_Order[i]>0) {
      int nodePoint=C_Order[i]-1; //index that node point is pointing towards
      C_Order[i]= Merge(nodePoint,0);
      C_Order.insert(C_Order.begin()+i, Merge(nodePoint,1));
      if (C_Order[i]>0) i--; // if the node point contains a node point in the first slot will revisit that point
    }
  }
  for (int i=0; i<nl; i++){ Order(i)=C_Order[i]*-1;}



  //combine all pieces into dataframe

  // Rcpp::List pieces= Rcpp::List::create(Merge, Height, Order);
  Rcpp::List pieces= Rcpp::List::create(Merge, Height, Order);

  return (pieces);

}

int contains(Rcpp::NumericMatrix merge, int node, int pointA, int pointB){
  int results=0;
  if (node<0){ //if leaf ->
    int leafA=(pointA+1)*-1;
    int leafB=(pointB+1)*-1;
    if (node==leafA) results+=1;
    if (node==leafB) results+=2;
  }else {
    //left
    int nodeL=merge((node-1),0);
    results+=contains(merge, nodeL, pointA, pointB );
    //right
    int nodeR=merge((node-1),1);
    results+=contains(merge, nodeR, pointA, pointB );
  }
  return (results);
}

RcppExport SEXP snip(SEXP margeTable, SEXP removed){
  Rcpp::NumericMatrix Merge(margeTable);
  int nl=Merge.nrow()+1; //number of data points will be #lines in merge+1
  int Frow; //row of the furthest pair point
  int point=Rcpp::as<int>(removed); //pull in data to usable form
  int Npairs=((nl-1)*(nl-1)-(nl-1))*.5;
  Rcpp::NumericVector MCC_Size(Npairs); //final MCC sizes  should be nl_C_2- skipped
  int MCi=0;

  for (int a=0; a<nl-1; a++){
    if (a==point-1) {continue;} //skip if removed point
    for(int b=a+1;b<nl;b++){
      if (b==point-1) {continue;} //skip if removed point
      int RowA;
      int RowB;

      // Find the rows of points in pair and find largest one

      for (int i=0;i<(nl-1);i++){ // find the furthest row that contains a member of the pair
        for (int j=0; j<2 ; j++){
          if (Merge(i,j)==-(a+1)) {RowA=i;} //cout<<i<<","<<j<<endl; cout<< ">>"<<Merge(i,j)<<endl;
          if (Merge(i,j)==-(b+1)) {RowB=i;} //cout<<i<<","<<j<<endl; cout<<">>"<<Merge(i,j)<<endl;
        }
      }

      if (RowA>RowB){Frow=RowA;} else {Frow=RowB;}

      vector<int> MCC(nl); // reusable vector to determine size of MCC
      int flagg=-10;
      for (int i=Frow;i<(nl-1);i++){ // skip to the line where the last of a or b appears -need both in subtree
        for (int j=0; j<2 ; j++){
          int results=contains(Merge, Merge(i,j), a , b); // check to see if subtree has both
          if (results==3){ // if subtree has both points
            flagg=1;
            MCC[0]=Merge(i,j);
            for (int k=0; k<nl; k++){//go through the vector from start to finish
              if (MCC[k]>0) {
                int nodePoint=MCC[k]-1; //index that node point is pointing towards
                MCC[k]= Merge(nodePoint,0);
                MCC.insert(MCC.begin()+k, Merge(nodePoint,1));
                if (MCC[k]>0) k--; // if the node point contains a node point in the first slot will revisit that point
              }
            }
          }
          if (flagg==1) break;
        }
        if (flagg==1) break;
        flagg=0;
      }


      if (flagg==0){
        MCC_Size[MCi]=nl-1; // if not a sub tree is whole tree - removed point
      } else {
        //remove all empty spaces in vector and removed point
        MCC.erase(std::remove(MCC.begin(), MCC.end(), 0), MCC.end());
        MCC.erase(std::remove(MCC.begin(), MCC.end(), -point), MCC.end());

        MCC_Size[MCi]=MCC.size();
      }
      MCi++;

      for (int i=0; i<(int)MCC.size(); i++){
      }
    }
  }

  return MCC_Size;

}

RcppExport SEXP rebuild(SEXP margeTable){
  Rcpp::NumericMatrix Merge(margeTable);
  int nl=Merge.nrow()+1; //number of data points will be #lines in merge+1
  int Frow; //row of the furthest pair point
  int Npairs=(nl*nl-nl)*.5;
  Rcpp::NumericVector MCC_Size(Npairs); //final MCC sizes  should be nl_C_2- skipped
  int MCi=0;

  //Rcpp::NumericVector error1(1);

  for (int a=0; a<nl-1; a++){
    for(int b=a+1;b<nl;b++){
      int RowA;
      int RowB;

      // Find the rows of points in pair and find largest one

      for (int i=0;i<(nl-1);i++){ // find the furthest row that contains a member of the pair
        for (int j=0; j<2 ; j++){
          if (Merge(i,j)==-(a+1)) {RowA=i;} //cout<<i<<","<<j<<endl; cout<< ">>"<<Merge(i,j)<<endl;
          if (Merge(i,j)==-(b+1)) {RowB=i;} //cout<<i<<","<<j<<endl; cout<<">>"<<Merge(i,j)<<endl;
        }

      }

      if (RowA>RowB){Frow=RowA;} else {Frow=RowB;}

      vector<int> MCC(nl); // reusable vector to determine size of MCC

      int flagg=0;
      for (int i=Frow;i<(nl-1);i++){ // skip to the line where the last of a or b appears -need both in subtree
        for (int j=0; j<2 ; j++){
          int results=contains(Merge, Merge(i,j), a , b); // check to see if subtree has both

          if (results==3){ // if subtree has both points
            flagg=1;
            MCC[0]=Merge(i,j);
            for (int k=0; k<nl; k++){//go through the vector from start to finish
              if (MCC[k]>0) {
                int nodePoint=MCC[k]-1; //index that node point is pointing towards
                MCC[k]= Merge(nodePoint,0);
                MCC.insert(MCC.begin()+k, Merge(nodePoint,1));
                if (MCC[k]>0) k--; // if the node point contains a node point in the first slot will revisit that point
              }
            }
          }
          if (flagg==1) break;
        }
        if (flagg==1) break;
      }

      if (flagg==0) MCC_Size[MCi]=nl; // if not a sub tree is whole tree - removed point
      else {
        //remove all empty spaces in vector and removed point

        MCC.erase(std::remove(MCC.begin(), MCC.end(), 0), MCC.end());

        MCC_Size[MCi]=MCC.size();
      }
      MCi++;

      for (int i=0; i<(int)MCC.size(); i++){
      }
    }
  }

  return MCC_Size;

}


RcppExport SEXP MLmerge(SEXP ML){

  Rcpp::NumericMatrix MLout(ML);
  int nl=MLout.nrow();

  Rcpp::NumericMatrix Merge(nl-1 , 2);
  Rcpp::NumericVector Height(nl-1); //empty height vector to be returned to R
  Rcpp::NumericVector Order(nl); //empty order vector to be returned to R

  int mergeRow[nl]; // contains last level of merge table to contain item


  // fill created vectors with intializing values
  for (int i=0; i<nl; i++) {
    mergeRow[i]=0;
  }

  int numal=0;
  for (int k=nl-1; k>0; k--){
    int a=0;
    int b=a;
    int tb=0; // this is to hold the post-join cluster number, needs to be translated to get other leaf
    int first=0;
    for (int j=0; j<nl; j++){ // find the first row where there is a difference this will be the join
      //cout<<MLout(j,k)-MLout(j,k-1)<<endl;
      if (MLout(j,k)-MLout(j,k-1)>=1&&first==0){
        a=j;
        tb=MLout(j,k-1);

        first=1;
        break;

      }
    }
    // get b from MLout[,k]
    for (int j=0; j<nl; j++){
      if (MLout(j,k)==tb){
        b=j;
        break;

      }
    }
    //cout<<a<<"\t"<<tb<<"\t"<<b<<endl;


    int rowA=0;
    int rowB=0;

    for(int i=0; i<nl; i++) {
      if (mergeRow[a]>0) rowA=mergeRow[a];
      else rowA=-1;
    }

    for(int i=0; i<nl; i++) {
      if (mergeRow[b]>0) rowB=mergeRow[b];
      else rowB=-1;
    }

    //add the appropriate values to the merge table, a or b already present include +row else -a and -b
    if (rowA<0) Merge(numal,0)=-(a+1);
    else Merge(numal,0)= (rowA);
    if (rowB<0) Merge(numal,1)=-(b+1);
    else Merge(numal,1)= (rowB);

    //rearange pairs so that interior node always on the right

    if (Merge(numal,0)>Merge(numal,1)){
      int placehold=Merge(numal,1);
      Merge(numal,1)=Merge(numal,0);
      Merge(numal,0)=placehold;
    }

    for(int i=0;i<nl;i++){ // update mergeRow
      if (MLout(i,k)==MLout(b,k)) mergeRow[i]=numal+1;
      if (MLout(i,k)==MLout(a,k)) mergeRow[i]=numal+1;
      //cout<<mergeRow[i]<<"\t";
    }
    //cout<<endl<<endl;
    Height[numal]= numal+1;
    numal++;
  }

  vector<int> C_Order(nl);
  for (int i=0; i<2; i++) C_Order[i]=Merge(nl-2,i);  //grab the last line of the merge table -1 for index -1 because n-1 # of joins
  for (int i=0; i<nl; i++){//go through the vector from start to finish
    if (C_Order[i]>0) {
      int nodePoint=C_Order[i]-1; //index that node point is pointing towards
      C_Order[i]= Merge(nodePoint,0);
      C_Order.insert(C_Order.begin()+i, Merge(nodePoint,1));
      if (C_Order[i]>0) i--; // if the node point contains a node point in the first slot will revisit that point
    }
  }
  for (int i=0; i<nl; i++){ Order(i)=C_Order[i]*-1;}



  //combine all pieces into dataframe

  // Rcpp::List pieces= Rcpp::List::create(Merge, Height, Order);
  Rcpp::List pieces= Rcpp::List::create(Merge, Height, Order);

  return (pieces);

}

RcppExport SEXP tripleComp(SEXP data_matrix, SEXP samples, SEXP columns, bool display_progress){//main routine
  //extract the information that is passed to C++ and put it in usable forms

  Rcpp::NumericVector error1(1);

  Rcpp::NumericMatrix data(data_matrix); //data matrix
  Rcpp::NumericVector dim1(samples);
  int lines=dim1[0]; //extract #lines from R data
  Rcpp::NumericVector dim2(columns);
  int features=dim2[0]; //extract # of features from R data
  float ij, jk, ik; // distances within triplets
  //generate all variables that require information from R
  Rcpp::NumericMatrix A(lines, lines); //make R data varaibale
  Rcpp::NumericMatrix DistM(lines, lines); //make R data varaibale
  // initialize Associator matrix with small + value
  fill(A.begin(), A.end(), 0.000001);
  int i,j, k;        //loop indices

  Progress p(lines, display_progress);

  //Make Distance Matrix

for(i=0;i<lines;i++){//loop over samples
  for(j=i+1;j<lines ;j++){
    DistM(i,j)=sqrt(MyD(i,j, features, data));

  }
}
  //Generate Associator Matrix from triplet comparisons

  for(i=0;i<lines-2;i++){//loop over samples
    for(j=i+1;j<lines-1 ;j++){
      ij=DistM(i,j);
      for (k=j+1; k<lines;k++){
        ik=DistM(i,k);
        jk=DistM(j,k);
        if ((ij<ik) && (ij<jk))
          {A(i,j)+=(ik+jk)/ij;}
        else if ((ik<ij) && (ik<jk))
        {A(i,k)+=(ij+jk)/ik;}
        else if ((jk<ij) && (jk<ik))
        {A(j,k)+=(ij+ik)/jk;}
      }
    }

    if (Progress::check_abort() ){ return error1;}
    p.increment(); // update progress

  }
  cout<<endl;


  return A;     //pass the associator matrix back to R
}


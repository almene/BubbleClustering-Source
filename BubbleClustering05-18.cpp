
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

void findminMAX(int&lines, int& features, Rcpp::NumericMatrix& data){//find the distance parameters

  int i,j;     //loop index variables
  double cp;   //comparison variables

  minn=MAXX=MyD(0,1, features, data);  //initialize distances
  for(i=0;i<lines-1;i++)for(j=i+1;j<lines;j++){
    cp=MyD(i,j,features,data); //get next number for comparison
    if(minn>cp)minn=cp;  //new minimum?
    if(MAXX<cp)MAXX=cp;  //new maximum?
  }

}

void memberate(int v,double r, int& lines, int& features, Rcpp::NumericMatrix& data, int* member){//find membership

  int i;  //loop index

  Z=0;  //number of members
  //find sphere members
  for(i=0;i<lines;i++){
    double dist=MyD(v,i, features, data);
   {if(dist<=r){member[i]=1;Z++;} else member[i]=0;}}

  reward=1.0/((double)Z); //compute the reward
}

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




// Main R program
RcppExport SEXP bubble(SEXP data_matrix, SEXP samples, SEXP columns, SEXP bubbleLevel, SEXP bubbleSize, bool display_progress){//main routine
  //extract the information that is passed to C++ and put it in usable forms

  Rcpp::NumericVector error1(1);

  Rcpp::NumericMatrix data(data_matrix); //data matrix
  Rcpp::NumericVector dim1(samples);
  int lines=dim1[0]; //extract #lines from R data
  Rcpp::NumericVector dim2(columns);
  int features=dim2[0]; //extract # of features from R data

  //generate all variables that require information from R
  Rcpp::NumericMatrix A(lines, lines); //make R data varaibale

  // initialize Associator matrix with small + value
  fill(A.begin(), A.end(), 0.000001);

  int member[lines];             //characteristic function
  int i,j;        //loop indices
  double radius;  //critical radius
  Rcpp::NumericVector dim3(bubbleLevel);
  Rcpp::NumericVector dim4(bubbleSize);

  int numB=dim3[0]; //extract Bubble level from R data
  long double sizeB=dim4[0]; // extract Bubble max size

  int SAM=lines*pow(10, numB);

  findminMAX(lines, features, data);   //find the minimum and maximum values
  MAXX-=minn;     //normalize the maximum
  //cout<<MAXX<<endl;
  Progress p(SAM, display_progress);

  for(i=0;i<SAM;i++){//loop over samples
    double randNum=unif_rand();
    radius=randNum*MAXX*sizeB+minn; //generate a random critical radius
    //if(i%100==0){cout<<radius<<endl;}
    j=i%lines;          //generate core vertex
    memberate(j,radius, lines, features, data, member);        //find the membership fuction
    Associate(lines, member, A);                //update the associator
    //print progress bar and check for abort
    if (Progress::check_abort() ){ return error1;}
    p.increment(); // update progress

  }
  cout<<endl;


  return A;     //pass the associator matrix back to R
}

RcppExport SEXP bubbleOrder(SEXP data_matrix, SEXP samples, SEXP columns, SEXP bubbleLevel, bool display_progress){//main routine
  //extract the information that is passed to C++ and put it in usable forms

  Rcpp::NumericVector error1(1);

  Rcpp::NumericMatrix data(data_matrix); //data matrix
  Rcpp::NumericVector dim1(samples);
  int lines=dim1[0]; //extract #lines from R data
  Rcpp::NumericVector dim2(columns);
  int features=dim2[0]; //extract # of features from R data

  //generate all variables that require information from R
  Rcpp::NumericMatrix A(lines, lines); //make R data varaibale

  // initialize Associator matrix with small + value
  fill(A.begin(), A.end(), 0.1);

  int member[lines];             //characteristic function
  int i,j;        //loop indices
  double radius;  //critical radius
  Rcpp::NumericVector dim3(bubbleLevel);

  int numB=dim3[0]; //extract Bubble level from R data
  int divisions=pow(10, numB);
  int SAM=lines*divisions;

  findminMAX(lines, features, data);   //find the minimum and maximum values
  MAXX-=minn;     //normalize the maximum
  //cout<<MAXX<<endl;
  Progress p(SAM, display_progress);

  for(i=0;i<divisions;i++){//loop over samples
    double sizeB=(i%divisions)*1.0/divisions;

    radius=MAXX*sizeB+minn; //generate a random critical radius
    //if(i%10==0){cout<<sizeB<<"\t"<<radius<<endl;}
    //if(i%100==0){cout<<radius<<endl;}
    for (int s=0;s<lines;s++){
      j=i%lines;//generate core vertex
    memberate(j,radius, lines, features, data, member);        //find the membership fuction
    Associate(lines, member, A);                //update the associator
    }
    //print progress bar and check for abort
    if (Progress::check_abort() ){ return error1;}
    p.increment(); // update progress

  }
  cout<<endl;


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

    //add the appropriate values to the merge table, a or b alread present include +row else -a and -b
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
      for (int i=Frow;i<(nl-2);i++){ // skip to the line where the last of a or b appears -need both in subtree
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

      if (flagg==0) MCC_Size[MCi]=nl-1; // if not a sub tree is whole tree - removed point
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


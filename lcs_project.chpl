use LinearAlgebra;
use IO;
use Math;
var string1="gatttatgcagg";
var string2="tcaggatt";

//correct ans="tcagg"

string2="ab";
string1="ab";
//correct ans=ab
string2="abc";
string1="ab";
//correct  ans="ab"
string2="abc";
string1="abc";
//correct ans="abc"
string2="abcd";
string1="abcd";
//correct ans="abcd"
string2="abcdcdb";
string1="abcdbe";
//correct ans="abcdb"
string2="spaceship";
string1="spakehip";
//correct ans="spaehip"
string2="abcdgh";
string1="aedfhr";
//correct ans="adh"
//string2="abcdghdr";
//string1="aedfhrd";
//correct
//string2="workattech";
//string1="branch";
//correct
//string2="helloworld";
//string1="playword";
//correct;
//string2="hello";
//string1="hello";
//correct
//string2="abc";
//string1="def";
//correct
//string2="abba";
//string1="baab";
//correct ans="ab"
//string2="aabbccbbaa";
//string1="abcb";
//correct ans="abcb"
//string2="abbababbbba";
//string1="bbaaaabba";
//correct ans="bbaabba";
//string2="bcdaacd";
//string1="acdbac";
//correct ans="cdac"
//string1="abcadb";
//string2="bfdb";
//correct; ans="bdb"
//string1="azbzc";
//string2="zq";
//correct; ans="z"
//string1="adf";
//string2="h";
//correct ans=""
string2="absdfesccbad";
string1="cbdaserthpsd";
//correct ans="bdesd";
//string2="bcdsaecvbzasdrf";
//string2="aazxdcvbsbaxzds";
//correct ans="cbsd"
if(string1.size>string2.size){
  var u= string2;
  string2=string1;
  string1=u;
  
}
const infin:int=9999;
const m:int=string1.size;
const n:int=string2.size;
var mat:[0..m][0..n] int;
//var j1s:[0..m-1] int;
forall i in 0..m-1{
    
  //  var maxi=1000*1000;
    forall j in 0..n-1{
      if string1[i]==string2[j]{
  //      maxi=min(maxi,j);
        mat[i][j]=1;
      }
    }
//    j1s[i]=infin;
}
proc printMatrix(matrix) {
  var nn=matrix.size;
  var mm=matrix[matrix.dim(0).low].size;
  for i in matrix.dim(0).low..matrix.dim(0).high{
    for j in matrix[matrix.dim(0).low].dim(0).low..matrix[matrix.dim(0).low].dim(0).high{
      write(matrix[i][j], " ");
    }
    writeln();
  }
}

proc prefsum(list:[0..n] int){
    var exp: int = 1;
    var expm1, expnot: int;

    while (exp < n+1) {
        expm1 = exp - 1;
        expnot = ~exp;
        forall j in 0..n {
          if (j&exp != 0){
               list[j] = list[j]+list[j&expnot|expm1];
           //   prefix[j] = min(prefix[j],prefix[j&expnot|expm1]);
          }
        }
      exp = exp << 1;
    //  writeln(prefix);
  }
//  writeln(list);
}
var jk:[0..m][0..n] int;

var jkimax:[0..m] int;
forall i in 0..m{
  var c=mat[i];
  prefsum(c);
//  writeln(c);
  forall j in 0..n{
    if (((j==0)&&(c[0]==1))||((j>0)&&(c[j]!=c[j-1]))){
      jk[i][c[j]-1]=j;
      jkimax[i]=max(jkimax[i],j+1);
    }

  }
}
//writeln(jk[2]);
//writeln(mat[0]);
var dgh:[0..m][0..0][0..n] int;
var dgs:[0..m][0..n][0..0] int;
forall i in 0..m{
  dgh[i][0][0]=jk[i][0]+1;
  forall j in 1..n{
    if (jk[i][j]>0)&&(jk[i][j-1]<n){
      dgh[i][0][jk[i][j-1]+1]=jk[i][j]-jk[i][j-1];
    }
  }
  prefsum(dgh[i][0]);
  forall j in 0..n{
    if(j>=jkimax[i]){
      dgh[i][0][j]=infin;
    }
  }
}
forall i in 0..m{
  forall j in 0..n{
    dgs[i][j][0]=dgh[i][0][j];
  }
}
//writeln((dgs[1]));
var rowdg:[0..2*m][0..n][0..1] int;
forall i in 0..m{
  forall j in 0..n{
    rowdg[i][j][1]=dgs[i][j][0];
    rowdg[i][j][0]=j;

  }
}
forall i in m+1..2*m{
  forall j in 0..n{
    rowdg[i][j][0]=j;
    rowdg[i][j][1]=infin;
  }
}
/*var firstone:[0..n][0..3*m];

for i in 0..m{
  printMatrix(rowdg[i]);
  writeln();
}*/
/*writeln("dgh[1][0]");
for i in 0..n-1{
  write(dgh[0][0][i],"  ");

}
writeln(" ");*/
//writeln(rowdg[0][0]);

//printMatrix(dgs[0],2,n);

//var dgu=dgs[0];
//var x = dgu.domain.dim(0).size;
//var k=2*x;
//writeln(x,k);
//writeln(size(dgs[1][0]));
//writeln(dgs[0].domain.dim(0).size);
//writeln(dgs[0].domain.dim(0..1).size);

proc find_cell(dgu,dgl,vertex,p,q){
  if p>q{
    return infin;
  }
 /* printMatrix(dgu);
  printMatrix(dgl);
  writeln(num_bre);
  writeln(p);
  writeln(q);*/
  //var x=dgu.dim(1).size();
  //the combination of them is a k*n matrix such that the minima of column j is the index of the num_bre th breakout
  // of vertex j.
  //  var f=dgu.domain.dims();
 //   x=f(0).size;
     /*   if(q<p){
            return infin;
          }
          else{*/
            if  ((p<dgu[0].size)&&(dgu[vertex][p]<n) && (q-p<dgl[0].size)&&(q-p>=0)){
              return dgl[dgu[vertex][p]][q-p];
    //  writeln("i== ",i," j== ",j, "  dgu[j][i]=  ",dgu[j][i],  "  dgl[dgu[j][i]][k-i]= ",dgl[dgu[j][i]][k-i]);
              }
              else{
                return infin;
              }
    //  writeln(ans[0][0]);
  }
//}

// Parallelize this function
proc findMinIndex( dgu,dgl,vertex:int, col: int, top: int, bottom: int) {
  //writeln("dgu =",dgu);
  //writeln("dgl =",dgl);
  var listsize = bottom - top + 1;
  var exp: int = 1;
  var expm1, expnot: int;

  var prefix: [0..listsize-1] int;
  var minIndex: [0..listsize-1] int;

  forall i in 0..listsize-1{
    prefix[i] = find_cell(dgu,dgl,vertex,i+top,col);//matrix[i + top,col];
    minIndex[i] = i;
  }

  //writeln("Input to findMinIndex: ", col," ", top," ", bottom," prefix: ",prefix);

  while (exp < listsize) {
    expm1 = exp - 1;
    expnot = ~exp;
    forall j in 0..listsize-1 {
      if (j&exp != 0){
        // prefix[j] = prefix[j]+prefix[j&expnot|expm1];
        // prefix[j] = min(prefix[j],prefix[j&expnot|expm1]);
        if (prefix[j&expnot|expm1] <= prefix[j]) {
          prefix[j] = prefix[j&expnot|expm1];
          minIndex[j] = minIndex[j&expnot|expm1];
        }
      }
    }
    exp = exp << 1;


}
  //writeln("prefix: ", prefix," minIndex: ", minIndex);
  return minIndex[listsize-1]+top;
}
var mins:[0..4*(n+1)*(m+1)]int;
proc findColMins(dgu,dgl,vertex:int, left: int, right: int, top: int, bottom: int, mins,firstind:int) {

  
  var cols = right - left + 1;
  if(cols < 1){
    return;
  }
  //writeln("num =",num,"  left= ", left, "right = ",right);
  var midCol = ceil((right+left)/2): int;
  var minIndex = findMinIndex(dgu,dgl,vertex:int, midCol, top, bottom);
  //writeln(firstind," ",firstind+midCol);
  if midCol==5{
   // writeln(" the minima of column ",midCol," is index ",findMinIndex(dgu,dgl,vertex:int, midCol, top, bottom)," but variable MinIndex is ",minIndex );
  }
  mins[firstind+midCol-left] = minIndex;
  //writeln("firstind= ",firstind," midcol= ",midCol," left= ",left," right= ",right," mins index= ",firstind+midCol-left,"  minIndex= ",minIndex," top= ",top," bottom = ",bottom);
  if find_cell(dgu,dgl,vertex,minIndex,midCol)!=infin{

   // cobegin{ //Parallel step
      findColMins(dgu,dgl,vertex:int, left, midCol-1, top, minIndex, mins,firstind); //left submatrix
      findColMins(dgu,dgl,vertex:int, midCol+1, right, minIndex, bottom, mins,firstind+midCol-left+1); //right submatrix
  //  }
  }
  else{
    findColMins(dgu,dgl,vertex:int, left, midCol-1, top, bottom, mins,firstind);
  }
}

/*var fortwo:[0..n-1][0..4] int;
proc compute2(dgu,dgl, left: int, right: int, top: int, bottom: int, mins,x:int){
  forall i in 0..x{
    findColMins(dgu,dgl,i,1,n-1,1,1,mins[i]);
    forall j in 0..n-1{
      if dgu[0].size>i{
        fortwo[j][i]=min(dgu[j][i],min(dgl[j][i],find_cell(dgu,dgl,i,mins[i][j],j)));
        if(j==3&&i==1){
          writeln("dgu[j][i] = ",dgu[j][i],"  dgl[j][i]= ",dgl[j][i]," tarkibi= ",find_cell(dgu,dgl,i,mins[j],j)," mins[j]=  ",mins[j]);
        }
      }
      else{
        if((j==0)&&(i==2)){
         // writeln(" for j=0 and i=2: ",mins[i][j]);
        }
        fortwo[j][i]=find_cell(dgu,dgl,i,mins[i][j],j);
      }

    }

  }
}*/
//writeln("minindex[1]");
//writeln(findMinIndex(rowdg[6],rowdg[7],2,0,0,1));
//findColMins(rowdg[6],rowdg[7],2,1,n-1,1,1,mins[0]);
//writeln("mins");
/*for i in 0..n-1{
  write(mins[0][i]," ");
}*/
//writeln();
//writeln(find_cell(rowdg[6],rowdg[7],2,0,2));
//compute2(rowdg[6],rowdg[7],0,1,0,n-1,mins,2);
//printMatrix(fortwo);

var javab:[0..50][0..n][0..4*m+1] int;
forall k in 0..50{
  forall i in 0..n{
    forall j in 0..4*m+1{
        javab[k][i][j]=infin;
    }
  }
}
forall i in 0..m-1{
  var sot=2*i;
  //printMatrix(rowdg[i]);
  forall j in 0..n{
    javab[0][j][sot]=rowdg[i][j][0];
    javab[0][j][sot+1]=rowdg[i][j][1];

  }
}
//printMatrix(javab);

/*proc test(dgu, dgl,vertex,num_bre){
//  var x = dgu.domain.dim(0).size;
    
    
    for i in 0..num_bre{
      for j in 0..n{
    //finding the ith breakout of the jth vertex

            write(find_cell(dgu,dgl,vertex,i,j)," ");
    //  writeln(ans[0][0]);
          }
      writeln();

    }
}*/
var nummat=0;
var lastbre=-1;
//writeln("matrix number 0 :");
//printMatrix(javab[0]);
proc computedg(left: int, right: int, top: int, bottom: int, mins){ 
  var x=2;
  var intmat=0;
  while((x)<=2*m-1){
    {
     // writeln("x = ",x);
      intmat=intmat+1;
      lastbre=x;
      nummat=intmat;
      var cp:[0..n][0..4*m+1] int ;
      if(intmat>0){
        cp=javab[intmat-1];
      }
      /*if(x>m){
        x=m;
      }*/
      forall f in 0..m/x {
        
       // writeln("n= ",n);
        forall i in 0..n{
          //writeln("i= ",i," x= ",x);
          var dgu:[0..n][0..((x+2)/2)-1] int;
          var dgl:[0..n][0..((x+2)/2)-1] int;
          if(x==2){
            dgu=rowdg[2*f];
            dgl=rowdg[2*f+1];
          //  writeln("reached step 0");
          }
          else{
           // writeln("x = ",x," f= ",f,"  ",(f*(x+2))+(x+2)/2-1);
           // writeln(f*(x+2)," ",(f*(x+2))+(x+2)/2-1);
            [u in 0..n][v in 0..(x+2)/2-1]dgu[u][v]=cp[u][(f*(x+2))+v];
            [u in 0..n][v in 0..(x+2)/2-1]dgl[u][v]=cp[u][f*(x+2)+(x+2)/2+v];
            
          }
        //  writeln("reached step 1");
        /*  writeln("nummat= ",nummat," i= ",i,"\n"," dgu= ");
          printMatrix(dgu);
          writeln();
          writeln("dgl= ");
          printMatrix(dgl);
          writeln();*/
       // writeln(i);
          if x==2{
           // writeln(" reached step2! i= ",i);
            findColMins(dgu,dgl,i,0,x,0,x/2,mins,(f*n+i)*(x+1));
          }
          else{
            findColMins(dgu,dgl,i,0,x,0,x/2,mins,(f*n+i)*(x+1));
          }
         // writeln("first stop x = ",x);
          forall j in 0..x{
          //  writeln("j = ",j," (f*n+i)*(x+1)+j = ",(f*n+i)*(x+1)+j," mins.size= ", mins.size," i = ",i,"   f=  ",f);
            if dgu[0].size>j{
             // writeln( i, "=i ", j, " = j", " index=  ",(f*n+i)*(x+1)+j);
              javab[intmat][i][f*(x+1)+j]=min(dgu[i][j],min(dgl[i][j],find_cell(dgu,dgl,i,mins[(f*n+i)*(x+1)+j],j)));
              if((f==0)&&(j==4)&&(i==0)&&(x==4)){
               // writeln(find_cell(dgu,dgl,i,mins[i][j],j));
            //    writeln("x==",x,"  dgu[j][i] = ",dgu[j][i],"  dgl[j][i]= ",dgl[j][i]," tarkibi= ",find_cell(dgu,dgl,i,mins[(f*n+i)*(x+1)+j],j),"  mins[i][j]= ",mins[i]);
            //    printMatrix(dgu);
            //    writeln();
              //  printMatrix(dgl);
              }
            }
            else{
              /* if((f==0)&&(j==5)&&(i==0)&&(x==8)){
               // writeln(find_cell(dgu,dgl,i,mins[i][j],j));
                writeln("x==",x," tarkibi= ",find_cell(dgu,dgl,i,mins[(f*n+i)*(x+1)+j],j)," index of min= ",mins[(f*n+i)*(x+1)+j]);
                printMatrix(dgu);
                writeln();
                printMatrix(dgl);
                writeln();
                test(dgu,dgl,i,5);
                writeln();
                for l in 0..x{
                  write(mins[(f*n+i)*(x+1)+l]," ");
                  
                }
                writeln();
                for l in 0..x{
                  write(findMinIndex(dgu,dgl,0,l,0,x)," ");
               
                }
                writeln();
                writeln("lets recheck");
                findColMins(dgu,dgl,i,0,x,0,x/2,mins,0);
                for l in 0..x{
                  write(mins[l]," ");
                  //writeln();
                }
                writeln();
                coforall l in 0..x{
                  write(l," : ",findMinIndex(dgu,dgl,0,l,0,x),"   ");
               
                }
                writeln();

              }*/
              //writeln(findMinIndex(dgu,dgl,0,5,0,x));
              //find_cell(dgu,dgl,i,mins[(f*n+i)*(x+1)+j],j);
              javab[intmat][i][f*(x+1)+j]=find_cell(dgu,dgl,i,mins[(f*n+i)*(x+1)+j],j);
            }
          }

          }

        }
    }
    //if(x==4){
    
   // writeln("matrix number ",intmat," :");
   // printMatrix(javab[intmat]);
    //}
   // writeln();
    x=x*2;
  }
}

  
var ans:[0..3][0..n-1] int;

//var wshift:[0..1][0..n-1] int;


//var z;
//test(rowdg[z],rowdg[z+1],2);
//printMatrix(mins);
//writeln("test");
//printMatrix(ans);
//findColMins(rowdg[z],rowdg[z+1],2,0,n-1,0,1,mins[0]);
/*for i in 0..n-1{
  write(mins[0][i]," ");
}*/
//writeln();
//for i in 0..n-1{
 // writeln(findMinIndex(rowdg[6],rowdg[7],2,i,0,1));
//}
/*for z in 0..m{
  writeln("For rows ",z," and ",z+1,":");
  raw(rowdg[z],rowdg[z+1],2);
  test(rowdg[z],rowdg[z+1],2);
}*/

 // printMatrix(ans);
 // writeln();
 //printMatrix(wshift);
 // printMatrix((rowdg[z]));
 // writeln();
 // printMatrix(rowdg[z+1]);
 // writeln();
 //printMatrix(rowdg[z+1]);
/*for z in 0..6{
  writeln("z= ",z);
  test(rowdg[z],rowdg[z+1],2);*/
//writeln(find_cell(rowdg[1],rowdg[2],2,1,9,2));
 // for j in 0..n-1{
//    writeln(findMinIndex(rowdg[z],rowdg[z+1],2,j,0,1,2));
 // }
 // findColMins(rowdg[z],rowdg[z+1],1,0,n-1,0,1,mins);
  //writeln(mins);
  //for i in 0..3{
   // writeln();
  //}
//}
/*for i in 0..5{
 // writeln();
}*/
//test(rowdg[0],rowdg[1],2);
//writeln();
computedg(0,n,0,2,mins);
lastbre=lastbre/2;
//printMatrix(javab);
//writeln();

/*var dgu:[0..n-1][0..2] int;
if(intmat>0){
    [u in 0..n-1][v in 0..2]dgu[u][v]=javab[intmat][u][v];
}
var dgl:[0..n-1][0..2]int ;
if(intmat>0){
  [u in 0..n-1][v in 0..2]dgl[u][v]=javab[u][v+3];
}*/
//printMatrix(dgu);
//writeln();
////printMatrix(dgl);
//writeln();
//test(dgu,dgl,2);
//writeln();
//printMatrix(ans);
//printMatrix(javab);
//writeln();
//findColMins(dgu,dgl,2,0,1,0,n-1,mins[2]);
//for i in 0..n-1{
//  printMatrix(mins[2]);
//  write(mins[2][i],"  ");
//}
//writeln();
forall i in 0..n-1{
 // writeln("for vertex i= ",i);
 // test(rowdg[0],rowdg[1],i,2);
 // writeln();
  findColMins(rowdg[0],rowdg[1],i,0,2,0,1,mins,0);
 // for j in 0..2{
  //  write(mins[i][j]," ");
 // }
  //writeln();
 // writeln();
}
proc findlast(){
  var left=0;
  var right=m;
  
  while(1){
    var mid=(left+right)/2;
    //writeln(left," ",right," ",mid);
    if((javab[nummat][0][mid]!=infin)&&((mid==right)||(javab[nummat][0][mid+1]==infin))){
      return mid;
    }
    if(javab[nummat][0][mid]==infin){
      right=mid-1;
    }
    else{
      left=mid+1;
    }
  }
  return left;
}
//writeln("findlast = ",findlast());
var cross:[0..m] int;
proc cross_finder(leftgr,rightgr,topgr,bottomgr,breakout,mat_num){
  
  if(rightgr-leftgr<breakout){
    return;
  }
  if((leftgr>rightgr)){
    return;
  }
  if((topgr>bottomgr)){
    return;
  }
  if(breakout<0){
    return;
  }
  if topgr==bottomgr{
    cross[topgr]=leftgr;
  //  writeln("kind 0 :  cross[",topgr,"]= ",leftgr);
    return;
  }
   
  if (leftgr==rightgr)||(breakout==0){
    forall i in topgr+1..bottomgr{
      cross[i]=leftgr;
    //  writeln("for kind 1:  ,leftgr= ",leftgr, " rightgr= ",rightgr,"  topgr = ",topgr, " ", "bottomgr= ",bottomgr, "  number of breakout= ",breakout," mat_num= ",mat_num," number of col = ",num_of_col);
    //  writeln("kind 1 : cross[",i,"]= ",leftgr);
    }
    return;
  }
  if (breakout==1)&&(bottomgr-topgr==1){
    cross[topgr]=leftgr;
    cross[bottomgr]=javab[0][leftgr][2*topgr+1];
  //  writeln("kind 2 :cross[",topgr,"] = ",leftgr);
  //  writeln("kind 3 :cross[",bottomgr,"] = ",javab[0][leftgr][2*topgr+1]);
    //if(bottomgr==4 && javab[0][topgr][2*topgr+1]==8){
    //    writeln("for kind 3:  ,leftgr= ",leftgr," leftmat= ",leftmat, " rightgr= ",rightgr," rightmat =  ",rightmat,"  topgr = ",topgr, " topmat= ",topmat," ", "bottomgr= ",bottomgr, " bottommat= ",bottommat,"  number of breakout= ",breakout," mat_num= ",mat_num," number of col = ",num_of_col);
    //}
    return;
  }
 
  var topmat=0;
  var bottommat=n;
  var the_related_matrix=logBasePow2(bottomgr-topgr-1,1);
  var num_of_col= 1<<(the_related_matrix);
  var u=(topgr)/(num_of_col);
  var leftmat= topgr+u;
  var rightmat=leftmat+num_of_col;

  var midcol=(leftgr+rightgr)/2;
  
  var which=-1;
  var dgu:[topmat..bottommat][leftmat..rightmat] int;
  
  if(mat_num>0){
    [u in topmat..bottommat][v in leftmat..rightmat]dgu[u][v]=javab[the_related_matrix][u][v];
  }
  
  var dgl:[topmat..bottommat][rightmat+1..2*rightmat-leftmat+1]int ;
  //writeln(dgl.dim(0).low," ",dgl.dim(0).high," ",dgl[6].dim(0).low," ",dgl[6].dim(0).high);
  if(mat_num>0){
      [u in topmat..bottommat][v in rightmat+1..2*rightmat-leftmat+1]dgl[u][v]=javab[the_related_matrix][u][v];
  }
  
  var a:[0..2*n] int;

  var midrow=(topgr+bottomgr+1)/2;
  //if topgr==0{
      var lg=logBasePow2(midrow,1);
      if((1<<lg)!=midrow){
        lg=lg+1;
      }
      if(topgr+(1<<lg)<bottomgr){
        midrow=topgr+1<<(lg);
      }
  //}
  /*writeln("leftgr= ",leftgr," leftmat= ",leftmat, " rightgr= ",rightgr," rightmat =  ",rightmat,"  topgr = ",topgr, " topmat= ",topmat," ", "bottomgr= ",bottomgr, " bottommat= ",bottommat," ", " midcol= ",midcol,"  number of breakout= ",breakout," mat_num= ",mat_num," number of col = ",num_of_col);
  writeln("dgu = ");
  printMatrix(dgu);
  writeln();
  writeln("dgl = ");
  printMatrix(dgl);
  writeln();
  writeln();
  writeln("leftgr= ",leftgr," leftmat= ",leftmat, " rightgr= ",rightgr," rightmat =  ",rightmat,"  topgr = ",topgr, " topmat= ",topmat," ", "bottomgr= ",bottomgr, " bottommat= ",bottommat," ", " midcol= ",midcol,"  number of breakout= ",breakout," mat_num= ",mat_num," number of col = ",num_of_col);*/
  //readLine();
  var f=leftmat+breakout;
  forall i in leftmat..min(rightmat,f) {
    
    if (dgu[leftgr][i]<=n){
      /*  if((leftgr==0)&&(midrow==1)){
          writeln("first one passed");
        }*/
        if((f-i+dgl[0].dim(0).low<=dgl[dgu[leftgr][i]].dim(0).high)){
         /* if((leftgr==0)&&(midrow==1)){
            writeln("second one passed");
          }*/
      
          /*if((leftgr==0)&&(midrow==1)){
            writeln("for leftgr= 0   i== ",i," dgu[leftgr][i]== ",dgu[leftgr][i]," f-i+dgl[dgl.dim(0).low].dim(0).low = ",f-i+dgl[dgu[leftgr][i]][f-i+dgl[dgl.dim(0).low].dim(0).low],"  dgl[dgu[leftgr][i]][f-i+dgl[0].dim(0).low] = ", dgl[dgu[leftgr][i]][f-i+dgl[0].dim(0).low]," but should be equal to ",rightgr);
          }*/
          if((dgl[dgu[leftgr][i]][f-i+dgl[0].dim(0).low]==rightgr)){
           // if((leftgr==11)&&(midrow==6)){
          //    writeln("third one passed! i= ",i-leftmat);
           // }
            a[i]=i-leftmat;
           // writeln("dgl[dgu[",leftgr,"][",i,"]][",f-i,"]= ",javab[mat_num][leftgr][breakout],")");
         }
      }
    }
  }
 // if((leftgr==0)&& (midrow==4)){

    
 // }
  
  which = max reduce a;
  /*if a[which]==0{
    which=0;
  }*/
  //writeln(a);
  //which=a[which];
  //writeln("which = ",which);
  //begin{
   //   writeln(" for line 579:");
   //   writeln("leftgr= ",leftgr," leftmat= ",leftmat, " rightgr= ",rightgr," rightmat =  ",rightmat,"  topgr = ",topgr, " topmat= ",topmat," ", "bottomgr= ",bottomgr, " bottommat= ",bottommat," ", " midcol= ",midcol,"  number of breakout= ",breakout," mat_num= ",mat_num," number of col = ",num_of_col);
  
      cross[midrow]=dgu[leftgr][which+leftmat];
 // }
  /*writeln("kind 4:  cross [",midrow,"]= ",dgu[leftgr][which+leftmat]);
  //readLine();
  //if(midrow==2){
      writeln("for kind 4:  leftgr= ",leftgr," leftmat= ",leftmat, " rightgr= ",rightgr," rightmat =  ",rightmat,"  topgr = ",topgr, " topmat= ",topmat," ", "bottomgr= ",bottomgr, " bottommat= ",bottommat,"  number of breakout= ",breakout," mat_num= ",mat_num," number of col = ",num_of_col," which= ",which);
      writeln();
    writeln("a = ",a);
    writeln("dgu = ");
    printMatrix(dgu);
    writeln();
    writeln("dgl = ");
    printMatrix(dgl);
    writeln(); */  

  //}
  //writeln();
  //begin{
   // if mat_num==nummat-1{
    /*cobegin{
      writeln("the below is for the left part of the recursive function!!");
      writeln();
      writeln();
      writeln("it is : leftgr= ",leftgr," leftmat= ",leftmat, " rightgr= ",dgu[leftgr][which+leftmat]," rightmat =  ",rightmat/2,"  topgr = ",topgr, " topmat= ",topmat," ", "bottomgr= ",midrow, " bottommat= ",bottommat," ", "  number of breakout= ",which," mat_num= ",mat_num-1," number of col = ",num_of_col/2);
      writeln();
      writeln("ended!");
      writeln();
    }*/
   // }
   // cobegin{
        cross_finder(leftgr,dgu[leftgr][which+leftmat],topgr,midrow,which,mat_num-1);
   // if mat_num==nummat-1{
    /*cobegin{
      writeln("the below is for the right part of the recursive function!!");
      writeln();
      writeln();
      writeln("it is : leftgr= ",dgu[leftgr][which+leftmat]," leftmat= ",dgu[leftmat][which+leftmat]/(num_of_col+1), " rightgr= ",rightgr," rightmat =  ",dgu[leftmat][which+leftmat]/(num_of_col+1)+num_of_col,"  topgr = ",midrow, " topmat= ",topmat," ", "bottomgr= ",bottomgr, " bottommat= ",bottommat," ", "  number of breakout= ",breakout-which," mat_num= ",mat_num-1," number of col = ",num_of_col/2);
      writeln();
      writeln("ended!");
      writeln();
    }*/
   // }
      cross_finder(dgu[leftgr][which+leftmat],rightgr,midrow,bottomgr,breakout-which,mat_num-1);
  //  }
  }

//}
//writeln("lastbre= ",lastbre);
cross_finder(0,javab[nummat][0][findlast()],0,m,findlast(),nummat);

var common:[0..m-1] int;
forall i in 1..m{
  if (cross[i]>0)&&(string1[i-1]==string2[cross[i]-1])&&(cross[i]!=cross[i-1]){
    common[i-1]=1;
  }
}
//writeln(common);
for i in 0..m-1{
  if(common[i]){
    write(string1[i]);
  }
}
writeln();
//printMatrix(rowdg[0]);
//test(rowdg[0],rowdg[1],1,2);
//compute2(dgu,dgl,0,4,0,n-1,mins,4);
//printMatrix(fortwo);

//writeln(find_cell(rowdg[0],rowdg[1],1,mins[1],1));
//var matrices :[3*n][m][n] int;
//}
/*const m=7;
const n=13;
var input:[1..m] (int,int)=[(1,1),(2,10),(3,10),(4,10),(5,11),(6,12),(7,13)];*/

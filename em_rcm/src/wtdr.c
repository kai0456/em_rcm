#include<R.h>
#include<Rmath.h>
#include<math.h>

void wtdr(double *x, int *d, double *weight, double *result)
{
 int nr=d[0];
 int nc=d[1];
 int i,j,m;
 double data[nr][nc],s=0.0,w[nr],ranks[nr][nc];
 double norm,invnorm,temp[nc];

 /*******
  Initialize data[][] and ranks[][] as nr*nc matrix
********/
 for (j=0;j<nc;j++){
     for (i=0;i<nr;i++){
         data[i][j]=x[j*nr+i];
         ranks[i][j]=0.0;
     }
 }
 
 for (i=0;i<nr;i++)
     s+=weight[i];
 for (i=0;i<nr;i++)
     w[i]=weight[i]/s;

 for (i=0;i<nr;i++){
     for (j=0;j<nr;j++){
         norm=0.0;
         for (m=0;m<nc;m++){
             temp[m]=data[i][m]-data[j][m];
             norm+=temp[m]*temp[m];
         }
         norm=sqrt(norm);
         (norm==0) ? (invnorm=norm) : (invnorm=1/norm);
         for (m=0;m<nc;m++){
             temp[m]=temp[m]*invnorm*w[j];
             ranks[i][m]+=temp[m];
         }
     }
 }

 for (j=0;j<nc;j++){
     for (i=0;i<nr;i++){
         result[j*nr+i]=ranks[i][j];
     }
 }
}



#include<stdio.h>
#include<stdlib.h>
#include<math.h>

//--------------------定数--------------------//

#define Z ( 2048 )
#define TMAX ( 5.5 )
#define RANGE_CHECK( x, xmin, xmax ) ( x = ( x < xmin ? xmin : ( x < xmax ?  x : xmax ) ) );


//--------------------関数--------------------//

double* make_vector( int n ); // 領域の確保(ベクトル)
double** make_matrix( int ROW, int COL ); // 領域の確保(行列)
void connect( int N, double *x ); // 閉曲線のつなぎ
void connect_double( int N, double* x, double *y ); // つなぎ2個

double ip( double x1, double x2, double y1, double y2 ); // 内積
double DIST( double x1, double x2, double y1, double y2 ); // |y - x|
double DET( double x1, double x2, double y1, double y2 ); // 行列式
void initial_condition( int N, double *x1, double *x2 );
void quantities( double t, int N, double *x1, double *x2, double *l, double *t1, double *t2, double *n1, double *n2, double *T1, double *T2, double *N1, double *N2, double *nu, double *phi, double *kappa );

//--------------------main--------------------//

int main(void){
  
  //---------------------変数--------------------//
  
  int i,j,z,N = 128;
  double t,dt = 0.1 / ( N * N );
  double L,A,L_tmp = 2 * M_PI;
  
  double *X1,*X2;
  double *t1,*t2;
  double *n1,*n2;
  double *l,*nu,*phi;
  double *T1,*T2;
  double *N1,*N2;
  double *kappa;
  X1 = make_vector(Z + 2);
  X2 = make_vector(Z + 2);
  t1 = make_vector(Z + 2);
  t2 = make_vector(Z + 2);
  n1 = make_vector(Z + 2);
  n2 = make_vector(Z + 2);
  l = make_vector(Z + 2);
  nu = make_vector(Z + 2);
  phi = make_vector(Z + 2);
  T1 = make_vector(Z + 2);
  T2 = make_vector(Z + 2);
  N1 = make_vector(Z + 2);
  N2 = make_vector(Z + 2);
  kappa = make_vector(Z + 2);

  char file[5000];
  FILE *fp;

  t = 0.0;
  z = 0;

  
  //--------------------初期条件--------------------//
  
  initial_condition(N,X1,X2);

  quantities(t,N,X1,X2,l,t1,t2,n1,n2,T1,T2,N1,N2,nu,phi,kappa);
  
  sprintf(file, "./data/yoko_kuro%06d.dat", z);
  fp = fopen(file, "w");

  //  measure(t,N,X1,X2,&L,&A);
  
  for( i = 0; i <= N; i++ ){
    
    //printf("%f %f\n", x[i], y[i]);
    fprintf(fp, "%f %f %f\n", X1[i], X2[i], nu[i]);
    
  }
  
  fclose(fp);
  
  printf("z = %d\n", z);

  
  //--------------------回す--------------------//
  
  
  free(X1);
  free(X2);
  free(t1);
  free(t2);
  free(n1);
  free(n2);
  free(l);
  free(nu);
  free(phi);
  free(T1);
  free(T2);
  free(N1);
  free(N2);
  free(kappa);
  
  return 0;
  
}


//--------------------関数--------------------//

// 領域の確保
double* make_vector(int n){
  
  double *a;
  
  // メモリの確保
  if( ( a = malloc( sizeof( double ) * n ) ) == NULL ){
    
    printf( "LACK of AVAILABLE MEMORY!" );
    exit(1);
    
  }
  
  return a;
  
}

// 領域の確保
double** make_matrix( int ROW, int COL ){
  
  int i;
  double **b;
  
  if( ( b = malloc( ROW * sizeof(double *) ) ) == NULL ){
    
    printf("メモリが確保できません\n");
    exit(1);
    
  }
  
  for( i = 0; i <= ROW; i++ ){
    
    b[i] = malloc( COL * sizeof(double) );
    
  }
  
  return b;
  
}

// 閉曲線のつなぎ
void connect( int N, double *x ){
  
  x[0] = x[N];
  x[N + 1] = x[1];
  
}

// つなぎ2個
void connect_double( int N, double *x, double *y ){
  
  connect(N,x);
  connect(N,y);
  
}

// 内積
double ip( double x1, double x2, double y1, double y2 ){
  
  return ( x1 * y1 + x2 * y2 );

}

// |y-x|
double DIST( double x1, double x2, double y1, double y2 ){
  
  return ( sqrt((y1 - x1) * (y1 - x1) + (y2 - x2) * (y2 - x2)) );
  
}

// 行列式
double DET( double x1, double x2, double y1, double y2 ){
  
  return ( x1 * y2 - y1 * x2 );
  
}

// 初期条件
void initial_condition( int N, double *x1, double *x2 ){
  
  int i;
  double u;
  double a1,a2,a3;

  /*
  for( i = 1; i <= N; i++ ){
    
    u = i * 1.0 / N;
    
    a1 = 1.8 * cos(2.0 * M_PI * u);
    a2 = 0.2 + sin(M_PI * u) * sin(6.0 * M_PI * u) * sin(2.0 * a1);
    a3 = 0.5 * sin(2.0 * M_PI * u) + sin(a1) + a2 * sin(2.0 * M_PI * u);
    
    x1[i] = 2.0 * 0.5 * a1;
    x2[i] = 2.0 * 0.54 * a3;
  }
  connect_double(N,x1,x2);
  */

  //----------初期条件(単位円)----------//
  
  for( i = 1; i <= N; i++ ){
    
    u = i * 2 * M_PI / N;
    
    x1[i] = cos(u);
    x2[i] = sin(u);

  }
  connect_double(N,x1,x2);
  
}

// x --> t,n,T,N,phi,kappa
void quantities( double t, int N, double *x1, double *x2, double *l, double *t1, double *t2, double *n1, double *n2, double *T1, double *T2, double *N1, double *N2, double *nu, double *phi, double *kappa ){
  
  int i;
  double D,I;
  
  for(i = 1; i <= N; i++){
    
    l[i] = DIST(x1[i-1],x2[i-1],x1[i],x2[i]);
    
    t1[i] = (x1[i] - x1[i-1]) / l[i];
    t2[i] = (x2[i] - x2[i-1]) / l[i];
    
    n1[i] = t2[i];
    n2[i] = -t1[i];
    
  }
  connect(N,l);
  connect_double(N,t1,t2);
  connect_double(N,n1,n2);
  
  /*
  if( t2[1] >= 0 ){
  
    nu[1] = acos(t1[1]);
    
  }

  else{

    nu[1] = -acos(t1[1]);
    
  }
  */
  
  for( i = 1; i <= N; i++ ){
    
    D = DET(t1[i],t2[i],t1[i+1],t2[i+1]);
    I = ip(t1[i],t2[i],t1[i+1],t2[i+1]);
    
    RANGE_CHECK(I,-1.0,1.0);
    
    //nu[i + 1] = nu[i] + D * acos(I);

    /*
    if( D >= 0.0 ){
      
      phi[i] = acos(I);
      
    }

    else{
      
      phi[i] = -acos(I);

    }
    
    
    }
    connect(N,phi);
    

    for( i = 1; i <= N; i++ ){
    
    nu[i + 1] = phi[i] + nu[i];
    
    }
    nu[0] = nu[1] - ( nu[N + 1] - nu[N] );
    */

    if( t2[1] < 0 ){
    
    nu[1] = -acos(t1[1]);
    
  }
  
  else{
    
    nu[1] = acos(t1[1]);
    
  }
  
  for( i = 1; i <= N; i++ ){
    
    nu[i + 1] = nu[i] + D * acos(I);
    
  }
  
  nu[0] = nu[1] - ( nu[N + 1] - nu[N] );

  }

  for( i = 1; i <= N; i++ ){
    
    phi[i] = nu[i + 1] - nu[i];
    
  }
  connect(N,phi);
  
  
  for( i = 1; i <= N; i++ ){
    
    T1[i] = ( t1[i] + t1[i + 1] ) / ( 2.0 * cos(phi[i] / 2.0) );
    T2[i] = ( t2[i] + t2[i + 1] ) / ( 2.0 * cos(phi[i] / 2.0) );
    
    N1[i] = T2[i];
    N2[i] = -T1[i];
    
    kappa[i] = ( tan(phi[i]) + tan(phi[i - 1]) ) / l[i];
    
  }
  connect_double(N,T1,T2);
  connect_double(N,N1,N2);
  connect(N,kappa);
  
}

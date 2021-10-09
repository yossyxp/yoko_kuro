
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

//--------------------定数--------------------//

#define NUM ( 240 )
#define TMAX ( 600.0 )
#define J ( 30 )
//#define dt ( 0.1 / ( NUM * NUM ) ) // 0.000006
//#define dt ( 1.0e-4 ) // 0.000006
#define dt ( 8.0e-5 ) // 0.000006
#define RANGE_CHECK(x, xmin, xmax) ( x = ( x < xmin ? xmin : ( x < xmax ?  x : xmax)));

//----------定数(結晶)----------//

#define T ( 281.15 ) // 絶対温度[K]
#define p_e ( 1.66e+3 ) // 平衡蒸気圧[dyn/cm^2]
#define v_c ( 3.25e-23 ) // 結晶相での水分子の体積[cm^3]
#define f_0 ( 8.3e-16 ) // 界面において1分子あたりが占める表面積[cm^2]
#define m ( 3.0e-23 ) // 水分子の質量[g]
#define e ( 4.5e-8 ) // ステップの高さ[m]
#define x_s ( 400 * e ) // 吸着分子がステップ滞在中に表面拡散する平均距離[m]
#define k_B ( 1.38e-16 ) // ボルツマン定数[erg/K]
#define energy ( 2.0e-6 ) // ステップの単位長さあたりの自由エネルギー[erg/cm]
#define E ( 40.0 ) // 拡散係数[m^2/s]
#define alpha_1 ( 0.1 ) // 凝縮定数

#define gamma ( v_c * p_e * E / ( k_B * T ) )
#define beta_max ( alpha_1 * v_c * p_e / sqrt(2 * M_PI * m * k_B * T) ) // 0.001995
#define sigma_star ( 9.5 * f_0 * energy / ( k_B * T * x_s ) )

//#define sigma_infty ( 17 ) // 初期値
#define sigma_infty ( 0.035 ) // 初期値
#define delta_beta ( 0.12 )

//----------定数(eps)----------//

#define epsl ( 1.0e-13 ) // 場合分け回避に用いるeps
#define eps_sim ( 1.0e-40 ) // 連立方程式に用いるeps


//--------------------関数--------------------//

double* make_vector( int N ); // 領域の確保(ベクトル)
double** make_matrix( int ROW, int COL ); // 領域の確保(行列)

void connect( double *x ); // 点のつなぎ(1個)
void connect_double( double* x, double *y ); // 点のつなぎ(2個)

double ip( double x1, double x2, double y1, double y2 ); // 内積
double DIST( double x1, double x2, double y1, double y2 ); // |y-x|
double DET( double x1, double x2, double y1, double y2 ); // 行列式
double EE( double x1, double x2, double y1, double y2 ); // 基本解
double E1( double x1, double x2, double y1, double y2 );
double E2( double x1, double x2, double y1, double y2 );

void gauss(int n, double **A, double *b, double *x); // Ax = b(0〜n-1までの行列)
void gauss2( int N, double **U, double *q, double *u );

void runge_kutta( double t, double *X1, double *X2 ); // ルンゲクッタ
void euler( double t, double *X1, double *X2 );
void F( double t, double *x1, double *x2, double *F1, double *F2 ); // 右辺
void pre( double t, double *x1, double *x2, double *T1, double *T2, double *N1, double *N2, double *V, double *W ); // 準備

void initial_condition( double *x1, double *x2 ); // 初期曲線

void quantities( double t, double *x1, double *x2, double *r, double *t1, double *t2, double *n1, double *n2, double *T1, double *T2, double *N1, double *N2, double *nu, double *phi, double *kappa ); //x --> t,n,T,N,phi,kappa
void measure( double t, double *x1, double *x2, double *L, double *A ); // x,l --> L,A

double omega( int n );

void supersaturation( double t, double *x1, double *x2, double *t1, double *t2, double *n1, double *n2, double *l, double *nu, double A, double *beta, double *u );
void velocity( double t, double *x1, double *x2, double *t1, double *t2, double *n1, double *n2, double *l, double *nu, double *phi, double *kappa, double L, double A, double *V, double *W ); // V,W
void normal_speed( double t, double *phi, double *beta, double *u, double *v, double *V ); // v,V
void tangent_speed( double t, double *l, double *phi, double *kappa, double *v, double *V, double L, double *W ); // W

double ff( double t, int i, int j, double *l, double *x1, double *x2, double *t1, double *t2, double *n1, double *n2, double beta, double a );
double gg( double t, int i, int j, double *l, double *x1, double *x2, double *t1, double *t2, double *n1, double *n2, double a );
double hh( double t, int i, int j, double *l, double *x1, double *x2, double *t1, double *t2, double a, double *beta );
double ii( double t, int j, double *l, double *x1, double *x2, double a, double R, double r_c );

//--------------------main-------------------//

int main(void){
  
  int i,z;
  double t;
  double *X1,*X2;
  X1 = make_vector(NUM + 2);
  X2 = make_vector(NUM + 2);
  
  double L,A;

  char file[5000];
  FILE *fp,*fp2;
  
  //fp = fopen("mfs.dat", "w");
  fp2 = fopen("L_A.dat", "w");
  
  t = 0.0;
  z = 0;

  initial_condition(X1,X2);

  sprintf(file, "./data/yoko_kuro%06d.dat", z);
  fp = fopen(file, "w");
  
  for( i = 0; i <= NUM; i++ ){
    
    fprintf(fp, "%f %f %f\n", X1[i], X2[i], t);
    
  }
  fprintf(fp,"\n");
  
  measure(t,X1,X2,&L,&A);
  
  fprintf(fp2, "t = %f, L = %f, A = %f\n", t, L, A );

  fclose(fp2);

  fclose(fp);
  
  while( t < TMAX ){
    
    runge_kutta(t,X1,X2);
    //euler(t,X1,X2);
    
    z++;
    t = z * dt;

    measure(t,X1,X2,&L,&A);
    printf("t = %f A = %.15f\n", t, A);
    
    if( z % 10000 == 0 ){
      
      sprintf(file, "./data/yoko_kuro%06d.dat", z / 10000 );
      fp = fopen(file, "w");

      measure(t,X1,X2,&L,&A);
      
      for( i = 0; i <= NUM; i++ ){
	
	fprintf(fp, "%f %f %f\n", X1[i], X2[i], t);

      }
      fprintf(fp,"\n");

      fclose(fp);
      fp2 = fopen("L_A.dat", "a");

      fprintf(fp2, "t = %f, L = %f, A = %f\n", t, L, A );
      fclose(fp2);

    }

  }

  free(X1);
  free(X2);
  
  return 0;
  
}

//--------------------関数--------------------//

// 領域の確保(ベクトル)
double* make_vector( int N ){
  
  double *a;
  
  // メモリの確保
  if( ( a = malloc( sizeof( double ) * N ) ) == NULL )
    {
      printf( "LACK of AVAILABLE MEMORY!" );
      exit(1);
    }
  
  return a;

}

// 領域の確保(行列)
double** make_matrix(int ROW, int COL){
  
  int i;
  double **b;
  
  // メモリの確保
  if( ( b = malloc( sizeof( double* ) * ROW ) ) == NULL ){
    
    printf( "LACK of AVAILABLE MEMORY!" );
    exit(1);

  }
  
  for( i = 0; i < ROW; i++){
    
    if( ( b[i] = malloc( sizeof( double ) * COL ) ) == NULL ){
      
      printf( "LACK of AVAILABLE MEMORY!" );
      free(b);
      exit(1);

    }

  }
  
  return b;

}

// 点のつなぎ(1個)
void connect( double *x ){
  
  x[0] = x[NUM];
  x[NUM + 1] = x[1];

}

// 点のつなぎ(2個)
void connect_double( double *x, double *y ){
  
  connect(x);
  connect(y);

}

// 内積
double ip( double x1, double x2, double y1, double y2 ){
  
  return ( x1 * y1 + x2 * y2 );

}

// |y-x|
double DIST(double x1, double x2, double y1, double y2){
  
  return ( sqrt( (y1 - x1) * (y1 - x1) + (y2 - x2) * (y2 - x2) ) );

}

// 行列式
double DET( double x1, double x2, double y1, double y2 ){
  
  return ( x1 * y2 - y1 * x2 );
  
}

// 基本解
double EE( double x1, double x2, double y1, double y2 ){
  
  return ( -log(DIST(y1,y2,x1,x2)) / ( 2.0 * M_PI ) );

}

// 基本解の勾配
double E1( double x1, double x2, double y1, double y2 ){
  
  double ry = DIST(y1,y2,x1,x2);
  
  return ( -( x1 - y1 ) / ( 2.0 * M_PI * ry * ry ) );
  
}

// 基本解の勾配
double E2( double x1, double x2, double y1, double y2 ){
  
  double ry = DIST(y1,y2,x1,x2);
  
  return ( -( x2 - y2 ) / ( 2.0 * M_PI * ry * ry ) );
  
}

// Ax = b(0〜n-1までの行列)
void gauss( int n, double **A, double *b, double *x ){
  
  int i,j,k,l,row_max;
  double max,temp,c;
  double temp_vec;
  double *temp_mat;
  
  for( i = 0; i < n - 1; i++ ){
    
    row_max = i;
    max = A[i][i];
    
    for( l = i + 1; l < n; l++ ){
      
      if( max < A[l][i] ){
	
	row_max = l;
	max = A[l][i];

      }		    

    }
    
    if( row_max > i ){
      
      temp_mat = A[i];
      A[i] = A[row_max];
      A[row_max] = temp_mat;
      
      temp_vec = b[i];
      b[i] = b[row_max];
      b[row_max] = temp_vec;

    }
      
    for( k = i + 1; k < n; k++ ){
      
      c = A[k][i] / A[i][i];

      for(j = i; j < n; j++){
	
	A[k][j] -= A[i][j] * c;

      }

      b[k] -= b[i] * c;

    }

  }
  
  for( i = n - 1; i >= 0; i--){
    
    temp = b[i];
    
    for( j = n - 1; j > i; j-- ){
      
      temp -= A[i][j] * x[j];

    }

    x[i] = temp / A[i][i];
    
  }

}

void gauss2( int n, double **U, double *q, double *u ){


  int i,j,k,ip;
  double eta, tmp, amax; 
  
  //----------ガウスの消去法----------//
  
  //----------消去----------//
    
    for( k = 1; k <= n - 1; k++ ){
      
      amax = fabs(U[k][k]);
      ip = k;
      
      for( i = k + 1; i <= n; i++ ){
	
	if( fabs(U[i][k]) > amax ){
	  
	  amax = fabs(U[i][k]);
	  ip = i;
	  
	}
	
      }
      
      if( amax < eps_sim ){
	
	printf("入力した行列は正則ではない\n");
	
      }
      
      if( ip != k ){
	
	for( j = k; j <= n; j++ ){
	  
	  tmp = U[k][j];
	  U[k][j] = U[ip][j];
	  U[ip][j] = tmp;
	  
	}
	
	tmp = q[k];
	q[k] = q[ip];
	q[ip] = tmp;
	
      }
      
      for( i = k + 1; i <= n; i++ ){
	
	eta = - U[i][k] / U[k][k];
	
	for( j = k + 1; j <= n; j++ ){
	  
	  U[i][j] = U[i][j] + eta * U[k][j];
	  
	}
	
	q[i] = q[i] + eta * q[k];
	
      }
      
    }
    
    u[n] = q[n] / U[n][n];
    
    for( k = n - 1; k >= 1; k-- ){
      
      tmp = q[k];
      
      for( j = k + 1; j <= n; j++ ){
	
	tmp = tmp - U[k][j] * u[j];
	
      }
      
      u[k] = tmp / U[k][k];
      
    }

}

// ルンゲクッタ
void runge_kutta( double t, double *X1, double *X2 ){
  
  int i;
  double t_temp;
  double *x_temp1,*x_temp2,*F1,*F2;
  double *k11,*k12,*k21,*k22,*k31,*k32,*k41,*k42;
  
  k11 = make_vector(NUM + 2);
  k12 = make_vector(NUM + 2);
  k21 = make_vector(NUM + 2);
  k22 = make_vector(NUM + 2);
  k31 = make_vector(NUM + 2);
  k32 = make_vector(NUM + 2);
  k41 = make_vector(NUM + 2);
  k42 = make_vector(NUM + 2);

  x_temp1 = make_vector(NUM + 2);
  x_temp2 = make_vector(NUM + 2);
  F1 = make_vector(NUM + 2);
  F2 = make_vector(NUM + 2);

  F(t,X1,X2,F1,F2);
  
  for( i = 1; i <= NUM; i++ ){
    
      k11[i] = F1[i];
      k12[i] = F2[i];
      x_temp1[i] = X1[i] + dt * k11[i] / 2.0;
      x_temp2[i] = X2[i] + dt * k12[i] / 2.0;

  }
  connect_double(x_temp1,x_temp2);

  t_temp = t + dt/2.0;

  F(t_temp,x_temp1,x_temp2,F1,F2);      

  for( i = 1; i <= NUM; i++ ){
    
      k21[i] = F1[i];
      k22[i] = F2[i];
      x_temp1[i] = X1[i] + dt * k21[i] / 2.0;
      x_temp2[i] = X2[i] + dt * k22[i] / 2.0;

  }
  connect_double(x_temp1,x_temp2);

  F(t_temp,x_temp1,x_temp2,F1,F2);
  
  for( i = 1; i <= NUM; i++ ){
    
    k31[i] = F1[i];
    k32[i] = F2[i];
    x_temp1[i] = X1[i] + k31[i] * dt;
    x_temp2[i] = X2[i] + k32[i] * dt;

  }
  connect_double(x_temp1,x_temp2);
  
  t_temp = t + dt;
      
  F(t_temp,x_temp1,x_temp2,F1,F2);
  
  for( i = 1; i <= NUM; i++ ){
    
    k41[i] = F1[i];
    k42[i] = F2[i];
    
    X1[i] = X1[i] + dt * ( k11[i] + 2.0 * k21[i] + 2.0 * k31[i] + k41[i] ) / 6.0;
    X2[i] = X2[i] + dt * ( k12[i] + 2.0 * k22[i] + 2.0 * k32[i] + k42[i] ) / 6.0;

  }
  connect_double(X1,X2);

  free(k11);
  free(k12);
  free(k21);
  free(k22);
  free(k31);
  free(k32);
  free(k41);
  free(k42);
  free(x_temp1);
  free(x_temp2);
  free(F1);
  free(F2);
  
}

void euler( double t, double *X1, double *X2 ){

  int i;
  double *F1,*F2;

  F1 = make_vector(NUM + 2);
  F2 = make_vector(NUM + 2);

  F(t,X1,X2,F1,F2);

  for( i = 1; i <= NUM; i++ ){
    
    X1[i] = X1[i] + dt * F1[i];
    X2[i] = X2[i] + dt * F2[i];

  }
  connect_double(X1,X2);
  
  free(F1);
  free(F2);
  
}

// 右辺
void F( double t, double *x1, double *x2, double *F1, double *F2 ){
  
  int i;
  
  double *T1,*T2,*N1,*N2;
  T1 = make_vector(NUM + 2);
  T2 = make_vector(NUM + 2);
  N1 = make_vector(NUM + 2);
  N2 = make_vector(NUM + 2);

  double *V,*W;
  V = make_vector(NUM + 2);
  W = make_vector(NUM + 2);
  
  pre(t,x1,x2,T1,T2,N1,N2,V,W);
  
  for( i = 1; i <= NUM; i++ ){
    
    F1[i] = V[i]*N1[i] + W[i]*T1[i];
    F2[i] = V[i]*N2[i] + W[i]*T2[i];
    //printf("V = %f\n", V[i]);

  }
  connect_double(F1,F2);

  free(W);
  free(V);
  free(T1);
  free(T2);
  free(N1);
  free(N2);
  
}

//x --> T,N,V,W
void pre( double t, double *x1, double *x2, double *T1, double *T2, double *N1, double *N2, double *V, double *W ){
  
  double *l;
  double *t1,*t2,*n1,*n2;
  double *nu;
  double *phi;
  double *kappa;
  double L,A;

  l = make_vector(NUM + 2);
  nu = make_vector(NUM + 2);
  phi = make_vector(NUM + 2);
  kappa = make_vector(NUM + 2);

  t1 = make_vector(NUM + 2);
  t2 = make_vector(NUM + 2);
  n1 = make_vector(NUM + 2);
  n2 = make_vector(NUM + 2);
  
  // T,N
  quantities(t,x1,x2,l,t1,t2,n1,n2,T1,T2,N1,N2,nu,phi,kappa);
  
  // L,A
  measure(t,x1,x2,&L,&A);
  
  // V,W
  velocity(t,x1,x2,t1,t2,n1,n2,l,nu,phi,kappa,L,A,V,W);

  free(t1);
  free(t2);
  free(n1);
  free(n2);
  
  free(kappa);
  free(phi);
  free(nu);
  free(l);
  
}

// 初期曲線
void initial_condition( double *x1, double *x2 ){

  int i, k;
  double u,lambda;

  /*
  for (i = 1; i <= NUM; i++)
  {
    u = i * 2 * M_PI / NUM;
  }

  for (k = 0; k < 6; k++)
  {
    for (i = k * (NUM / 6) + 1; i <= (k + 1) * (NUM / 6); i++)
    {
      lambda = (i - k * (NUM / 6)) * 1.0 / (NUM / 6);
      if (k == 0)
      {
        x1[i] = (1.0 - lambda) * 1.0 + lambda * 1.0 / 2.0;
        x2[i] = (1.0 - lambda) * 0.0 + lambda * sqrt(3.0) / 2.0;
      }
      else if (k == 1)
      {
        x1[i] = (1.0 - lambda) * 1.0 / 2.0 + lambda * (-1.0 / 2.0);
        x2[i] = sqrt(3.0) / 2.0;
      }
      else if (k == 2)
      {
        x1[i] = (1.0 - lambda) * (-1.0 / 2.0) + lambda * (-1.0);
        x2[i] = (1.0 - lambda) * (sqrt(3.0) / 2.0) + lambda * 0.0;
      }
      else if (k == 3)
      {
        x1[i] = (1.0 - lambda) * (-1.0) + lambda * (-1.0 / 2.0);
        x2[i] = (1.0 - lambda) * 0.0 + lambda * (-sqrt(3.0) / 2.0);
      }
      else if (k == 4)
      {
        x1[i] = (1.0 - lambda) * (-1.0 / 2.0) + lambda * (1.0 / 2.0);
        x2[i] = -sqrt(3.0) / 2.0;
      }
      else
      {
        x1[i] = (1.0 - lambda) * (1.0 / 2.0) + lambda * (1.0);
        x2[i] = (1.0 - lambda) * (-sqrt(3.0) / 2.0) + lambda * 0.0;
      }

      x1[i] = 1.0e-2 * x1[i];
      x2[i] = 1.0e-2 * x2[i];
    }
  }
  connect_double(x1, x2);
  */

  
  for( i = 1; i <= NUM; i++ ){
    
    u = i * 2 * M_PI / NUM;
    
    x1[i] = 5.0e-3 * cos(u);
    x2[i] = 5.0e-3 * sin(u);
    
  }
  connect_double(x1,x2);
  
  
}

//x --> t,n,T,N,phi,kappa
void quantities( double t, double *x1, double *x2, double *l, double *t1, double *t2, double *n1, double *n2, double *T1, double *T2, double *N1, double *N2, double *nu, double *phi, double *kappa ){
  
  int i;
  double D,I,D_sgn;
  
  for( i = 1; i <= NUM; i++ ){
    
    l[i] = DIST(x1[i - 1],x2[i - 1],x1[i],x2[i]);
    
    t1[i] = ( x1[i] - x1[i - 1] ) / l[i];
    t2[i] = ( x2[i] - x2[i - 1] ) / l[i];
    
    n1[i] = -t2[i];
    n2[i] = t1[i];

  }
  connect(l);
  connect_double(t1,t2);
  connect_double(n1,n2);
  

  RANGE_CHECK(t1[1],-1.0,1.0);
  
  if( t2[1] >= 0 ){
    
    nu[1] = acos(t1[1]);
    
  }
  
  else{
    
    nu[1] = -acos(t1[1]);
    
  }
  
  for( i = 1; i <= NUM; i++){
    
    D = DET(t1[i],t2[i],t1[i + 1],t2[i + 1]);
    I = ip(t1[i],t2[i],t1[i + 1],t2[i + 1]);
    
    RANGE_CHECK(I,-1.0,1.0);
    
    if( D < 0 ){
      
      D_sgn = -1;
      
    }
    
    else if( D > 0 ){
      
      D_sgn = 1;
      
    }
    
    else{
      
      D_sgn = 0;
      
    }
    
    nu[i + 1] = nu[i] + D_sgn * acos(I);
    
  }
  nu[0] = nu[1] - ( nu[NUM + 1] - nu[NUM] );
  
  
  for( i = 1; i <= NUM; i++ ){
    
    phi[i] = nu[i + 1] - nu[i];
    
  }
  connect(phi);
  

  for( i = 1; i <= NUM; i++){
    
    T1[i] = ( t1[i] + t1[i + 1] ) / ( 2.0 * cos(phi[i] / 2.0) );
    T2[i] = ( t2[i] + t2[i + 1] ) / ( 2.0 * cos(phi[i] / 2.0) );

    N1[i] = -T2[i];
    N2[i] = T1[i];
    
    kappa[i] = ( tan(phi[i] / 2.0) + tan(phi[i - 1] / 2.0) ) / l[i];
    
  }
  connect_double(T1,T2);
  connect_double(N1,N2);
  connect(kappa);

  //printf("%f %f\n", N1[0], N2[0]);

}

// x,l --> L,A
void measure( double t, double *x1, double *x2, double *L, double *A ){
  
  int i;
  
  *L = 0.0;
  *A = 0.0;
  
  for( i = 1; i <= NUM; i++){
    
    *L += DIST(x1[i],x2[i],x1[i - 1],x2[i - 1]);
    *A += DET(x1[i - 1],x2[i - 1],x1[i],x2[i]);

  }
  *A = *A / 2.0;

}

// x,nu --> beta,u
void supersaturation( double t, double *x1, double *x2, double *t1, double *t2, double *n1, double *n2, double *l, double *nu, double A, double *beta, double *u ){
  
  int i,j,k;
  double r_c,R;
  double dx_pi,dx_sim;
  double z1;
  double **U;
  double *q;

  U = make_matrix(NUM + 2,NUM + 2);
  q = make_vector(NUM + 2);

  
  r_c = sqrt(A / M_PI);
  
  R = 6.5 * r_c;
  
  /*
  if( t == 0.0 ){
    
    for( i = 1; i <= NUM; i++ ){
      
      //beta[i] = beta_max * cos(6 * nu[i]) + beta_max + beta_max;
      //beta[i] = beta_max * 2 * x_s * tan(3 * nu[i]) * tanh(e / ( 2 * x_s * tan(3 * nu[i]) )) / e;
      beta[i] = beta_max * ( sin(( atan(tan(3*nu[i]-M_PI/2.0)) + M_PI/2.0 ) / 3.0) + sin(M_PI/3.0 - ( atan(tan(3*nu[i]-M_PI/2.0)) + M_PI/2.0 ) / 3.0) ) / sin(M_PI/3.0);
      //beta[i] = beta_max * tan(3 * nu[i]) * tanh(1.0 / tan(3 * nu[i])) + beta_max;
      //beta[i] = ( 1 + 0.1 * ( fabs(cos(3 * nu[i] + M_PI / 2.0)) - 0.5 ) ) * beta_max;
      //beta[i] = beta_max * ( 1 + ( 0.99 / 35.0 ) * cos(6 * nu[i]) );
      //beta[i] = beta_max * cos(4 * ( nu[i] - M_PI / 4.0 )) + beta_max + beta_max;
      //beta[i] = beta_max;

      if( nu[i] > -1.0e-5 ){
	
	if( nu[i] < 1.0e-5 ){
	  
	  printf("a");
	  
	  beta[i] = ( 1 - delta_beta ) * beta_max; 
	  
	}
	
      }

      if( nu[i] > M_PI / 3.0 - 1.0e-5 ){
	
	if( nu[i] < M_PI / 3.0 + 1.0e-5 ){
	  
	  printf("b");
	  
	  beta[i] = ( 1 - delta_beta ) * beta_max; 
	  
	}
	
      }

      if( nu[i] > 2 * M_PI / 3.0 - 1.0e-5 ){
	
	if( nu[i] < 2 * M_PI / 3.0 + 1.0e-5 ){
	  
	  printf("c");
	  
	  beta[i] = ( 1 - delta_beta ) * beta_max; 
	  
	}
	
      }
      
      if( nu[i] > M_PI - 1.0e-5 ){
	
	if( nu[i] < M_PI + 1.0e-5 ){
	  
	  printf("d");
	  
	  beta[i] = ( 1 - delta_beta ) * beta_max; 
	  
	}
	
      }

      if( nu[i] > 4 * M_PI / 3.0 - 1.0e-5 ){
	
	if( nu[i] < 4 * M_PI / 3.0 + 1.0e-5 ){
	  
	  printf("e");
	  
	  beta[i] = ( 1 - delta_beta ) * beta_max; 
	  
	}
	
      }
      
      if( nu[i] > 5 * M_PI / 3.0 - 1.0e-5 ){
	
	if( nu[i] < 5 * M_PI / 3.0 + 1.0e-5 ){
	  
	  printf("f");
	  
	  beta[i] = ( 1 - delta_beta ) * beta_max; 
	  
	}
	
      }
      
      if( nu[i] > 2 * M_PI - 1.0e-5 ){
	
	if( nu[i] < 2 * M_PI + 1.0e-5 ){
	  
	  printf("g");
	  
	  beta[i] = ( 1 - delta_beta ) * beta_max; 
	  
	}
	
      }
      
      if( nu[i] > 7 * M_PI / 3.0 - 1.0e-5 ){
	
	if( nu[i] < 7 * M_PI / 3.0 + 1.0e-5 ){
	  
	  printf("h");
	  
	  beta[i] = ( 1 - delta_beta ) * beta_max; 
	  
	}
	
      }

      if( nu[i] > 8 * M_PI / 3.0 - 1.0e-5 ){
	
	if( nu[i] < 8 * M_PI / 3.0 + 1.0e-5 ){
	  
	  printf("i");
	  
	  beta[i] = ( 1 - delta_beta ) * beta_max; 
	  
	}
	
      }					  
      
    }
    connect(beta);
    
  }

  else{
    
    for( i = 1; i <= NUM; i++ ){
      
      //beta[i] = beta_max * cos(6 * nu[i]) + beta_max + beta_max;
      //beta[i] = beta_max * 2 * x_s * tan(3 * nu[i]) * tanh(e / ( 2 * x_s * tan(3 * nu[i]) )) / e;
      //beta[i] = beta_max * tan(3 * nu[i]) * tanh(1.0 / tan(3 * nu[i]));
      beta[i] = beta_max * ( sin(( atan(tan(3*nu[i]-M_PI/2.0)) + M_PI/2.0 ) / 3.0) + sin(M_PI/3.0 - ( atan(tan(3*nu[i]-M_PI/2.0)) + M_PI/2.0 ) / 3.0) ) / sin(M_PI/3.0);
      //beta[i] = ( 1 + 0.1 * ( fabs(cos(3 * nu[i] + M_PI / 2.0)) - 0.5 ) ) * beta_max;
      //beta[i] = beta_max * ( 1 + ( 0.99 / 35.0 ) * cos(6 * nu[i]) );
      //beta[i] = beta_max * cos(4 * ( nu[i] - M_PI / 4.0 )) + beta_max + beta_max;
      //beta[i] = beta_max;
      
      if( nu[i] > -1.0e-5 ){
	
	if( nu[i] < 1.0e-5 ){
	  
	  printf("a");
	  
	  beta[i] = beta_max * u[i] * tanh(sigma_star / u[i]) / sigma_star;
	  
	}
	
      }

      if( nu[i] > M_PI / 3.0 - 1.0e-5 ){
	
	if( nu[i] < M_PI / 3.0 + 1.0e-5 ){
	  
	  printf("b");
	  
	  beta[i] = beta_max * u[i] * tanh(sigma_star / u[i]) / sigma_star;
	  
	}
	
      }

      if( nu[i] > 2 * M_PI / 3.0 - 1.0e-5 ){
	
	if( nu[i] < 2 * M_PI / 3.0 + 1.0e-5 ){
	  
	  printf("c");
	  
	  beta[i] = beta_max * u[i] * tanh(sigma_star / u[i]) / sigma_star;
	  
	}
	
      }
      
      if( nu[i] > M_PI - 1.0e-5 ){
	
	if( nu[i] < M_PI + 1.0e-5 ){
	  
	  printf("d");
	  
	  beta[i] = beta_max * u[i] * tanh(sigma_star / u[i]) / sigma_star;
	  
	}
	
      }

      if( nu[i] > 4 * M_PI / 3.0 - 1.0e-5 ){
	
	if( nu[i] < 4 * M_PI / 3.0 + 1.0e-5 ){
	  
	  printf("e");
	  
	  beta[i] = beta_max * u[i] * tanh(sigma_star / u[i]) / sigma_star;
	  
	}
	
      }
      
      if( nu[i] > 5 * M_PI / 3.0 - 1.0e-5 ){
	
	if( nu[i] < 5 * M_PI / 3.0 + 1.0e-5 ){
	  
	  printf("f");
	  
	  beta[i] = beta_max * u[i] * tanh(sigma_star / u[i]) / sigma_star;
	  
	}
	
      }
      
      if( nu[i] > 2 * M_PI - 1.0e-5 ){
	
	if( nu[i] < 2 * M_PI + 1.0e-5 ){
	  
	  printf("g");
	  
	  beta[i] = beta_max * u[i] * tanh(sigma_star / u[i]) / sigma_star;
	  
	}
	
      }
      
      if( nu[i] > 7 * M_PI / 3.0 - 1.0e-5 ){
	
	if( nu[i] < 7 * M_PI / 3.0 + 1.0e-5 ){
	  
	  printf("h");
	  
	  beta[i] = beta_max * u[i] * tanh(sigma_star / u[i]) / sigma_star;
	  
	}
	
      }

      if( nu[i] > 8 * M_PI / 3.0 - 1.0e-5 ){
	
	if( nu[i] < 8 * M_PI / 3.0 + 1.0e-5 ){
	  
	  printf("i");
	  
	  beta[i] = beta_max * u[i] * tanh(sigma_star / u[i]) / sigma_star;
	  
	}
	
      }	
      
    }
    connect(beta);
    
  }
  */

  
  for( i = 1; i <= NUM; i++ ){
      
      //beta[i] = beta_max * cos(6 * nu[i]) + beta_max + beta_max;
      //beta[i] = beta_max * 2 * x_s * tan(3 * nu[i]) * tanh(e / ( 2 * x_s * tan(3 * nu[i]) )) / e;
      //beta[i] = beta_max * tan(3 * nu[i]) * tanh(1.0 / tan(3 * nu[i]));
      //beta[i] = beta_max * tan(3 * nu[i]) * tanh(0.2 / tan(3 * nu[i])) / 0.2;
      beta[i] = beta_max * ( sin(( atan(tan(3*nu[i]-M_PI/2.0)) + M_PI/2.0 ) / 3.0) + sin(M_PI/3.0 - ( atan(tan(3*nu[i]-M_PI/2.0)) + M_PI/2.0 ) / 3.0) ) / sin(M_PI/3.0);
      //beta[i] = ( 1 + 0.1 * ( fabs(cos(3 * nu[i] + M_PI / 2.0)) - 0.5 ) ) * beta_max;
      //beta[i] = beta_max * cos(4 * ( nu[i] - M_PI / 4.0 )) + beta_max + beta_max;
      //beta[i] = beta_max;

  }
  connect(beta);
  
  
  for( i = 1; i <= NUM; i++ ){
    
    dx_pi = 2 * M_PI / 100.0;
    
    q[i] = dx_pi * ( ii(t,i,l,x1,x2,0,R,r_c)
		     + ii(t,i,l,x1,x2,2*M_PI,R,r_c) ) / 2.0;
    
    for( k = 1; k < 100; k++ ){
      
      z1 = k * dx_pi;
      
      q[i] = q[i] + dx_pi * ii(t,i,l,x1,x2,z1,R,r_c);
      
    }

  }
  connect(q);

  //printf("%.15f\n", q[1]);
  
  for( i = 1; i <= NUM; i++ ){
    
    for( j = 1; j <= NUM; j++ ){
      
      if( j == i ){
	
	U[i][j] = 0.5 - ( beta[i] / ( 2 * M_PI * gamma ) ) * l[i] * ( 1 - log(l[i] / 2.0) );
	
      }

      else{
	
	//----------数値積分----------//
	
	dx_sim = l[j] / J;
	
	U[i][j] = dx_sim * ( ( ff(t,i,j,l,x1,x2,t1,t2,n1,n2,beta[j],0)
			     + ff(t,i,j,l,x1,x2,t1,t2,n1,n2,beta[j],l[j]) ) ) / 2.0;
	
	for( k = 1; k < J; k++ ){
	  
	  z1 = k * dx_sim;
	  
	  U[i][j] = U[i][j] + dx_sim * ff(t,i,j,l,x1,x2,t1,t2,n1,n2,beta[j],z1);
	  
	}

      /*
      else{
	
	//----------数値積分----------//
	
	dx_sim = l[j] / J;
	
	U[i][j] = dx_sim * (
			    ( gg(t,i,j,l,x1,x2,t1,t2,n1,n2,0)
			    + gg(t,i,j,l,x1,x2,t1,t2,n1,n2,l[j]) )
			    + ( hh(t,i,j,l,x1,x2,t1,t2,0,beta)
			    + hh(t,i,j,l,x1,x2,t1,t2,l[j],beta) )
	 ) / 2.0;
	
	for( k = 1; k < J; k++ ){
	  
	  z1 = k * dx_sim;
	  
	  U[i][j] = U[i][j] + dx_sim *
	    ( gg(t,i,j,l,x1,x2,t1,t2,n1,n2,z1)
	      + hh(t,i,j,l,x1,x2,t1,t2,z1,beta)
	      );
					 
	}
      */
      }
      
    }
    
  }
  
  /*
  for( j = 1; j <= NUM; j++ ){

    U[0][j] = U[NUM][j];
    
  }

  for( i = 0; i <= NUM; i++ ){

    U[i][0] = U[i][NUM];
    
  }
  

  gauss(NUM,U,q,u);
  u[NUM] = u[0];
  u[NUM + 1] = u[1];
  */

  gauss2(NUM,U,q,u);
  connect(u);

  
  for( i = 1; i <= NUM; i++ ){
    
    //printf("u[%d] = %.15f\n", i, u[i]);
    // 戻せ
  }
  

  for( i = 0; i <= NUM + 1; i++ ){
    
    free((void *)U[i]);
    
  }
  free((void *)U);
  
  free(q);
  
}

//緩和項
double omega( int n ){
  
  return ( 10.0 * n );

}

// x,n,r,t,ohi,kappa,L --> V,W
void velocity( double t, double *x1, double *x2, double *t1, double *t2, double *n1, double *n2, double *l, double *nu, double *phi, double *kappa, double L, double A, double *V, double *W ){
  
  int i;
  double *beta;
  double *u;
  double *v;
  
  beta = make_vector(NUM + 2);
  u = make_vector(NUM + 2);
  v = make_vector(NUM + 2);
  
  supersaturation(t,x1,x2,t1,t2,n1,n2,l,nu,A,beta,u);
  normal_speed(t,phi,beta,u,v,V);
  tangent_speed(t,l,phi,kappa,v,V,L,W);


  free(beta);
  free(u);
  free(v);
  
}

// v,V
void normal_speed( double t, double *phi, double *beta, double *u, double *v, double *V ){
  
  int i;
  
  for( i = 1; i <= NUM; i++ ){
    
    v[i] = beta[i] * u[i];

  }
  connect(v);

  for( i = 1; i <= NUM; i++ ){
    
    V[i] = ( v[i] + v[i + 1] ) / ( 2.0 * cos(phi[i] / 2.0) );
  }
  connect(V);
  
}

// r,phi,kappa,v,V,L --> W
void tangent_speed( double t, double *l, double *phi, double *kappa, double *v, double *V, double L, double *W ){
  
  int i;
  double *psi,*PSI;
  double L_dot;
  double a,b,c;
  
  psi = make_vector(NUM + 1);
  PSI = make_vector(NUM + 1);
  
  psi[1] = 0.0;
  L_dot = 0.0;
  
  for(i = 1; i <= NUM; i++ ){
    
    L_dot += kappa[i] * v[i] * l[i];

  }
  
  for( i = 2; i <= NUM; i++ ){
    
    psi[i] = ( L_dot / NUM ) - V[i] * sin(phi[i] / 2.0) - V[i-1] * sin(phi[i - 1] / 2.0) + ( ( L / NUM ) - l[i] ) * omega(NUM);

  }
  
  PSI[1] = psi[1];
  
  for( i = 2; i <= NUM; i++ ){
    
    PSI[i] = PSI[i-1] + psi[i];

  }

  a = 0.0;
  b = 0.0;

  for( i = 1; i <= NUM; i++ ){
    
    a += PSI[i] / cos(phi[i] / 2.0);
    b += 1.0 / cos(phi[i] / 2.0);

  }
  c = -a / b;

  for( i = 1; i <= NUM; i++ ){
    
    W[i] = ( PSI[i] + c ) / cos(phi[i] / 2.0);

  }
  connect(W);

  free(PSI);
  free(psi);

}

double gg( double t, int i, int j, double *l, double *x1, double *x2, double *t1, double *t2, double *n1, double *n2, double a ){

  double x1_star,x2_star;

  x1_star = ( x1[i - 1] + x1[i] ) / 2.0;
  x2_star = ( x2[i - 1] + x2[i] ) / 2.0;
  
  return (
	  /*
	  -( ( x1[j - 1] - ( x1[i - 1] + a * t1[i] ) ) * n1[j] + ( x2[j - 1] - ( x2[i - 1] + a * t2[i] ) ) * n2[j] ) / ( 2 * M_PI * ( ( x1[j - 1] - ( x1[i - 1] + a * t1[i] ) ) * ( x1[j - 1] - ( x1[i - 1] + a * t1[i] ) ) + ( x2[j - 1] - ( x2[i - 1] + a * t2[i] ) ) * ( x2[j - 1] - ( x2[i - 1] + a * t2[i] ) ) + epsl * epsl ) )
	  */
	  
	  /*
	  ( ( x1_star - ( x1[j] - a * t1[j] ) ) * n1[j]
	  + ( x2_star - ( x2[j] - a * t2[j] ) ) * n2[j]
	  )
	  / ( 2 * M_PI * (
			  ( x1_star - ( x1[j] - a * t1[j] ) ) * ( x1_star - ( x1[j] - a * t1[j] ) )
			+ ( x2_star - ( x2[j] - a * t2[j] ) ) * ( x2_star - ( x2[j] - a * t2[j] ) )
			  + epsl * epsl
			 )
	    )
	  */

	  
	  ( E1(x1_star,x2_star,( x1[j - 1] + a * t1[j] ),( x2[j - 1] + a * t2[j] )) * n1[j]
	    + E2(x1_star,x2_star,( x1[j - 1] + a * t1[j] ),( x2[j - 1] + a * t2[j] )) * n2[j] )
	  
	  );

}

double hh( double t, int i, int j, double *l, double *x1, double *x2, double *t1, double *t2, double a, double *beta ){

  double x1_star,x2_star;

  x1_star = ( x1[i - 1] + x1[i] ) / 2.0;
  x2_star = ( x2[i - 1] + x2[i] ) / 2.0;
  
  return (
	  /*
	  ( ( k_B * T * beta[i] ) / ( 2 * M_PI * E * p_e * v_c ) ) * log( ( x1[j - 1] - ( x1[i - 1] + a * t1[i] ) ) * ( x1[j - 1] - ( x1[i - 1] + a * t1[i] ) ) + ( x2[j - 1] - ( x2[i - 1] + a * t2[i] ) ) * ( x2[j - 1] - ( x2[i - 1] + a * t2[i] ) ) + epsl * epsl )
	  */

	  /*
	  ( ( k_B * T * beta[j] ) / ( 2 * M_PI * E * p_e * v_c ) )
	  * log( sqrt( ( x1_star - ( x1[j] - a * t1[j] ) ) * ( x1_star - ( x1[j] - a * t1[j] ) )
		     + ( x2_star - ( x2[j] - a * t2[j] ) ) * ( x2_star - ( x2[j] - a * t2[j] ) )
		     )
	         + epsl * epsl
	       )
	  */

	   ( beta[j] / gamma ) * EE(x1_star,x2_star,( x1[j - 1] + a * t1[j] ),( x2[j - 1] + a * t2[j] ))
	  
	  );

}

double ff( double t, int i, int j, double *l, double *x1, double *x2, double *t1, double *t2, double *n1, double *n2, double beta, double a ){

  double x1_star,x2_star;
  double y1,y2;
  
  x1_star = ( x1[i - 1] + x1[i] ) / 2.0;
  x2_star = ( x2[i - 1] + x2[i] ) / 2.0;
  y1 = x1[j] - a * t1[j];
  y2 = x2[j] - a * t2[j];

  return(

	 ( ( E1(x1_star,x2_star,y1,y2) * n1[j] )
	  + ( E2(x1_star,x2_star,y1,y2) * n2[j] ) )
	 -
	 ( ( beta / gamma ) * EE(x1_star,x2_star,y1,y2) )
	 
	 
	 );
  
}

double ii( double t, int i, double *l, double *x1, double *x2, double a, double R, double r_c ){

  double x1_star,x2_star;
  double y1,y2;

  x1_star = ( x1[i - 1] + x1[i] ) / 2.0;
  x2_star = ( x2[i - 1] + x2[i] ) / 2.0;
  y1 = R * cos(a);
  y2 = R * sin(a);

  return (
	  /*	  
	  ( -( log( sqrt( ( x1[j] - R * cos(a) ) * ( x1[j] - R * cos(a) ) + ( x2[j] - R * sin(a) ) * ( x2[j] - R * sin(a) ) ) + epsl * epsl ) / ( 2 * M_PI * R * log(R / r_c) ) )
	    + ( ( x1[j] * cos(a) + x2[j] * sin(a) - R ) / ( 2 * M_PI * ( ( x1[j] - R * cos(a) ) * ( x1[j] - R * cos(a) ) + ( x2[j] - R * sin(a) ) * ( x2[j] - R * sin(a) ) + epsl * epsl ) ) ) ) * R * sigma_infty
	  */

	  /*
	  ( ( log( sqrt( ( x1_star - R * cos(a) ) * ( x1_star - R * cos(a) ) + ( x2_star - R * sin(a) ) * ( x2_star - R * sin(a) ) ) + epsl * epsl ) / ( 2 * M_PI * R * log(R / r_c) ) )
	    + ( ( x1_star * cos(a) + x2_star * sin(a) - R ) / ( 2 * M_PI * ( ( x1_star - R * cos(a) ) * ( x1_star - R * cos(a) ) + ( x2_star - R * sin(a) ) * ( x2_star - R * sin(a) ) + epsl * epsl ) ) ) ) * R * sigma_infty
	  */
	  
	  
	  (
	   ( EE(x1_star,x2_star,y1,y2) / log(R / r_c) )	   
	   -
	    ( ( ( E1(x1_star,x2_star,y1,y2) * cos(a) )
	      + ( E2(x1_star,x2_star,y1,y2) * sin(a) ) ) * R )
	   ) * sigma_infty
	  
	  
	  );

}

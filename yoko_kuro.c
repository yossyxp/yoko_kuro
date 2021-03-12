
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

//--------------------定数--------------------//

#define Z ( 2048 )
#define TMAX ( 10.1 )
#define RANGE_CHECK( x, xmin, xmax ) ( x = ( x < xmin ? xmin : ( x < xmax ?  x : xmax ) ) );


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
#define E ( 40 ) // 拡散係数[m^2/s]
#define alpha_1 ( 0.1 ) // 凝縮定数

//#define beta_max ( alpha_1 * v_c * p_e / sqrt(2 * M_PI * m * k_B * T) )
#define beta_max ( 0.01 )
#define sigma_star ( 9.5 * f_0 * energy / ( k_B * T * x_s ) )

//#define sigma_infty ( 17 ) // 初期値
#define sigma_infty ( 0.035 ) // 初期値
#define delta_beta ( 0.4 )


//----------定数(eps)----------//

#define epsl ( 1.0e-15 ) // 場合分け回避に用いるeps
#define eps_sim ( 1.0e-40 ) // 連立方程式に用いるeps


//--------------------関数--------------------//

double* make_vector( int n ); // 領域の確保(ベクトル)
double** make_matrix( int ROW, int COL ); // 領域の確保(行列)
void connect( int N, double *x ); // 閉曲線のつなぎ
void connect_double( int N, double* x, double *y ); // つなぎ2個

double ip( double x1, double x2, double y1, double y2 ); // 内積
double DIST( double x1, double x2, double y1, double y2 ); // |y - x|
double DET( double x1, double x2, double y1, double y2 ); // 行列式

void runge_qutta( double t, double dt, int N, double *X1, double *X2 ); // ルンゲクッタ
void F( double t, int N, double *x1, double *x2, double *F1, double *F2 ); // 右辺
void ODE_pre( double t, int N, double *x1, double *x2, double *T1, double *T2, double *N1, double *N2, double *V, double *W ); // x --> T,N,V,W
void initial_condition( int N, double *x1, double *x2 ); // 初期条件
void quantities( double t, int N, double *x1, double *x2, double *l, double *t1, double *t2, double *n1, double *n2, double *T1, double *T2, double *N1, double *N2, double *nu, double *phi, double *kappa ); //x --> t,n,T,N,phi,kappa
void measure( double t, int N, double *x1, double *x2, double *L, double *A ); // x,l --> L,A
double omega( int n ); // 緩和項
void velocity( double t, int N, double *x1, double *x2, double *n1, double *n2, double *l, double *phi, double *kappa, double L, double A, double *V, double *W ); // x,n,l,t,phi,kappa,L --> V,W
void normal_speed( double t, int N, double *kappa, double *phi, double *beta, double *u, double *v, double *V ); // n,phi --> v,V
void tangent_speed( double t, int N, double *l, double *phi, double *kappa, double *v, double *V, double L, double *W ); // l,phi,kappa,v,V,L --> W

void supersaturation( double t, int N, double *x1, double *x2, double *l, double *t1, double *t2, double *n1, double *n2, double *nu, double A, double r_c, double R, double *beta, double **U, double *q, double *u );

double gg( double t, int i, int j, double *l, double *x1, double *x2, double *t1, double *t2, double *n1, double *n2, double a );
double hh( double t, int i, int j, double *l, double *x1, double *x2, double *t1, double *t2, double a, double *beta );
double ii( double t, int j, double *l, double *x1, double *x2, double a, double R, double r_c);

//--------------------main--------------------//

int main(void){

  //---------------------変数--------------------//
  
  int i,j,z,N = 64;
  double t,dt = 0.1 / ( N * N );
  double L,A,L_tmp = 2 * M_PI;
  
  double *X1,*X2;
  double *T1,*T2;
  double *N1,*N2;
  double *kappa;
  X1 = make_vector(Z + 2);
  X2 = make_vector(Z + 2);
  T1 = make_vector(Z + 2);
  T2 = make_vector(Z + 2);
  N1 = make_vector(Z + 2);
  N2 = make_vector(Z + 2);
  kappa = make_vector(Z + 2);


  //----------点を増やすのに使うだけ----------//
  
  double *a0,*a1,*a2,*a3,*a4,*a5;
  double *b0,*b1,*b2,*b3,*b4,*b5;
  
  a0 = make_vector(Z + 2);
  a1 = make_vector(Z + 2);
  a2 = make_vector(Z + 2);
  a3 = make_vector(Z + 2);
  a4 = make_vector(Z + 2);
  a5 = make_vector(Z + 2);
  b0 = make_vector(Z + 2);
  b1 = make_vector(Z + 2);
  b2 = make_vector(Z + 2);
  b3 = make_vector(Z + 2);
  b4 = make_vector(Z + 2);
  b5 = make_vector(Z + 2);

  char file[5000];
  FILE *fp;

  t = 0.0;
  z = 0;

  
  //--------------------初期条件--------------------//
  
  initial_condition(N,X1,X2);
  
  sprintf(file, "./data/yoko_kuro%06d.dat", z);
  fp = fopen(file, "w");

  measure(t,N,X1,X2,&L,&A);
  
  for( i = 0; i <= N; i++ ){
    
    //printf("%f %f\n", x[i], y[i]);
    fprintf(fp, "%f %f %f\n", X1[i], X2[i], A);
    
  }
  
  fclose(fp);
  
  printf("z = %d\n", z);

  
  //--------------------回す--------------------//
  
  while(t < TMAX){
    
    runge_qutta(t,dt,N,X1,X2);
    
    measure(t,N,X1,X2,&L,&A);
    
    printf("L = %f, L_tmp = %f\n", L,L_tmp);

    /*
    //----------点を増やす操作----------//
    
    if( L > 2 * L_tmp ){
      
      for( i = 1; i <= N; i++ ){
	
	a0[i] = X1[i - 1];
	b0[i] = X2[i - 1];
	a1[i] = T1[i - 1];
	b1[i] = T2[i - 1];
	a2[i] = -kappa[i - 1] * N1[i - 1] / 2.0;
	b2[i] = -kappa[i - 1] * N2[i - 1] / 2.0;
	a5[i] = 6 * X1[i] - 6 * X1[i - 1] - 6 * T1[i] + 12 * T1[i - 1] - 11 * kappa[i - 1] * N1[i - 1] / 2.0 - kappa[i] * N1[i] / 2.0;
	b5[i] = 6 * X2[i] - 6 * X2[i - 1] - 6 * T2[i] + 12 * T2[i - 1] - 11 * kappa[i - 1] * N2[i - 1] / 2.0 - kappa[i] * N2[i] / 2.0;
	a4[i] = -2 * a5[i] - 3 * X1[i] + 3 * X1[i - 1] + T1[i] - 4 * T1[i - 1] + 5 * kappa[i - 1] * N1[i - 1];
	b4[i] = -2 * b5[i] - 3 * X2[i] + 3 * X2[i - 1] + T2[i] - 4 * T2[i - 1] + 5 * kappa[i - 1] * N2[i - 1];
	a3[i] = -a4[i] - a5[i] + X1[i] - X1[i - 1] + T1[i - 1] - kappa[i - 1] * N1[i - 1] / 2.0;
	b3[i] = -b4[i] - b5[i] + X2[i] - X2[i - 1] + T2[i - 1] - kappa[i - 1] * N2[i - 1] / 2.0;
	
      }
      connect_double(N,a0,b0);
      connect_double(N,a1,b1);
      connect_double(N,a2,b2);
      connect_double(N,a3,b3);
      connect_double(N,a4,b4);
      connect_double(N,a5,b5);

      
      N = N * 2; // 点の数を2倍に
      
      dt = 0.1 / ( N * N );
      

      for( i = N; i >= 2; i -= 2 ){
	
	X1[i] = X1[i / 2];
	X2[i] = X2[i / 2];
	
      }

      j = 0;

      for( i = 1; i <= N; i += 2 ){

	j += 1;
	
	X1[i] = a0[j] + a1[j] / 2.0 + a2[j] / 4.0 + a3[j] / 8.0 + a4[j] / 16.0 + a5[j] / 32.0;
	X2[i] = b0[j] + b1[j] / 2.0 + b2[j] / 4.0 + b3[j] / 8.0 + b4[j] / 16.0 + b5[j] / 32.0;
	
      }
      connect_double(N,X1,X2);
      
      L_tmp = L;
      
    }
    */
    
    z++;
    t = z * dt;

    //-----------ファイル格納----------//
    
    if( z % 1000 == 0 ){
      
      sprintf(file, "./data/yoko_kuro%06d.dat", z / 1000 );
      fp = fopen(file, "w");

      measure(t,N,X1,X2,&L,&A);
      
      for(i = 0; i <= N; i++){
	
	fprintf(fp, "%f %f %f\n", X1[i], X2[i], L);
	
      }
      
      fclose(fp);
      
    }
    
    printf("z = %d, N = %d\n", z, N);

  }
  
  free(X1);
  free(X2);
  free(T1);
  free(T2);
  free(N1);
  free(N2);
  free(kappa);

  free(a0);
  free(a1);
  free(a2);
  free(a3);
  free(a4);
  free(a5);
  free(b0);
  free(b1);
  free(b2);
  free(b3);
  free(b4);
  free(b5);
  
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

// ルンゲクッタ
void runge_qutta( double t, double dt, int N, double *X1, double *X2 ){
  
  int i;
  double t_temp;
  double *x_temp1,*x_temp2,*F1,*F2;
  double *k11,*k12,*k21,*k22,*k31,*k32,*k41,*k42;
  
  k11 = make_vector(Z + 2);
  k12 = make_vector(Z + 2);
  k21 = make_vector(Z + 2);
  k22 = make_vector(Z + 2);
  k31 = make_vector(Z + 2);
  k32 = make_vector(Z + 2);
  k41 = make_vector(Z + 2);
  k42 = make_vector(Z + 2);
  
  x_temp1 = make_vector(Z + 2);
  x_temp2 = make_vector(Z + 2);
  F1 = make_vector(Z + 2);
  F2 = make_vector(Z + 2);
  
  F(t,N,X1,X2,F1,F2);
  
  for( i = 1; i <= N; i++ ){
    
    k11[i] = F1[i];
    k12[i] = F2[i];
    x_temp1[i] = X1[i] + k11[i] * dt / 2.0;
    x_temp2[i] = X2[i] + k12[i] * dt / 2.0;

  }
  connect_double(N,x_temp1,x_temp2);
  
  t_temp = t + dt / 2.0;
  
  F(t_temp,N,x_temp1,x_temp2,F1,F2);
  
  for( i = 1; i <= N; i++ ){
    
    k21[i] = F1[i];
    k22[i] = F2[i];
    x_temp1[i] = X1[i] + k21[i] * dt / 2.0;
    x_temp2[i] = X2[i] + k22[i] * dt / 2.0;

  }
  connect_double(N,x_temp1,x_temp2);

  F(t_temp,N,x_temp1,x_temp2,F1,F2);
  
  for( i = 1; i <= N; i++ ){
    
    k31[i] = F1[i];
    k32[i] = F2[i];
    x_temp1[i] = X1[i] + k31[i] * dt;
    x_temp2[i] = X2[i] + k32[i] * dt;
    
  }
  connect_double(N,x_temp1,x_temp2);

  t_temp = t + dt;
  
  F(t_temp,N,x_temp1,x_temp2,F1,F2);
  
  for( i = 1; i <= N; i++ ){
    
    k41[i] = F1[i];
    k42[i] = F2[i];
      
    X1[i] = X1[i] + dt * ( k11[i] + 2.0 * k21[i] + 2.0 * k31[i] + k41[i] ) / 6.0;
    X2[i] = X2[i] + ( k12[i] + 2.0 * k22[i] + 2.0 * k32[i] + k42[i] ) * dt / 6.0;

  }
  connect_double(N,X1,X2);

  
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

// 右辺
void F( double t, int N, double *x1, double *x2, double *F1, double *F2 ){
  
  int i;
  
  double *T1,*T2,*N1,*N2;
  T1 = make_vector(Z + 2);
  T2 = make_vector(Z + 2);
  N1 = make_vector(Z + 2);
  N2 = make_vector(Z + 2);
  
  double *V,*W;
  V = make_vector(Z + 2);
  W = make_vector(Z + 2);
  
  ODE_pre(t,N,x1,x2,T1,T2,N1,N2,V,W);
  
  for( i = 1; i <= N; i++ ){
    
    F1[i] = V[i] * N1[i] + W[i] * T1[i];
    F2[i] = V[i] * N2[i] + W[i] * T2[i];
    
  }
  connect_double(N,F1,F2);
  
  free(W);
  free(V);
  free(T1);
  free(T2);
  free(N1);
  free(N2);
  
}

// x --> T,N,V,W
void ODE_pre( double t, int N, double *x1, double *x2, double *T1, double *T2, double *N1, double *N2, double *V, double *W ){

  int i;
  double *l;
  double *t1,*t2,*n1,*n2;
  double *nu,*phi;
  double *kappa;
  double *beta,*u,**U,*q;
  double L,A,r_c,R;
  
  l = make_vector(Z + 2);
  nu = make_vector(Z + 2);
  phi = make_vector(Z + 2);
  kappa = make_vector(Z + 2);
  beta = make_vector(Z + 2);
  u = make_vector(Z + 2);
  q = make_vector(Z + 2);
  U = make_matrix(Z + 2,Z + 2);

  t1 = make_vector(Z + 2);
  t2 = make_vector(Z + 2);
  n1 = make_vector(Z + 2);
  n2 = make_vector(Z + 2);
  
  // T,N
  quantities(t,N,x1,x2,l,t1,t2,n1,n2,T1,T2,N1,N2,nu,phi,kappa);

  // beta,u
  supersaturation(t,N,x1,x2,l,t1,t2,n1,n2,nu,A,r_c,R,beta,U,q,u);
  
  // L,A
  measure(t,N,x1,x2,&L,&A);
  
  // V,W
  velocity(t,N,x1,x2,n1,n2,l,phi,kappa,L,A,V,W);
  
  free(t1);
  free(t2);
  free(n1);
  free(n2);

  for( i = 0; i <= Z; i++ ){
    
    free((void *)U[i]);
    
  }
  free((void *)U);
  free(q);
  free(u);
  free(beta);
  free(kappa);
  free(phi);
  free(nu);
  free(l);
  
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
    
    u = i * 2.0 * M_PI / N;
    
    x1[i] = cos(u);
    x2[i] = sin(u);

  }
  connect_double(N,x1,x2);
  
}

// x --> t,n,T,N,phi,kappa
void quantities( double t, int N, double *x1, double *x2, double *l, double *t1, double *t2, double *n1, double *n2, double *T1, double *T2, double *N1, double *N2, double *nu, double *phi, double *kappa ){
  
  int i;
  double D,I;
  
  for( i = 1; i <= N; i++ ){
    
    l[i] = DIST(x1[i-1],x2[i-1],x1[i],x2[i]);
    
    t1[i] = ( x1[i] - x1[i - 1] ) / l[i];
    t2[i] = ( x2[i] - x2[i - 1] ) / l[i];
    
    n1[i] = t2[i];
    n2[i] = -t1[i];
    
  }
  connect(N,l);
  connect_double(N,t1,t2);
  connect_double(N,n1,n2);

  if( t2[1] >= 0 ){
  
    nu[1] = acos(t1[1]);
    
  }

  else{

    nu[1] = -acos(t1[1]);
    
  }
  
  for( i = 1; i <= N; i++ ){
    
    D = DET(t1[i],t2[i],t1[i+1],t2[i+1]);
    I = ip(t1[i],t2[i],t1[i+1],t2[i+1]);
    
    RANGE_CHECK(I,-1.0,1.0);

    nu[i + 1] = nu[i] + D * acos(I);
    
    if( D >= 0.0 ){
      
      phi[i] = acos(I);
      
    }

    else{
      
      phi[i] = -acos(I);

    }
    
  }
  connect(N,phi);
  nu[0] = nu[1] - ( nu[N + 1] - nu[N] );
  
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

// x,l --> L,A
void measure( double t, int N, double *x1, double *x2, double *L, double *A ){
  
  int i;
  
  *L = 0.0;
  *A = 0.0;
  
  for( i = 1; i <= N; i++ ){
    
    *L += DIST(x1[i],x2[i],x1[i - 1],x2[i - 1]);
    *A += DET(x1[i - 1],x2[i - 1],x1[i],x2[i]);
    
  }
  
  *A = *A / 2.0;
  
}

//緩和項
double omega( int n ){
  
  return 10.0 * n;

}

// x,n,l,t,phi,kappa,L --> V,W
void velocity( double t, int N, double *x1, double *x2, double *n1, double *n2, double *l, double *phi, double *kappa, double L, double A, double *V, double *W ){
  
  int i;
  double *v;
  double *beta;
  double *u;
  
  v = make_vector(Z + 2);
  beta = make_vector(Z + 2);
  u = make_vector(Z + 2);
  
  normal_speed(t,N,kappa,phi,beta,u,v,V);
  tangent_speed(t,N,l,phi,kappa,v,V,L,W);
  
  free(v);
  free(beta);
  free(u);
  
}

// n,phi --> v,V
void normal_speed( double t, int N, double *kappa, double *phi, double *beta, double *u, double *v, double *V ){
  
  int i;

  for( i = 1; i <= N; i++ ){
    
    v[i] = beta_max * u[i];

  }
  connect(N,v);
  
  for( i = 1; i <= N; i++ ){
    
    V[i] = ( v[i] + v[i + 1] ) / ( 2.0 * cos(phi[i] / 2.0) );

  }
  connect(N,V);

}

// l,phi,kappa,v,V,L --> W
void tangent_speed( double t, int N, double *l, double *phi, double *kappa, double *v, double *V, double L, double *W ){
  
  int i;
  double *psi, *PSI;
  double L_dot;
  double a,b,c;
  
  psi = make_vector(Z + 2);
  PSI = make_vector(Z + 2);
  
  psi[1] = 0.0;
  L_dot = 0.0;
  
  for( i = 1; i <= N; i++ ){
    
    L_dot += kappa[i] * v[i] * l[i];
    
  }
  
  for( i = 2; i <= N; i++ ){
    
    psi[i] = ( L_dot / N ) - V[i] * sin(phi[i] / 2.0) - V[i - 1] * sin(phi[i - 1] / 2.0) + ( ( L / N ) - l[i] ) * omega(N);

  }
  
  PSI[1] = psi[1];
  
  for( i = 2; i <= N; i++ ){
    
    PSI[i] = PSI[i - 1] + psi[i];

  }
  
  a = 0.0;
  b = 0.0;
  
  for( i = 1; i <= N; i++ ){
    
    a += PSI[i] / cos(phi[i] / 2.0);
    b += 1.0 / cos(phi[i] / 2.0);
    
  }
  
  c = -a / b;
  
  for( i = 1; i <= N; i++ ){
    
    W[i] = ( PSI[i] + c ) / cos(phi[i] / 2.0);

  }
  connect(N,W);

  free(PSI);
  free(psi);
  
}

// 過飽和度
void supersaturation( double t, int N, double *x1, double *x2, double *l, double *t1, double *t2, double *n1, double *n2, double *nu, double A, double r_c, double R, double *beta, double **U, double *q, double *u ){

  int i,j,k,ep;
  double z1, z2;
  double gamma, tmp, amax;
  double dx_sim; // 数値積分の分割幅
  double dx_pi; // 数値積分の分割幅

  r_c = sqrt(A / M_PI);
  
  R = 6.5 * r_c;
  
  if( t == 0.0 ){

    // printf("a\n");

    for( i = 1; i <= N; i++ ){
      
      if( nu[i] == 0 || nu[i] == ( M_PI / 3.0 ) || nu[i] == ( 2 * M_PI / 3.0 ) || nu[i] == M_PI || nu[i] == ( 4 * M_PI / 3.0 ) || nu[i] == ( 5 * M_PI / 3.0 ) || nu[i] == 2 * M_PI ){
	
	beta[i] = ( 1 - delta_beta ) * beta_max;
	
      }
      
      else{
	
	beta[i] = beta_max * 2 * x_s * tan(nu[i]) * tanh(e / ( 2 * x_s * tan(nu[i]) )) / e;
	
      }
      
    }
    
  }
  
  else{
    
    // printf("b\n");
    
    for( i = 1; i <= N; i++ ){
      
      if( nu[i] == 0 || nu[i] == ( M_PI / 3.0 ) || nu[i] == ( 2 * M_PI / 3.0 ) || nu[i] == M_PI || nu[i] == ( 4 * M_PI / 3.0 ) || nu[i] == ( 5 * M_PI / 3.0 ) || nu[i] == 2 * M_PI ){
	
	beta[i] = beta_max * u[i] * tanh(sigma_star / u[i]) / sigma_star;

      }
      
      else{
	
	beta[i] = beta_max * 2 * x_s * tan(nu[i]) * tanh(e / ( 2 * x_s * tan(nu[i]) )) / e;
	
      }
      
    }
    
  }
  connect(N,beta);
  
  for( i = 1; i <= N; i++ ){
    
    for( j = 1; j <= N; j++ ){
      
      if( j == i ){
	
	U[i][j] = 0.5 - ( ( k_B * T * beta_max ) / ( 2 * M_PI * v_c * p_e * E ) ) * l[i] * ( 1 - log(l[i] / 2.0) );
	
      }

      else{
	
	//----------数値積分----------//
	
	dx_sim = l[i] / N;
	
	U[i][j] = dx_sim * ( gg(t,i,j,l,x1,x2,t1,t2,n1,n2,0) + gg(t,i,j,l,x1,x2,t1,t2,n1,n2,dx_sim) ) / 2.0;
	
	for( k = 1; k < N; k++ ){
	  
	  z1 = k * dx_sim;
	  z2 = ( k + 1 ) * dx_sim;
	  
	  U[i][j] = U[i][j] + dx_sim * ( gg(t,i,j,l,x1,x2,t1,t2,n1,n2,z1) + gg(t,i,j,l,x1,x2,t1,t2,n1,n2,z2) ) / 2.0;
	  
	}
	
	U[i][j] = U[i][j] + dx_sim * ( hh(t,i,j,l,x1,x2,t1,t2,0,beta) + hh(t,i,j,l,x1,x2,t1,t2,dx_sim,beta) ) / 2.0;
	
	for( k = 1; k < N; k++ ){
	  
	  z1 = k * dx_sim;
	  z2 = ( k + 1 ) * dx_sim;
	  
	  U[i][j] = U[i][j] + dx_sim * ( hh(t,i,j,l,x1,x2,t1,t2,z1,beta) + hh(t,i,j,l,x1,x2,t1,t2,z2,beta) ) / 2.0;
	  
	}
	
      }
      
    }
    
  }

  
  for( j = 1; j <= N; j++ ){
    
      dx_pi = 2 * M_PI / N;
      
      q[j] = dx_pi * ( ii(t,j,l,x1,x2,0,R,r_c) + ii(t,j,l,x1,x2,dx_pi,R,r_c) ) / 2.0;
      
      for( k = 1; k < N; k++ ){
	
	z1 = k * dx_pi;
	z2 = ( k + 1 ) * dx_pi;
	
	q[j] = q[j] + dx_pi * ( ii(t,j,l,x1,x2,z1,R,r_c) + ii(t,j,l,x1,x2,z2,R,r_c) ) / 2.0;
	
      }
      
  }
  
  /*
  for( j = 1; j <= N; j++ ){

    q[i] = -2 * sigma_infty;
    
  }
  */

  //----------ガウスの消去法----------//
  
  //----------消去----------//
  
  for( k = 1; k <= N - 1; k++ ){
    
    amax = fabs(U[k][k]);
    ep = k;
    
    for( i = k + 1; i <= N; i++ ){
      
      if( fabs(U[i][k]) > amax ){
	
	amax = fabs(U[i][k]);
	ep = i;
	
      }
      
    }
    
    if( amax < eps_sim ){
      
      printf("入力した行列は正則ではない\n");
      
    }
    
    if( ep != k ){
      
      for( j = k; j <= N; j++ ){
	
	tmp = U[k][j];
	U[k][j] = U[ep][j];
	U[ep][j] = tmp;
	
      }
      
      tmp = q[k];
      q[k] = q[ep];
      q[ep] = tmp;
      
    }
    
    for( i = k + 1; i <= N; i++ ){
      
      gamma = - U[i][k] / U[k][k];
      
      for( j = k + 1; j <= N; j++ ){
	
	U[i][j] = U[i][j] + gamma * U[k][j];
	
      }
      
      q[i] = q[i] + gamma * q[k];
      
    }
    
  }
  
  u[N] = q[N] / U[N][N];
  
  for( k = N - 1; k >= 1; k-- ){
    
    tmp = q[k];
    
    for( j = k + 1; j <= N; j++ ){
      
      tmp = tmp - U[k][j] * u[j];
      
    }
    
    u[k] = tmp / U[k][k];
    
  }
  connect(N,u);
  
}

double gg( double t, int i, int j, double *l, double *x1, double *x2, double *t1, double *t2, double *n1, double *n2, double a ){

  return (

	  ( ( x1[j - 1] - ( x1[i - 1] - a * t1[i] ) ) * n1[i] + ( x2[j - 1] - ( x2[i - 1] - a * t2[i] ) ) * n2[i] ) / ( 2 * M_PI * ( sqrt( ( x1[j - 1] - ( x1[i - 1] - a * t1[i] ) ) * ( x1[j - 1] - ( x1[i - 1] - a * t1[i] ) ) + ( x2[j - 1] - ( x2[i - 1] - a * t2[i] ) ) * ( x2[j - 1] - ( x2[i - 1] - a * t2[i] ) ) + epsl * epsl ) ) )
	  
	  );

}

double hh( double t, int i, int j, double *l, double *x1, double *x2, double *t1, double *t2, double a, double *beta ){

  return (

	  ( ( k_B * T * beta_max ) / ( 2 * M_PI * E * p_e * v_c ) ) * log( ( x1[j - 1] - ( x1[i - 1] - a * t1[i] ) ) * ( x1[j - 1] - ( x1[i - 1] - a * t1[i] ) ) + ( x2[j - 1] - ( x2[i - 1] - a * t2[i] ) ) * ( x2[j - 1] - ( x2[i - 1] - a * t2[i] ) ) + epsl * epsl )
	  
	  );

}

double ii( double t, int j, double *l, double *x1, double *x2, double a, double R, double r_c){

  return (

	  ( -( log( sqrt( ( x1[j - 1] - R * cos(a) ) * ( x1[j - 1] - R * cos(a) ) + ( x2[j - 1] - R * sin(a) ) * ( x2[j - 1] - R * sin(a) ) ) + epsl * epsl ) / ( 2 * M_PI * R * log(R / r_c) ) )
	  + ( ( x1[j - 1] * cos(a) + x2[j - 1] * sin(a) - R ) / ( 2 * M_PI * ( x1[j - 1] - R * cos(a) ) * ( x1[j - 1] - R * cos(a) ) + ( x2[j - 1] - R * sin(a) ) * ( x2[j - 1] - R * sin(a) ) ) ) ) * R * sigma_infty
	  
	  );

}

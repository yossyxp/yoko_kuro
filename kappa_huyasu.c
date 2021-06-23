
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

//--------------------定数--------------------//

// #define N ( 600 )
#define Z ( 2048 )
#define TMAX ( 4.4 )
// #define dt ( 0.1 / ( N * N ) )
#define RANGE_CHECK(x, xmin, xmax) ( x = ( x < xmin ? xmin : ( x < xmax ?  x : xmax ) ) );


//--------------------関数--------------------//

double* make_vector(int n); // 領域の確保

void connect(int N, double *x); // 閉曲線のつなぎ
void connect_double(int N, double* x, double *y); // つなぎ2個
double ip(double x1, double x2, double y1, double y2); // 内積
double DIST(double x1, double x2, double y1, double y2); // |y - x|
double DET(double x1, double x2, double y1, double y2); // 行列式
void runge_qutta(double t, double dt, int N, double *X1, double *X2); // ルンゲクッタ
void F(double t, int N, double *x1, double *x2, double *F1, double *F2);
void ODE_pre(double t, int N, double *x1, double *x2, double *T1, double *T2, double *N1, double *N2, double *V, double *W); // x --> T,N,V,W
void initial_condition(int N, double *x1, double *x2); // 初期条件
void quantities(double t, int N, double *x1, double *x2, double *l, double *t1, double *t2, double *n1, double *n2, double *T1, double *T2, double *N1, double *N2, double *nu, double *phi, double *kappa); //x --> t,n,T,N,phi,kappa
void measure(double t, int N, double *x1, double *x2, double *L, double *A); // x,l --> L,A
void increase( double t, int N, double *x1, double *x2, double *T1, double *T2, double *N1, double *N2, double *kappa, double L, double A, double L_tmp, double dt );
double omega(int n); // 緩和項
void velocity(double t, int N, double *x1, double *x2, double *n1, double *n2, double *l, double *phi, double *kappa, double L, double A, double *V, double *W); // x,n,l,t,phi,kappa,L --> V,W
void normal_speed(double t, int N, double *kappa, double *phi, double *v, double *V); // n,phi --> v,V
void tangent_speed(double t, int N, double *l, double *phi, double *kappa, double *v, double *V, double L, double *W); // l,phi,kappa,v,V,L --> W

//--------------------main--------------------//

int main(void){

  //-----------変数----------//
  
  int i,j,z,N = 128;
  double t,dt = 0.1 / ( N * N );
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
  
  double L,A,L_tmp = 2 * M_PI;

  char file[5000];
  FILE *fp;

  t = 0.0;
  z = 0;

  //----------初期条件----------//
  
  initial_condition(N,X1,X2);
  
  sprintf(file, "./data/yoko_kuro%06d.dat", z);
  fp = fopen(file, "w");

  measure(t,N,X1,X2,&L,&A);
  
  for( i = 0; i <= N; i++ ){
    
    //printf("%f %f\n", x[i], y[i]);
    fprintf(fp, "%f %f %f\n", X1[i], X2[i], L);
    
  }
  
  fclose(fp);
  
  printf("z = %d\n", z);

  
  //--------------------回そう--------------------//
  
  while( t < TMAX ){
    
    runge_qutta(t,dt,N,X1,X2);
    
    measure(t,N,X1,X2,&L,&A);
    
    printf("L = %f, L_tmp = %f\n", L,L_tmp);

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
      
      N = N * 2;
      
      dt = 0.1 / ( N * N );

      // ここが入らない
      for( i = N; i >= 2; i-=2 ){
	
	X1[i] = X1[i / 2];
	X2[i] = X2[i / 2];

	printf("babb\n");

	printf("%f %f %d\n",X1[i],X2[i],i);
	
      }

      j = 0;

      // ここが 1 --> 1 3 --> 2 5 --> 3
      for( i = 1; i <= N; i += 2 ){

	j += 1;
	
	X1[i] = a0[j] + a1[j] / 2.0 + a2[j] / 4.0 + a3[j] / 8.0 + a4[j] / 16.0 + a5[j] / 32.0;
	X2[i] = b0[j] + b1[j] / 2.0 + b2[j] / 4.0 + b3[j] / 8.0 + b4[j] / 16.0 + b5[j] / 32.0;

	printf("%f %f %d\n",X1[i],X2[i],i);
	
      }
      connect_double(N,X1,X2);
      
      L_tmp = L;
      
    }
    
    
    z++;
    t = z * dt;
    
    if( z % 10000 == 0 ){
      
      sprintf(file, "./data/yoko_kuro%06d.dat", z / 10000 );
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
  
  for(i = 1; i <= N; i++){
    
    k11[i] = F1[i];
    k12[i] = F2[i];
    x_temp1[i] = X1[i] + k11[i] * dt / 2.0;
    x_temp2[i] = X2[i] + k12[i] * dt / 2.0;

  }
  connect_double(N,x_temp1,x_temp2);
  
  t_temp = t + dt / 2.0;
  
  F(t_temp,N,x_temp1,x_temp2,F1,F2);
  
  for(i = 1; i <= N; i++){
    
    k21[i] = F1[i];
    k22[i] = F2[i];
    x_temp1[i] = X1[i] + k21[i] * dt / 2.0;
    x_temp2[i] = X2[i] + k22[i] * dt / 2.0;

  }
  connect_double(N,x_temp1,x_temp2);

  F(t_temp,N,x_temp1,x_temp2,F1,F2);
  
  for(i = 1; i <= N; i++){
    
    k31[i] = F1[i];
    k32[i] = F2[i];
    x_temp1[i] = X1[i] + k31[i] * dt;
    x_temp2[i] = X2[i] + k32[i] * dt;
    
  }
  connect_double(N,x_temp1,x_temp2);

  t_temp = t + dt;
  
  F(t_temp,N,x_temp1,x_temp2,F1,F2);
  
  for(i = 1; i <= N; i++){
    
    k41[i] = F1[i];
    k42[i] = F2[i];
      
    X1[i] = X1[i] + (k11[i] + 2.0*k21[i] + 2.0*k31[i] + k41[i])*dt/6.0;
    X2[i] = X2[i] + (k12[i] + 2.0*k22[i] + 2.0*k32[i] + k42[i])*dt/6.0;

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
  
  for(i = 1; i <= N; i++){
    
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
  
  double *l;
  double *t1,*t2,*n1,*n2;
  double *nu;
  double *phi;
  double *kappa;
  double L,A;
  
  l = make_vector(Z + 2);
  nu = make_vector(Z + 2);
  phi = make_vector(Z + 2);
  kappa = make_vector(Z + 2);

  t1 = make_vector(Z + 2);
  t2 = make_vector(Z + 2);
  n1 = make_vector(Z + 2);
  n2 = make_vector(Z + 2);
  
  // T,N
  quantities(t,N,x1,x2,l,t1,t2,n1,n2,T1,T2,N1,N2,nu,phi,kappa);
  
  // L,A
  measure(t,N,x1,x2,&L,&A);
  
  // V,W
  velocity(t,N,x1,x2,n1,n2,l,phi,kappa,L,A,V,W);
  
  free(t1);
  free(t2);
  free(n1);
  free(n2);
  
  free(kappa);
  free(nu);
  free(phi);
  free(l);
  
}

// 初期条件
void initial_condition( int N, double *x1, double *x2 ){
  
  int i;
  double u;
  double a1,a2,a3;

  /*
  for(i = 1; i <= N; i++){
    
    u = i * 1.0 / N;
    
    a1 = 1.8 * cos(2.0 * M_PI * u);
    a2 = 0.2 + sin(M_PI * u) * sin(6.0 * M_PI * u) * sin(2.0 * a1);
    a3 = 0.5 * sin(2.0 * M_PI * u) + sin(a1) + a2 * sin(2.0 * M_PI * u);
    
    x1[i] = 2.0 * 0.5 * a1;
    x2[i] = 2.0 * 0.54 * a3;

  }
  connect_double(N,x1,x2);
  */
  
  for(i = 1; i <= N; i++){
    
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

    printf("%f\n",nu[i]);
    
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

// 増やす
void increase( double t, int N, double *x1, double *x2, double *T1, double *T2, double *N1, double *N2, double *kappa, double L, double A, double L_tmp, double dt ){

  int i;
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
  
  printf("abc\n");
  
  for( i = 1; i <= N; i++ ){
      
    a0[i] = x1[i - 1];
    b0[i] = x2[i - 1];
    a1[i] = T1[i - 1];
    b1[i] = T2[i - 1];
    a2[i] = -kappa[i - 1] * N1[i - 1] / 2.0;
    b2[i] = -kappa[i - 1] * N2[i - 1] / 2.0;
    a5[i] = 6 * x1[i] - 6 * x1[i - 1] - 6 * T1[i] + 12 * T1[i - 1] - 11 * kappa[i - 1] * N1[i - 1] / 2.0 - kappa[i] * N1[i] / 2.0;
    b5[i] = 6 * x2[i] - 6 * x2[i - 1] - 6 * T2[i] + 12 * T2[i - 1] - 11 * kappa[i - 1] * N2[i - 1] / 2.0 - kappa[i] * N2[i] / 2.0;
    a4[i] = -2 * a5[i] - 3 * x1[i] + 3 * x1[i - 1] + T1[i] - 4 * T1[i - 1] + 5 * kappa[i - 1] * N1[i - 1];
    b4[i] = -2 * b5[i] - 3 * x2[i] + 3 * x2[i - 1] + T2[i] - 4 * T2[i - 1] + 5 * kappa[i - 1] * N2[i - 1];
    a3[i] = -a4[i] - a5[i] + x1[i] - x1[i - 1] + T1[i - 1] - kappa[i - 1] * N1[i - 1] / 2.0;
    b3[i] = -b4[i] - b5[i] + x2[i] - x2[i - 1] + T2[i - 1] - kappa[i - 1] * N2[i - 1] / 2.0;
    
  }
  connect_double(N,a0,b0);
  connect_double(N,a1,b1);
  connect_double(N,a2,b2);
  connect_double(N,a3,b3);
  connect_double(N,a4,b4);
  connect_double(N,a5,b5);
  
  N = N * 2;
  
  dt = 0.1 / ( N * N ); 
  
  for( i = N; i <= 2; i-- ){
    
    x1[i] = x1[i / 2];
    x2[i] = x2[i / 2];
    
  }
  
  for( i = 1; i <= N; i += 2 ){
    
    x1[i] = a0[i] + a1[i] / 2.0 + a2[i] / 4.0 + a3[i] / 8.0 + a4[i] / 16.0 + a5[i] / 32.0;
    x2[i] = b0[i] + b1[i] / 2.0 + b2[i] / 4.0 + b3[i] / 8.0 + b4[i] / 16.0 + b5[i] / 32.0;
    
  }
  connect_double(N,x1,x2);
  

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
  
}

// h = h(t)
double height( double t ){
  
  return exp(t);

}

double height_t( double t ){
  
  return exp(t);
  
}

//緩和項
double omega( int n ){
  
  return 10.0 * n;

}

// x,n,l,t,phi,kappa,L --> V,W
void velocity( double t, int N, double *x1, double *x2, double *n1, double *n2, double *l, double *phi, double *kappa, double L, double A, double *V, double *W ){
  
  int i;
  double *v;
  
  v = make_vector(Z + 2);
  
  normal_speed(t,N,kappa,phi,v,V);
  tangent_speed(t,N,l,phi,kappa,v,V,L,W);
  
  free(v);
  
}

// n,phi --> v,V
void normal_speed( double t, int N, double *kappa, double *phi, double *v, double *V )
{
  int i;

  for( i = 1; i <= N; i++ ){
    
    v[i] = 0.4;

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


#include<stdio.h>
#include<math.h>
#include<stdlib.h>

//--------------------定数--------------------//

#define N ( 1024 ) // 配列の要素数
#define alpha ( 2.0 / sqrt(3) )
#define dt ( 0.0005 ) // 時間刻み
#define eps ( 1.0e-2 ) // 辺の結合に用いるeps
#define K ( 1000.0 ) // パラメータ

//----------定数(ラプラス)----------//

#define s ( 250 ) // 空間分割数(x軸)
#define m ( 250 ) // 空間分割数(y軸)
#define dx ( 0.16 ) // 空間刻み幅(x軸)
#define dy ( 0.16 ) // 空間刻み幅(y軸)
#define epsl ( 1.0e-15 ) // 場合分け回避に用いるeps

//----------定数(連立)----------//

#define eps_sim ( 1.0e-40 )

//--------------------関数宣言--------------------//

void connect(int n, double *a); // 閉曲線の繋ぎ
void normal(int n, double *t1, double *t2, double *n1, double *n2); // 法線ベクトル
void height(int n, double *x, double *y, double *n1, double *n2, double *h); // 高さ関数
void tangent_angle_pre(int n, double *t1, double *t2, double *D); // 接線ベクトルの行列式の符号
void tangent_angle(int n, double *t1, double *t2, double *D, double *nu); // 接線角度
void outside_angle(int n, double *nu, double *phi); // 外角
void transition_number_pre(int n, double *n1, double *n2, double *sigma); // 法線ベクトルの行列式の符号
void transition_number(int n, double *sigma, double *chi); // 遷移数
void new_coordinate(int n, double *h, double *t1, double *t2, double *phi, double *x, double *y); // 新たな座標
void new_length(int n, double *h, double *phi, double *l); // 新たな辺の長さ
double ff( int i, double *p, double t ); // 高さ関数の発展方程式の右辺
double gg( int i, int j, double *l, double *x, double *y, double *t1, double *t2, double a );

//--------------------main文--------------------//

int main(void){
  
  
  //----------変数----------//
  
  // for文等に使う
  int i, j, k, n, z, ip;
  double t = 0.0;
  
  // 許容多角形
  double *x, *y; // 曲線上の点
  double *l; // 辺の長さ
  double l_max; // 辺の最大値
  double lc = 10.0; // critical length
  double *t1, *t2; // 接線ベクトル
  double *n1, *n2; // 法線ベクトル
  double *h; // 高さ関数
  double *nu; // 接戦角度
  double *phi; // 外角
  double *D; // 接線ベクトルの行列式の符号
  double *sigma; // 法線ベクトルの行列式の符号
  double *chi; // 遷移数
  double *kappa; // 曲率
  double k1, k2, k3, k4; // ルンゲクッタ
  double x_temp, y_temp; // 2分割する辺の真ん中の座標
  double x_temp1, y_temp1; // 3分割する辺の1:2の座標
  double x_temp2, y_temp2; // 3分割する辺の2:1の座標
  
  // uについての計算
  double x1, x2; // 数値積分
  double dx_sim; // 数値積分の分割幅
  double xl, yl; // 空間全体の座標
  double beta, tmp, amax; 
  double **T; // 連立方程式の行列
  double *p, *q; // 連立方程式の解と右辺
  double A, c; // 一時的な入れ物
  double *B, *d; // 一時的な入れ物
  double **u; // 過飽和度

  double r, v, w;
  
  // file出力
  char file[5000];
  char file2[5000];
  FILE *fp;
  FILE *fp2;
  
  
  //--------------------動的変数の確保--------------------//
  
  if( ( x = malloc( N * sizeof(double *) ) ) == NULL ){
    
    printf("メモリが確保できません\n");
    exit(1);
    
  }
  
  if( ( y = malloc( N * sizeof(double *) ) ) == NULL ){
    
    printf("メモリが確保できません\n");
    exit(1);
    
  }
  
  if( ( l = malloc( N * sizeof(double *) ) ) == NULL ){
    
    printf("メモリが確保できません\n");
    exit(1);
    
  }
  
  if( ( t1 = malloc( N * sizeof(double *) ) ) == NULL ){
    
    printf("メモリが確保できません\n");
    exit(1);
    
  }
  
  if( ( t2 = malloc( N * sizeof(double *) ) ) == NULL ){
    
    printf("メモリが確保できません\n");
    exit(1);
    
  }
  
  if( ( n1 = malloc( N * sizeof(double *) ) ) == NULL ){
    
    printf("メモリが確保できません\n");
    exit(1);
    
  }
  
  if( ( n2 = malloc( N * sizeof(double *) ) ) == NULL ){
    
    printf("メモリが確保できません\n");
    exit(1);
    
  }
  
  if( ( h = malloc( N * sizeof(double *) ) ) == NULL ){
    
    printf("メモリが確保できません\n");
    exit(1);
    
  }
  
  if( ( nu = malloc( N * sizeof(double *) ) ) == NULL ){
    
    printf("メモリが確保できません\n");
    exit(1);
    
  }
  
  if( ( phi = malloc( N * sizeof(double *) ) ) == NULL ){
    
    printf("メモリが確保できません\n");
    exit(1);
    
  }
  
  if( ( D = malloc( N * sizeof(double *) ) ) == NULL ){
    
    printf("メモリが確保できません\n");
    exit(1);
    
  }
  
  if( ( sigma = malloc( N * sizeof(double *) ) ) == NULL ){
    
    printf("メモリが確保できません\n");
    exit(1);
    
  }
  
  if( ( chi = malloc( N * sizeof(double *) ) ) == NULL ){
    
    printf("メモリが確保できません\n");
    exit(1);
    
  }
  
  if( ( kappa = malloc( N * sizeof(double *) ) ) == NULL ){
    
    printf("メモリが確保できません\n");
    exit(1);
    
  }
  
  if( ( T = malloc( N * sizeof(double *) ) ) == NULL ){
    
    printf("メモリが確保できません\n");
    exit(1);
    
  }
  
  for( i = 0; i <= N; i++ ){
    
    T[i] = malloc( N * sizeof(double) );
    
  }
  
  if( ( p = malloc( N * sizeof(double *) ) ) == NULL ){
    
    printf("メモリが確保できません\n");
    exit(1);
    
  }
  
  if( ( q = malloc( N * sizeof(double *) ) ) == NULL ){
    
    printf("メモリが確保できません\n");
    exit(1);
    
  }
  
  if( ( B = malloc( N * sizeof(double *) ) ) == NULL ){
    
    printf("メモリが確保できません\n");
    exit(1);
    
  }
  
  if( ( d = malloc( N * sizeof(double *) ) ) == NULL ){
    
    printf("メモリが確保できません\n");
    exit(1);
    
  }
  
  if( ( u = malloc( s * sizeof(double *) ) ) == NULL ){
    
    printf("メモリが確保できません\n");
    exit(1);
    
  }
  
  for( i = 0; i <= s; i++ ){
    
    u[i] = malloc( m * sizeof(double) );
    
  }
  
  
  //----------初期値(正六角形)---------//
  
  n = 6;
  
  x[0] = sqrt(3) / 2.0;
  x[1] = 0.0;
  x[2] = -sqrt(3) / 2.0;
  x[3] = -sqrt(3) / 2.0;
  x[4] = 0.0;
  x[5] = sqrt(3) / 2.0;
  x[n] = x[0];
  x[n + 1] = x[1];
  
  y[0] = 1.0 / 2.0;
  y[1] = 1.0;
  y[2] = 1.0 / 2.0;
  y[3] = -1.0 / 2.0;
  y[4] = -1.0;
  y[5] = -1.0 / 2.0;
  y[n] = y[0];
  y[n + 1] = y[1];
  
  
  //----------初期値の出力-----------//
  
  for( z = 0; z <= 0; z++ ){
    
    sprintf(file, "./data1/snow%06d.dat", z);
    fp = fopen(file, "w");
    
    for( i = 0; i <= n; i++ ){
      
      //printf("%f %f\n", x[i], y[i]);
      fprintf(fp, "%f %f %f\n", x[i], y[i], 0.0);
      
    }
    
    fclose(fp);
    
    printf("t = %d\n", z);
    
  }
  
  
  //----------準備----------//
  
  // 辺の長さ
  for( i = 1; i <= n; i++ ){
    
    l[i] = sqrt( ( x[i] - x[i - 1] ) * ( x[i] - x[i - 1] ) + ( y[i] - y[i - 1] ) * ( y[i] - y[i - 1] ) );
    
  }
  connect(n,l);
  
  // 接線ベクトル
  for( i = 1; i <= n; i++ ){
    
    t1[i] = ( x[i] - x[i - 1] ) / l[i];
    t2[i] = ( y[i] - y[i - 1] ) / l[i];
    
  }
  connect(n,t1);
  connect(n,t2);
  
  // 外向き法線ベクトル
  normal(n,t1,t2,n1,n2);
  connect(n,n1);
  connect(n,n2);
  
  // 高さ関数
  height(n,x,y,n1,n2,h);
  connect(n,h);
  
  //接線角度の準備
  tangent_angle_pre(n,t1,t2,D);
  connect(n,D);
  
  // 接線角度
  tangent_angle(n,t1,t2,D,nu);
  connect(n,nu);
  
  // 外角
  outside_angle(n,nu,phi);
  connect(n,phi);
  
  //遷移数の準備
  transition_number_pre(n,n1,n2,sigma);
  connect(n,sigma);
  
  //遷移数
  transition_number(n,sigma,chi);
  connect(n,chi);
  
  // 曲率
  for( i = 1; i <= n; i++ ){
    
    kappa[i] = alpha * chi[i] / l[i];
    
  }
  connect(n,kappa);
  

  //----------critical length----------//

  v = -1.0 / 6.0;
  
  r = ( l[1] * l[1] * log(l[1]) - ( ( 3 * l[1] * l[1] ) / 2.0 ) ) * v / ( 2 * M_PI * l[1] );
  
  w = -6 * ( ( ( l[1] / 3.0 ) * ( l[1] / 3.0 ) * log(( l[1] / 3.0 )) - ( ( 3 * ( l[1] / 3.0 ) * ( l[1] / 3.0 ) ) / 2.0 ) ) * v ) / ( 2 * M_PI ) + r;

  lc = -16 / ( sqrt(3) * K * w );

  printf("lc = %f\n", lc);
  
  //--------------------全体を回す--------------------//
  
  t = 0.0;
  
  for( z = 1; z <= 500000; z++ ){ //10000000
    
    t += dt;
    
    
    //--------------------連立方程式-------------------//
    
    //----------数値積分----------//
    
    for( i = 1; i <= n; i++ ){
      
      for( j = 1; j <= n; j++ ){
	
	if( j == i ){
	  
	  T[i][j] = l[j] * l[j] * log(l[j]) - ( ( 3 * l[j] * l[j] ) / 2.0 );
	  
	}
	
	else{
	  
	  dx_sim = l[i] / N;
	  
	  T[i][j] = dx_sim * ( gg(i,j,l,x,y,t1,t2,0) + gg(i,j,l,x,y,t1,t2,dx_sim) ) / 2.0;
	  
	  for( k = 1; k < N; k++ ){
	    
	    x1 = k * dx_sim;
	    x2 = ( k + 1 ) * dx_sim;
	    
	    T[i][j] = T[i][j] + dx_sim * ( gg(i,j,l,x,y,t1,t2,x1) + gg(i,j,l,x,y,t1,t2,x2) ) / 2.0;
	    
	  }
	  
	}
	
      }
      
    }
    
    
    //----------ガウスの消去法----------//
    
    //----------準備----------//
    
    for( i = 1; i <= n; i++ ){
      
      for( j = 1; j <= n; j++ ){
	
	T[i][j] = - T[i][j] / ( 2 * M_PI );
	
      }
      
    }
    
    for( i = 1; i <= n; i++ ){
      
      T[i][n + 1] = l[i];
      
    }
    
    for( j = 1; j <= n; j++ ){
      
      T[n + 1][j] = l[j];
      
    }
    
    T[n + 1][n + 1] = 0.0;
    
    
    for( i = 1; i <= n; i++ ){
      
      q[i] = l[i] * kappa[i];
      
    }
    
    q[n + 1] = -K;
    
    
    //----------消去----------//
    
    for( k = 1; k <= n; k++ ){
      
      amax = fabs(T[k][k]);
      ip = k;
      
      for( i = k + 1; i <= n + 1; i++ ){
	
	if( fabs(T[i][k]) > amax ){
	  
	  amax = fabs(T[i][k]);
	  ip = i;
	  
	}
	
      }
      
      if( amax < eps_sim ){
	
	printf("入力した行列は正則ではない\n");
	
      }
      
      if( ip != k ){
	
	for( j = k; j <= n + 1; j++ ){
	  
	  tmp = T[k][j];
	  T[k][j] = T[ip][j];
	  T[ip][j] = tmp;
	  
	}
	
	tmp = q[k];
	q[k] = q[ip];
	q[ip] = tmp;
	
      }
      
      for( i = k + 1; i <= n + 1; i++ ){
	
	beta = - T[i][k] / T[k][k];
	
	for( j = k + 1; j <= n + 1; j++ ){
	  
	  T[i][j] = T[i][j] + beta * T[k][j];
	  
	}
	
	q[i] = q[i] + beta * q[k];
	
      }
      
    }
    
    p[n + 1] = q[n + 1] / T[n + 1][n + 1];
    
    for( k = n; k >= 1; k-- ){
      
      tmp = q[k];
      
      for( j = k + 1; j <= n + 1; j++ ){
	
	tmp = tmp - T[k][j] * p[j];
	
      }
      
      p[k] = tmp / T[k][k];
      
    }
    
    
    //----------p_jとalphaからuを求める----------//
    
    for( i = 0; i <= s; i++ ){
      
      xl = i * dx - 20.0;
      
      for( j = 0; j <= m; j++ ){
	
	yl = j * dy - 20.0;
	
	A = ( x[0] - xl ) * t2[1] - ( y[0] - yl ) * t1[1];
	c = ( x[0] - xl ) * t1[1] + ( y[0] - yl ) * t2[1];
	
	u[i][j] = ( ( l[1] + c ) * log( sqrt( ( l[1] + c ) * ( l[1] + c ) + A * A + epsl * epsl ) ) - l[1] + fabs( A ) * atan( ( l[1] + c ) / sqrt( A * A + epsl * epsl ) ) - c * log( sqrt( c * c + A * A + epsl * epsl ) ) - fabs( A ) * atan( c / sqrt( A * A + epsl * epsl ) ) ) * p[1];
	
	for( k = 2; k <= n; k++ ){
	  
	  B[k] = ( x[k - 1] - xl ) * t2[k] - ( y[k - 1] - yl ) * t1[k];
	  d[k] = ( x[k - 1] - xl ) * t1[k] + ( y[k - 1] - yl ) * t2[k];
	  
	  u[i][j] = u[i][j] + ( ( l[k] + d[k] ) * log( sqrt( ( l[k] + d[k] ) * ( l[k] + d[k] ) + B[k] * B[k] + epsl * epsl ) ) - l[k] + fabs( B[k] ) * atan( ( l[k] + d[k] ) / sqrt( B[k] * B[k] + epsl * epsl ) ) - d[k] * log( sqrt( d[k] * d[k] + B[k] * B[k] + epsl * epsl ) ) - fabs( B[k] ) * atan( d[k] / sqrt( B[k] * B[k] + epsl * epsl ) ) ) * p[k];
	  
	}
	
      }
      
    }
    
    
    for( i = 0; i <= s; i++ ){
      
      for( j = 0; j <= m; j++ ){
	
	u[i][j] = -( u[i][j] / ( 2 * M_PI ) ) + p[n + 1];
	
      }
      
    }
    
    
    //----------出力(u)----------//
    
    if( z % 2500 == 0 ){ // 50000
      
      sprintf(file2, "./data2/sim%06d.dat", z / 2500 );
      //sprintf(file2, "./data2/sim%06d.dat", z );
      fp2 = fopen(file2,"w");
      
      
      for( i = 0; i <= s; i++ ){
	
	xl = dx * i - 20.0;
	
	for( j = 0; j <= m; j++ ){ 
	  
	  yl = dy * j - 20.0;
	  
	  fprintf(fp2, "%f %f %f\n", xl, yl, u[i][j] / K);
	  
	}
	
	fprintf(fp2, "\n");
	
      }
      
      fclose(fp2);
      
    }
    
    
    //--------------------ODEの計算--------------------//
    
    // ルンゲクッタ
    for( i = 1; i <= n; i++ ){
      
      k1 = dt * ff(i,p,t);
      p[i] = p[i] + k1 / 2.0;
      k2 = dt * ff(i,p,t + dt / 2);
      p[i] = p[i] + k2 / 2.0;
      k3 = dt * ff(i,p,t + dt / 2);
      p[i] = p[i] + k3;
      k4 = dt * ff(i,p,t + dt);
      
      h[i] = h[i] + dt * ( k1 + 2 * k2 + 2 * k3 + k4 ) / 6.0;
      
    }
    connect(n,h);
    
    
    // 新たな座標
    new_coordinate(n,h,t1,t2,phi,x,y);
    connect(n,x);
    connect(n,y);
    
    
    //----------出力(x,y)----------//
    
    if( z % 2500 == 0 ){ // 50000
      
      sprintf(file, "./data1/snow%06d.dat", z / 2500 );
      //sprintf(file, "./data/snow%06d.dat", z );
      fp = fopen(file, "w");
      
      for( i = 0; i <= n; i++ ){
	
	//printf("x[%d] = %f y[%d] = %f\n", i, x[i], i, y[i]);
	fprintf(fp, "%f %f %f\n", x[i], y[i], 0.0 );
	
      }
      
      fclose(fp);
      
    }
    
    // 新たな辺の長さ
    new_length(n,h,phi,l);
    connect(n,l);
    
    
    // 最大の辺と最小の辺を取り出す
    l_max = l[0];
    
    for( i = 0; i <= n; i++ ){
      
      if( l_max < l[i] ){  //配列i番目の数値がmaxよりも大きかったら
	
	l_max = l[i];   //maxに配列i番目の数値を格納
	
      }
      
    }
    
    
    //----------辺の判定----------//
    
    for( i = 1; i <= n; i++ ){
      
      // 2分割
      if( chi[i] == 0.0 && l[i] > lc ){
	
	printf("i1 = %d\n", i);
	
	n = n + 2;
	
	for( j = n + 1; j >= i + 3; j-- ){
	  
	  h[j] = h[j - 2];
	  l[j] = l[j - 2];
	  t1[j] = t1[j - 2];
	  t2[j] = t2[j - 2]; 
	  
	}
	
	if( t1[i] < 0 && t2[i] > 0){
	  
	  if( t1[i - 1] < 0 && t2[i - 1] < 0){
	    
	    t1[i + 1] = 0.0;
	    t2[i + 1] = 1.0;
	    t1[i + 2] = t1[i];
	    t2[i + 2] = t2[i];
	    h[i + 2] = h[i] + eps;
	    h[i] = h[i] - eps;
	    
	  }
	  
	  else{
	    
	    t1[i + 1] = -sqrt(3) / 2.0; 
	    t2[i + 1] = -1.0 / 2.0; 
	    t1[i + 2] = t1[i];
	    t2[i + 2] = t2[i];
	    h[i + 2] = h[i] - eps;
	    h[i] = h[i] + eps;
	    
	  }
	  
	} 
	
	else if( t1[i] < 0 && t2[i] < 0 ){
	  
	  if( t1[i - 1] < 0 && t2[i - 1] > 0){
	    
	    t1[i + 1] = 0.0;
	    t2[i + 1] = -1.0;
	    t1[i + 2] = t1[i];
	    t2[i + 2] = t2[i];
	    h[i + 2] = h[i] - eps;
	    h[i] = h[i] + eps;
	    
	  }
	  
	  else{
	    
	    t1[i + 1] = -sqrt(3) / 2.0; 
	    t2[i + 1] = 1.0 / 2.0; 
	    t1[i + 2] = t1[i];
	    t2[i + 2] = t2[i];
	    h[i + 2] = h[i] + eps;
	    h[i] = h[i] - eps;
	    
	  }
	  
	}
	
	else if( t1[i] == 0.0  && t2[i] == -1.0 ){
	  
	  if( t1[i - 1] > 0 && t2[i - 1] < 0){
	    
	    t1[i + 1] = -sqrt(3) / 2.0; 
	    t2[i + 1] = -1.0 / 2.0; 
	    t1[i + 2] = t1[i];
	    t2[i + 2] = t2[i];
	    h[i + 2] = h[i] + eps;
	    h[i] = h[i] - eps;
	    
	  }
      
	  else{
	    
	    t1[i + 1] = sqrt(3) / 2.0; 
	    t2[i + 1] = -1.0 / 2.0;
	    t1[i + 2] = t1[i];
	    t2[i + 2] = t2[i];
	    h[i + 2] = h[i] - eps;
	    h[i] = h[i] + eps;
	    
	  }
	  
	}
	
	else if( t1[i] > 0 && t2[i] < 0 ){
	  
	  if( t1[i - 1] > 0 && t2[i - 1] > 0){
	    
	    t1[i + 1] = 0.0;
	    t2[i + 1] = -1.0; 
	    t1[i + 2] = t1[i];
	    t2[i + 2] = t2[i];
	    h[i + 2] = h[i] + eps;
	    h[i] = h[i] - eps;
	    
	  }
	  
	  else{
	    
	    t1[i + 1] = sqrt(3) / 2.0; 
	    t2[i + 1] = 1.0 / 2.0;
	    t1[i + 2] = t1[i];
	    t2[i + 2] = t2[i];
	    h[i + 2] = h[i] - eps;
	    h[i] = h[i] + eps;
	    
	  }
	  
	}
	
	else if( t1[i] > 0 && t2[i] > 0){
	  
	  if( t1[i - 1] > 0 && t2[i - 1] < 0){
	    
	    t1[i + 1] = 0.0;
	    t2[i + 1] = 1.0; 
	    t1[i + 2] = t1[i];
	    t2[i + 2] = t2[i];
	    h[i + 2] = h[i] - eps;
	    h[i] = h[i] + eps;
	    
	  }
	  
	  else{
	    
	    t1[i + 1] = sqrt(3) / 2.0; 
	    t2[i + 1] = -1.0 / 2.0;
	    t1[i + 2] = t1[i];
	    t2[i + 2] = t2[i];
	    h[i + 2] = h[i] + eps;
	    h[i] = h[i] - eps;
	    
	  }
	  
	}
	
	else{
	  
	  if( t1[i - 1] < 0 && t2[i - 1] > 0){
	    
	    t1[i + 1] = sqrt(3) / 2.0; 
	    t2[i + 1] = 1.0 / 2.0;
	    t1[i + 2] = t1[i];
	    t2[i + 2] = t2[i];
	    h[i + 2] = h[i] + eps;
	    h[i] = h[i] - eps;
	    
	  }
	  
	  else{
	    
	    t1[i + 1] = -sqrt(3) / 2.0; 
	    t2[i + 1] = 1.0 / 2.0;
	    t1[i + 2] = t1[i];
	    t2[i + 2] = t2[i];
	    h[i + 2] = h[i] - eps;
	    h[i] = h[i] + eps;
	    
	  }
	  
	}
	
	connect(n,h);
	connect(n,l);
	connect(n,t1);
	connect(n,t2);
	
	
	// 外向き法線ベクトル
	normal(n,t1,t2,n1,n2);
	connect(n,n1);
	connect(n,n2);
	
	
	//----------短い辺を作ろう----------//
	
	x_temp = ( x[i] + x[i - 1] ) / 2.0;
	y_temp = ( y[i] + y[i - 1] ) / 2.0;
	
	h[i + 1] = x_temp * n1[i + 1] + y_temp * n2[i + 1];
	connect(n,h);
	
	
	//接線角度の準備
	tangent_angle_pre(n,t1,t2,D);
	connect(n,D);
	
	// 接線角度
	tangent_angle(n,t1,t2,D,nu);
	connect(n,nu);
	
	// 外角
	outside_angle(n,nu,phi);
	connect(n,phi);
	
	//遷移数の準備
	transition_number_pre(n,n1,n2,sigma);
	connect(n,sigma);
	
	//遷移数
	transition_number(n,sigma,chi);
	connect(n,chi);
	
	// 新たな座標
	new_coordinate(n,h,t1,t2,phi,x,y);
	connect(n,x);
	connect(n,y);
	
	// 新たな辺の長さ
	new_length(n,h,phi,l);
	connect(n,l);
	
      }
      
      // 3分割
      else if( ( chi[i] == 1 || chi[i] == -1 ) && l[i] > lc ){
	
	printf("i2 = %d\n", i);
	
	n = n + 4;
	
	for( j = n + 1; j >= i + 5; j-- ){
	  
	  h[j] = h[j - 4];
	  l[j] = l[j - 4];
	  t1[j] = t1[j - 4];
	  t2[j] = t2[j - 4]; 
	  
	}
	
	
	if( t1[i] < 0 && t2[i] > 0){
	  
	  if( t1[i - 1] < 0 && t2[i - 1] < 0 ){
	    
	    t1[i + 1] = 0.0;
	    t2[i + 1] = 1.0;
	    t1[i + 2] = t1[i];
	    t2[i + 2] = t2[i];
	    t1[i + 3] = -sqrt(3) / 2.0; 
	    t2[i + 3] = -1.0 / 2.0;
	    t1[i + 4] = t1[i];
	    t2[i + 4] = t2[i];
	    h[i + 4] = h[i] - eps;
	    h[i + 2] = h[i] + eps;
	    h[i] = h[i] - eps;
	    
	  }
	  
	  else{
	    
	    t1[i + 1] = -sqrt(3) / 2.0; 
	    t2[i + 1] = -1.0 / 2.0; 
	    t1[i + 2] = t1[i];
	    t2[i + 2] = t2[i];
	    t1[i + 3] = 0.0;
	    t2[i + 3] = 1.0;
	    t1[i + 4] = t1[i];
	    t2[i + 4] = t2[i];
	    h[i + 4] = h[i] + eps;
	    h[i + 2] = h[i] - eps;
	    h[i] = h[i] + eps;
	    
	  }
	  
	} 
	
	else if( t1[i] < 0 && t2[i] < 0 ){
	  
	  if( t1[i - 1] < 0 && t2[i - 1] > 0 ){
	    
	    t1[i + 1] = 0.0;
	    t2[i + 1] = -1.0;
	    t1[i + 2] = t1[i];
	    t2[i + 2] = t2[i];
	    t1[i + 3] = -sqrt(3) / 2.0; 
	    t2[i + 3] = 1.0 / 2.0;
	    t1[i + 4] = t1[i];
	    t2[i + 4] = t2[i];
	    h[i + 4] = h[i] + eps;
	    h[i + 2] = h[i] - eps;
	    h[i] = h[i] + eps;
	    
	  }
	  
	  else{
	    
	    t1[i + 1] = -sqrt(3) / 2.0; 
	    t2[i + 1] = 1.0 / 2.0; 
	    t1[i + 2] = t1[i];
	    t2[i + 2] = t2[i];
	    t1[i + 3] = 0.0;
	    t2[i + 3] = -1.0;
	    t1[i + 4] = t1[i];
	    t2[i + 4] = t2[i];
	    h[i + 4] = h[i] - eps;
	    h[i + 2] = h[i] + eps;
	    h[i] = h[i] - eps;
	    
	  }
	  
	}
	
	else if( t1[i] == 0.0  && t2[i] == -1.0 ){
	  
	  if( t1[i - 1] < 0 && t2[i - 1] < 0 ){
	    
	    t1[i + 1] = sqrt(3) / 2.0;
	    t2[i + 1] = -1.0 / 2.0;
	    t1[i + 2] = t1[i];
	    t2[i + 2] = t2[i];
	    t1[i + 3] = -sqrt(3) / 2.0; 
	    t2[i + 3] = -1.0 / 2.0;
	    t1[i + 4] = t1[i];
	    t2[i + 4] = t2[i];
	    h[i + 4] = h[i] + eps;
	    h[i + 2] = h[i] - eps;
	    h[i] = h[i] + eps;
	    
	  }
	  
	  else{
	    
	    t1[i + 1] = -sqrt(3) / 2.0;
	    t2[i + 1] = -1.0 / 2.0;
	    t1[i + 2] = t1[i];
	    t2[i + 2] = t2[i];
	    t1[i + 3] = sqrt(3) / 2.0; 
	    t2[i + 3] = -1.0 / 2.0;
	    t1[i + 4] = t1[i];
	    t2[i + 4] = t2[i];
	    h[i + 4] = h[i] - eps;
	    h[i + 2] = h[i] + eps;
	    h[i] = h[i] - eps;
	    
	  }
	  
	}
	
	else if( t1[i] > 0 && t2[i] < 0 ){
	  
	  if( t1[i - 1] > 0 && t2[i - 1] > 0 ){
	    
	    t1[i + 1] = 0.0;
	    t2[i + 1] = -1.0;
	    t1[i + 2] = t1[i];
	    t2[i + 2] = t2[i];
	    t1[i + 3] = sqrt(3) / 2.0; 
	    t2[i + 3] = 1.0 / 2.0;
	    t1[i + 4] = t1[i];
	    t2[i + 4] = t2[i];
	    h[i + 4] = h[i] - eps;
	    h[i + 2] = h[i] + eps;
	    h[i] = h[i] - eps;
	    
	  }
	  
	  else{
	    
	    t1[i + 1] = sqrt(3) / 2.0; 
	    t2[i + 1] = 1.0 / 2.0; 
	    t1[i + 2] = t1[i];
	    t2[i + 2] = t2[i];
	    t1[i + 3] = 0.0;
	    t2[i + 3] = -1.0;
	    t1[i + 4] = t1[i];
	    t2[i + 4] = t2[i];
	    h[i + 4] = h[i] + eps;
	    h[i + 2] = h[i] - eps;
	    h[i] = h[i] + eps;
	    
	  }
	  
	}
	
	else if( t1[i] > 0 && t2[i] > 0){
	  
	  if( t1[i - 1] > 0 && t2[i - 1] < 0 ){
	    
	    t1[i + 1] = 0.0;
	    t2[i + 1] = 1.0;
	    t1[i + 2] = t1[i];
	    t2[i + 2] = t2[i];
	    t1[i + 3] = sqrt(3) / 2.0; 
	    t2[i + 3] = -1.0 / 2.0;
	    t1[i + 4] = t1[i];
	    t2[i + 4] = t2[i];
	    h[i + 4] = h[i] + eps;
	    h[i + 2] = h[i] - eps;
	    h[i] = h[i] + eps;
	    
	  }
	  
	  else{
	    
	    t1[i + 1] = sqrt(3) / 2.0; 
	    t2[i + 1] = -1.0 / 2.0; 
	    t1[i + 2] = t1[i];
	    t2[i + 2] = t2[i];
	    t1[i + 3] = 0.0;
	    t2[i + 3] = 1.0;
	    t1[i + 4] = t1[i];
	    t2[i + 4] = t2[i];
	    h[i + 4] = h[i] - eps;
	    h[i + 2] = h[i] + eps;
	    h[i] = h[i] - eps;
	    
	  }
	  
	}
	
	else{
	  
	  if( t1[i - 1] > 0 && t2[i - 1] > 0 ){
	    
	    t1[i + 1] = -sqrt(3) / 2.0;
	    t2[i + 1] = 1.0 / 2.0;
	    t1[i + 2] = t1[i];
	    t2[i + 2] = t2[i];
	    t1[i + 3] = sqrt(3) / 2.0; 
	    t2[i + 3] = 1.0 / 2.0;
	    t1[i + 4] = t1[i];
	    t2[i + 4] = t2[i];
	    h[i + 4] = h[i] + eps;
	    h[i + 2] = h[i] - eps;
	    h[i] = h[i] + eps;
	    
	  }
	  
	  else{
	    
	    t1[i + 1] = sqrt(3) / 2.0;
	    t2[i + 1] = 1.0 / 2.0;
	    t1[i + 2] = t1[i];
	    t2[i + 2] = t2[i];
	    t1[i + 3] = -sqrt(3) / 2.0; 
	    t2[i + 3] = 1.0 / 2.0;
	    t1[i + 4] = t1[i];
	    t2[i + 4] = t2[i];
	    h[i + 4] = h[i] - eps;
	    h[i + 2] = h[i] + eps;
	    h[i] = h[i] - eps;
	    
	  }
	  
	}
	
	connect(n,h);
	connect(n,l);
	connect(n,t1);
	connect(n,t2);
	
	
	// 外向き法線ベクトル
	normal(n,t1,t2,n1,n2);
	connect(n,n1);
	connect(n,n2);
	
	
	//----------短い辺を作ろう---------//
	
	// 1:1:1の分割
	x_temp1 = ( x[i] + 2 * x[i - 1] ) / 3.0;
	y_temp1 = ( y[i] + 2 * y[i - 1] ) / 3.0;
	x_temp2 = ( 2 * x[i] + x[i - 1] ) / 3.0;
	y_temp2 = ( 2 * y[i] + y[i - 1] ) / 3.0;
	
	/*
	// 1:sqrt(2):1の分割
	x_temp1 = ( x[i] + ( 1 + sqrt(2) ) * x[i - 1] ) / ( 2 + sqrt(2) );
	y_temp1 = ( y[i] + ( 1 + sqrt(2) ) * y[i - 1] ) / ( 2 + sqrt(2) );
	x_temp2 = ( ( 1 + sqrt(2) ) * x[i] + x[i - 1] ) / ( 2 + sqrt(2) );
	y_temp2 = ( ( 1 + sqrt(2) ) * y[i] + y[i - 1] ) / ( 2 + sqrt(2) );
	*/
	
	h[i + 1] = x_temp1 * n1[i + 1] + y_temp1 * n2[i + 1];
	h[i + 3] = x_temp2 * n1[i + 3] + y_temp2 * n2[i + 3];
	connect(n,h);
	
	
	//接線角度の準備
	tangent_angle_pre(n,t1,t2,D);
	connect(n,D);
	
	// 接線角度
	tangent_angle(n,t1,t2,D,nu);
	connect(n,nu);
	
	// 外角
	outside_angle(n,nu,phi);
	connect(n,phi);
	
	//遷移数の準備
	transition_number_pre(n,n1,n2,sigma);
	connect(n,sigma);
	
	//遷移数
	transition_number(n,sigma,chi);
	connect(n,chi);
	
	// 新たな座標
	new_coordinate(n,h,t1,t2,phi,x,y);
	connect(n,x);
	connect(n,y);
	
	// 新たな辺の長さ
	new_length(n,h,phi,l);
	connect(n,l);
	
      }
      
      
      //辺を減らす
      else if( chi[i] == 0.0 && ( l[i] / l_max ) < 1.0e-3 ){
	
	printf("i3 = %d\n", i);
	
	h[i - 1] = ( h[i - 1] + h[i + 1] ) / 2.0;
	
	
	for( j = i; j <= n; j++ ){
	  
	  h[j] = h[j + 2];
	  t1[j] = t1[j + 2];
	  t2[j] = t2[j + 2];
	  
	}
	
	n = n - 2;
	
	connect(n,h);
	connect(n,t1);
	connect(n,t2);
	
	
	// 外向き法線ベクトル
	normal(n,t1,t2,n1,n2);
	connect(n,n1);
	connect(n,n2);
	
	//接線角度の準備
	tangent_angle_pre(n,t1,t2,D);
	connect(n,D);
	
	// 接線角度
	tangent_angle(n,t1,t2,D,nu);
	connect(n,nu);
	
	// 外角
	outside_angle(n,nu,phi);
	connect(n,phi);
	
	//遷移数の準備
	transition_number_pre(n,n1,n2,sigma);
	connect(n,sigma);
	
	//遷移数
	transition_number(n,sigma,chi);
	connect(n,chi);
	
	// 新たな座標
	new_coordinate(n,h,t1,t2,phi,x,y);
	connect(n,x);
	connect(n,y);
	
	// 新たな辺の長さ
	new_length(n,h,phi,l);
	connect(n,l);
	
      }
      
    }
    
   
    // 新たな辺の長さ
    new_length(n,h,phi,l);
    connect(n,l);
    
    // 曲率
    for( i = 1; i <= n; i++ ){
      
      kappa[i] = alpha * chi[i] / l[i];
      
    }
    
    connect(n,kappa);
    
    printf("t = %d n = %d\n", z, n);
    
  }
  
  
  //--------------------領域の解放--------------------//
  
  free(x);
  free(y);
  free(l);
  free(t1);
  free(t2);
  free(n1);
  free(n2);
  free(h);
  free(nu);
  free(phi);
  free(D);
  free(sigma);
  free(chi);
  free(kappa);

  for( i = 0; i <= N; i++ ){
    
    free((void *)T[i]);
    
  }
  
  free((void *)T);
  
  free(p);
  free(q);
  free(B);
  free(d);
  
  for( i = 0; i <= s; i++ ){
    
    free((void *)u[i]);
    
  }
  
  free((void *)u);
  
  
  //----------return----------//
  
  return 0;
  
}


//--------------------関数--------------------//

void connect(int n, double *a){
  
  a[0] = a[n];
  a[n + 1] = a[1];
  
}

void normal(int n, double *t1, double *t2, double *n1, double *n2){
  
  int i;
  
  for( i = 1; i <= n; i++ ){
    
    n1[i] = t2[i];
    n2[i] = -t1[i];
    
  }
  
}

void height(int n, double *x, double *y, double *n1, double *n2, double *h){
  
  int i;
  
  for( i = 1; i <= n; i++ ){
    
    h[i] = x[i] * n1[i] + y[i] * n2[i];
    
  }
  
}

void tangent_angle_pre(int n, double *t1, double *t2, double *D){
  
  int i;
  
  for( i = 1; i <= n; i++ ){
    
    if( t1[i] * t2[i + 1] - t2[i] * t1[i + 1] < 0 ){
      
      D[i] = -1;
      
    }
    
    else if( t1[i] * t2[i + 1] - t2[i] * t1[i + 1] > 0 ){
      
      D[i] = 1;
      
    }
    
    else{
      
      D[i] = 0;
      
    }
    
  }

}

void tangent_angle(int n, double *t1, double *t2, double *D, double *nu){
  
  int i;
  
  if( t2[1] < 0 ){
    
    nu[1] = -acos(t1[1]);
    
  }
  
  else{
    
    nu[1] = acos(t1[1]);
    
  }
  
  for( i = 1; i <= n; i++ ){
    
    nu[i + 1] = nu[i] + D[i] * acos( t1[i] * t1[i + 1] + t2[i] * t2[i + 1] );
    
  }
  
  nu[0] = nu[1] - ( nu[n + 1] - nu[n] );
  
}

void outside_angle(int n, double *nu, double *phi){
  
  int i;
  
  for( i = 1; i <= n; i++ ){
    
    phi[i] = nu[i + 1] - nu[i];
    
  }
  
}

void transition_number_pre(int n, double *n1, double *n2, double *sigma){
  
  int i;
  
  for( i = 1; i <= n; i++ ){
    
    if( n1[i - 1] * n2[i] - n2[i - 1] * n1[i] < 0 ){
      
      sigma[i] = -1;
      
    }
    
    else if( n1[i - 1] * n2[i] - n2[i - 1] * n1[i] > 0 ){
      
      sigma[i] = 1;
      
    }
    
    else{
      
      sigma[i] = 0;
      
    }
    
  }
  
}

void transition_number(int n, double *sigma, double *chi){
  
  int i;
  
  for( i = 1; i <= n; i++ ){
    
    chi[i] = ( sigma[i] + sigma[i + 1] ) / 2.0;
    
  }
  
}

void new_coordinate(int n, double *h, double *t1, double *t2, double *phi, double *x, double *y){
  
  int i;
  
  for( i = 1; i <= n; i++ ){
    
    x[i] = ( h[i + 1] * t1[i] - h[i] * t1[i + 1] ) / sin(phi[i]);
    y[i] = ( h[i + 1] * t2[i] - h[i] * t2[i + 1] ) / sin(phi[i]);
    
  }
  
}

void new_length(int n, double *h, double *phi, double *l){
  
  int i;
  
  for( i = 1; i <= n; i++ ){
    
    l[i] = ( h[i + 1] / sin(phi[i]) ) - ( h[i] * ( ( 1.0 / tan(phi[i]) ) + ( 1.0 / tan(phi[i - 1] ) ) ) ) + ( h[i - 1] / sin(phi[i - 1]) );
    
  }
  
}

double ff( int i, double *p, double t ){
  
  return ( -p[i] );
  
}

double gg( int i, int j, double *l, double *x, double *y, double *t1, double *t2, double a ){

  return ( (
     
	    ( l[j] + ( ( x[j - 1] - x[i - 1] - a * t1[i] ) * t1[j] + ( y[j - 1] - y[i - 1] - a * t2[i] ) * t2[j] ) )
	    * log( sqrt( ( l[j] + ( ( x[j - 1] - x[i - 1] - a * t1[i] ) * t1[j] + ( y[j - 1] - y[i - 1] - a * t2[i] ) * t2[j] ) ) * ( l[j] + ( ( x[j - 1] - x[i - 1] - a * t1[i] ) * t1[j] + ( y[j - 1] - y[i - 1] - a * t2[i] ) * t2[j] ) ) + ( ( ( x[j - 1] - x[i - 1] - a * t1[i] ) * t2[j] - ( y[j - 1] - y[i - 1] - a * t2[i] ) * t1[j] ) ) * ( ( ( x[j - 1] - x[i - 1] - a * t1[i] ) * t2[j] - ( y[j - 1] - y[i - 1] - a * t2[i] ) * t1[j] ) ) + epsl * epsl ) )
	    - l[j]
	    + fabs( ( ( x[j - 1] - x[i - 1] - a * t1[i] ) * t2[j] - ( y[j - 1] - y[i - 1] - a * t2[i] ) * t1[j] ) ) * atan( ( l[j] + ( ( x[j - 1] - x[i - 1] - a * t1[i] ) * t1[j] + ( y[j - 1] - y[i - 1] - a * t2[i] ) * t2[j] ) ) / sqrt( ( ( ( x[j - 1] - x[i - 1] - a * t1[i] ) * t2[j] - ( y[j - 1] - y[i - 1] - a * t2[i] ) * t1[j] ) * ( ( x[j - 1] - x[i - 1] - a * t1[i] ) * t2[j] - ( y[j - 1] - y[i - 1] - a * t2[i] ) * t1[j] ) + epsl * epsl ) ) )
	    - ( ( x[j - 1] - x[i - 1] - a * t1[i] ) * t1[j] + ( y[j - 1] - y[i - 1] - a * t2[i] ) * t2[j] ) * log( sqrt( ( ( ( x[j - 1] - x[i - 1] - a * t1[i] ) * t1[j] + ( y[j - 1] - y[i - 1] - a * t2[i] ) * t2[j] ) ) * ( ( ( x[j - 1] - x[i - 1] - a * t1[i] ) * t1[j] + ( y[j - 1] - y[i - 1] - a * t2[i] ) * t2[j] ) ) + ( ( ( x[j - 1] - x[i - 1] - a * t1[i] ) * t2[j] - ( y[j - 1] - y[i - 1] - a * t2[i] ) * t1[j] ) ) * ( ( ( x[j - 1] - x[i - 1] - a * t1[i] ) * t2[j] - ( y[j - 1] - y[i - 1] - a * t2[i] ) * t1[j] ) ) + epsl * epsl ) )
	    - fabs( ( ( x[j - 1] - x[i - 1] - a * t1[i] ) * t2[j] - ( y[j - 1] - y[i - 1] - a * t2[i] ) * t1[j] ) ) * atan( ( ( x[j - 1] - x[i - 1] - a * t1[i] ) * t1[j] + ( y[j - 1] - y[i - 1] - a * t2[i] ) * t2[j] ) / sqrt( ( ( ( x[j - 1] - x[i - 1] - a * t1[i] ) * t2[j] - ( y[j - 1] - y[i - 1] - a * t2[i] ) * t1[j] ) * ( ( x[j - 1] - x[i - 1] - a * t1[i] ) * t2[j] - ( y[j - 1] - y[i - 1] - a * t2[i] ) * t1[j] ) + epsl * epsl ) ) )
	    ) );

}
    

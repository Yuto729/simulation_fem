#ifndef __APPLE__
#include <GL/glut.h>
#else
#include <GLUT/glut.h>
#endif

//Mac (retina display)の場合は1の部分を2に変更する必要がある場合も
#define RETINA_SCALE 1

#include "MathTool.h"
#include "Mesh.h"
#include "SolidFEM.h"
#include "GLTool.h"

//ウィンドウのサイズ設定
#define WINDOW_SIZE_X	1024
#define WINDOW_SIZE_Y	768
//マウスによる回転角のゲイン
#define ANGLE_GAIN	0.2
//マウスのズームゲイン
#define ZOOM_GAIN	0.001

//インタフェースの状態を表す列挙型定数
enum{
	IM_NODE_SELECTION,
	IM_OBJECT_INTERACTIVE,
	IM_FEM_ANALYSIS
};

//パラメータ設定
#define POISSONRATIO	0.49//[arb. unit]
#define YOUNGMODULUS	10//[MPa]
#define MAX_MS			10.0//[N/mm^2]

//必要なグローバル変数を宣言
//物体
Mesh mesh;
//座標変換行列
Matd matrix_scale;
Matd matrix_rotate;
Matd matrix_object2world;
Matd matrix_world2object;
Matd matrix4opengl;
//ウィンドウの幅
int width;
//ウィンドウの高さ
int height;
//一フレーム前のマウス座標
Vec3d mouse_pre;
//ドラッグ開始前クリック点
Vec3d pos_clicked;
//ドラッグ中前クリック点
Vec3d pos_dragging;
//入力インタフェースのモード
int interface_mode;
//カメラのズーム量
double camera_zoom;
//マススのボタン
int mouse_button;
//ミーゼス応力を取得する座標
Vec3d pos;
//特定の座標のミーゼス応力
double ms_at;


//ディスプレイ・コールバック関数
void glDisplay( void )
{
	glClearColor( 1.0, 1.0, 1.0, 1.0 );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	//カメラの設定を行う
	setCamera( width, height );

	//物体を描画する
	glPushMatrix();
	glMultMatrixd( matrix4opengl.X );
	//renderGrid( 5 );
	renderFEMMesh( &mesh, MAX_MS );
	glPopMatrix();

	glutSwapBuffers();
}

//変換行列の更新
void updateMatrix()
{
	//回転とスケール行列から変換行列を計算
	multiMatandMat( &matrix_rotate, &matrix_scale, &matrix_object2world);
	//逆行列の計算（逆変換）
	invMat( &matrix_object2world, &matrix_world2object );
	//OpenGL用に転置行列を用意
	trMat( &matrix_object2world, &matrix4opengl );
}


void glIdle( void )
{
	glutPostRedisplay();
}

void glResize( int _width, int _height )
{
	width = _width;
	height = _height;
	setCamera( width, height );
}


//マウスボタン・コールバック関数
void glMouse( int _button, int _state, int _x, int _y )
{
	Vec3d pos_window;
	Vec3d pos_world;
	Vec3d pos_object;
	_x = RETINA_SCALE * _x;
	_y = RETINA_SCALE * _y;

	if ( _state == GLUT_DOWN ) {
		mouse_button = _button;
		pos_window.x = _x;
		pos_window.y = _y;
		//クリックした点のZバッファを再取得
		pos_window.z = getDepth( _x, _y );
		convertWindow2World( &pos_window, &pos_world );
		multiTransferMatandVec3( &matrix_world2object, &pos_world, &pos_object );
		if( interface_mode == IM_NODE_SELECTION){//ノード選択モードの場合
			switch ( _button ) {
				case GLUT_LEFT_BUTTON://左ボタンの場合
					selectMeshNodeNear( &mesh, &pos_object, NODE_FIXED, &pos_object );
					break;
				case GLUT_RIGHT_BUTTON://右ボタンの場合
					if( selectMeshNodeNear( &mesh, &pos_object, NODE_DEFORM, &pos_object ) != -1 ){
						//荷重領域を設定
						setDeformRegion( &mesh );
						memcpy( pos_clicked.X, pos_object.X, sizeof( double ) * 3 );
						memcpy( pos_dragging.X, pos_object.X, sizeof( double ) * 3 );
					}
					break;
				default:
					break;
			}
		}
		setVec3(&mouse_pre, _x, _y, 0);
	}else{
		mouse_button = -1;
	}
}

//マウス動作・コールバック関数
void glMotion( int _x, int _y)
{
	double output;
	Vec3d input;
	Vec3d pos_window;
	Vec3d pos_world;
	Vec3d pos_mouse;
	switch( interface_mode ){
		case IM_OBJECT_INTERACTIVE:
			//マウス動作による座標変換行列の設定
			if( mouse_button == GLUT_LEFT_BUTTON ){
				//左ボタンの場合は回転
				setVec3(&pos_mouse, (_x - mouse_pre.x)*ANGLE_GAIN, (_y - mouse_pre.y)*ANGLE_GAIN, 0);
				setMouseRotation( pos_mouse.x, pos_mouse.y , &matrix_rotate );
			}else if( mouse_button == GLUT_RIGHT_BUTTON ){
				//右ボタンの場合はズーム
				camera_zoom += ZOOM_GAIN * ((_x - mouse_pre.x) + (mouse_pre.y - _y));
				setMouseScroll( camera_zoom, &matrix_scale );
			}
			updateMatrix();
			break;
		case IM_FEM_ANALYSIS:
			//ドラッグ開始点のZバッファを再取得
			convertWorld2Window( &pos_clicked, &pos_window );
			//ドラッグ中の点の３次元世界座標を取得
			pos_window.x = _x;
			pos_window.y = _y;
			convertWindow2World( &pos_window, &pos_world );
			//入力はオブジェクト座標系で与える
			multiTransferMatandVec3( &matrix_world2object, &pos_world, &pos_dragging );
			//ドラッグ開始点と現在のマウスの座標の相対ベクトルを入力とする
			subVec3andVec3( &pos_dragging, &pos_clicked, &input );

			//[TODO5ヒント]
			//inputに変位ベクトルを設定すれば指定した変位を入力できる

			//境界条件が設定されていればマウスを動かすたびに解析を実行
			if( mesh.is_boundary_on == 1){
				//境界条件を設定
				setDeformCondition( &mesh, &input );
				//剛性方程式を解く
				solveStiffnessEquation( &mesh );
				//最大ミーゼス応力を計算
				output = calTotalMisessStress( &mesh );

				//[TODO5ヒント]
				//指定した座標のミーゼス応力を取得
				setVec3( &pos, 0, 0, 0 );
				ms_at = getMisessStressAt( &mesh, pos );
				//結果をコンソールに表示
				printf( "Max MS = %3.3f, MS at (%2.1f, %2.1f, %2.1f)= %3.3f\r", output, pos.x, pos.y, pos.z, ms_at);
			}
			break;
		default:
			break;
	}
	setVec3(&mouse_pre, _x, _y, 0);
}


//キーボード・コールバック関数
void glKeyboard( unsigned char _key, int _x, int _y )
{
	switch (_key) {//押されたキーを判定
		case 'q':
		case 'Q':
		case '\033':
			exit(0);
			break;
		case 'n'://ノード選択モード
			interface_mode = IM_NODE_SELECTION;
			break;
		case 'r'://物体回転・スケールモード
			interface_mode = IM_OBJECT_INTERACTIVE;
			break;
		case 'f'://有限要素解析モード
			interface_mode = IM_FEM_ANALYSIS;
			break;
		case 's'://前処理
			//拘束領域を設定
			setFixRegion( &mesh );
			//逆行列を計算する
			printf("Pre matrix calculation ....");
			calPreMatrix( &mesh );
			printf("[OK]\n");
			//荷重領域を設定
			setDeformRegion( &mesh );
			break;
		case 'i'://変位をクリア
			clearDeform( &mesh );
			break;
		case 'c'://ノード選択をクリア
			clearMeshNodeSelected( &mesh );
			break;
		default:
			break;
	}
}

int main( int _argc, char *_argv[] )
{
	char filename[256];
	char title[256];
	double current_time, last_time;

	//メッシュデータのロード
	if(_argc>1)//コマンドライン入力の場合
		sprintf( filename, "%s", _argv[1] );
	else//デフォルトパス
		sprintf( filename, "../data/armadillo.fem" );

	//メッシュの初期化
	initMesh( &mesh );
	//メッシュの読み込み
	printf("Loading mesh ...");
	if(loadFEMMesh( &mesh, filename ) == -1 ){
		printf( "[FAIL]\n");
		return -1;
	}
	printf("[OK]\n");

	//物性パラメータの設定
	setMaterialProperty( &mesh, POISSONRATIO, YOUNGMODULUS );
	//剛性行列の生成
	setTotalStiffnessMatrix( &mesh );

	//[TODO5ヒント]
	//固定境界を設定するにはmesh.node[ #ID ].state = NODE_FIXED;
	//変位境界を設定するにはmesh.node[ #ID ].state = NODE_DEFORM;

	//初期ウィンドウサイズの初期化
	width = WINDOW_SIZE_X;
	height = WINDOW_SIZE_Y;	
	//インタフェースモードの初期設定
	interface_mode = IM_OBJECT_INTERACTIVE;
	//カメラのズームを初期化
	camera_zoom = 1;
	//マウスボタンの初期化
	mouse_button = -1;


	//座標変換行列の初期化
	initMat( &matrix_scale );
	initMat( &matrix_rotate );
	initMat( &matrix4opengl );
	initMat( &matrix_object2world );
	initMat( &matrix_world2object );
	setMatDim( &matrix_scale, 4, 4 );
	setMatDim( &matrix_rotate, 4, 4 );
	setMatDim( &matrix4opengl, 4, 4 );
	setMatDim( &matrix_object2world, 4, 4 );
	setMatDim( &matrix_world2object, 4, 4 );
	identityMat( &matrix_scale );
	identityMat( &matrix_rotate );
	identityMat( &matrix4opengl );
	identityMat( &matrix_object2world );
	identityMat( &matrix_world2object );
	clearVec3( &mouse_pre );
	
	//キーボード・マウスの操作方法表示
	printf("[r] to change the view (L:rotate, R:scale)\n");
	printf("[n] to select/unselect a node as (L:fixed, R:deform)\n");
	printf("[s] to start matrix calculation\n");
	printf("[f] to input deformation and execute FEM analysis\n");
	printf("[i] to initialize deformation\n");
	printf("[c] to clear all boundary conditions\n");
	printf("[Esc/q] to terminate the program\n");

	//ウィンドウ生成と初期化
	glutInit( &_argc, _argv );
	glutInitDisplayMode( GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH );
	glutInitWindowSize( WINDOW_SIZE_X, WINDOW_SIZE_Y );
	glutCreateWindow( _argv[ 0 ] );

	//OpenGLコールバック関数
	glutDisplayFunc( glDisplay );
	glutReshapeFunc( glResize );
	glutMouseFunc( glMouse );
	glutMotionFunc( glMotion );
	glutKeyboardFunc( glKeyboard );
	glutIdleFunc( glIdle );
	glInit();
	glutMainLoop();

	//メッシュの開放
	releaseMesh( &mesh );

	//座標変換行列の開放
	releaseMat( &matrix_scale );
	releaseMat( &matrix_rotate );
	releaseMat( &matrix4opengl );
	releaseMat( &matrix_object2world );
	releaseMat( &matrix_world2object );

	return 0;
}
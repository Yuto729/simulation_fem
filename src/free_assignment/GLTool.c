#include "GLTool.h"

void calColorMap(double _value, Vec3d *_color)
{
	int Hi;
	double H, f;
	if( _value > 1 ){
		_value = 1;
	}
	H = 6 * 0.7 * ( 1 - _value );
	Hi = ( int )floor( H ) % 6;
	f = H - ( double ) Hi;
	switch( Hi ){
		case 0:
			_color->x = 1; _color->y = f; _color->z = 0;
			break;
		case 1:
			_color->x = 1 - f; _color->y = 1; _color->z = 0;
			break;
		case 2:
			_color->x = 0; _color->y = 1; _color->z = f;
			break;
		case 3:
			_color->x = 0; _color->y = 1 - f; _color->z = 1;
			break;
		case 4:
			_color->x = f; _color->y = 0; _color->z = 1;
			break;
		case 5:
			_color->x = 1; _color->y = 0; _color->z = 1 - f;
			break;
	}
}

void renderFEMMesh( Mesh *_mesh, double _max_mises_stress )
{
	Vec3d color;
	unsigned int i,j,k;
	//[TODO3]描画処理の書き写し
	//ノードの描画
	glPointSize(10);
	glBegin(GL_POINTS);
	for (i = 0; i < _mesh->num_node; i++) {
		//頂点の状態に応じて色を変える
		switch (_mesh->node[i].state) {
			case NODE_FREE:
				glColor3d(0, 1, 1);
				break;
			case NODE_FIXED:
				glColor3d(1, 0, 1);
				break;
			case NODE_DEFORM:
				glColor3d(1, 1, 0);
				break;
		}
		//頂点の座標を指定
		glVertex3dv(_mesh->node[i].new_position.X);
	}
	glEnd();

	//ノード間ラインの描画(変形前・変形後)
	for (i = 0; i < _mesh->num_tetrahedra; i++) {
		for (j = 0; j < 4; j++) {
			for (k = 0; k < j; k++) {
				glLineWidth(0.5);
				glColor3d(0.3, 0.3, 0.3);
				glBegin(GL_LINE_STRIP);
				glVertex3dv(_mesh->tetrahedra[i].position[j].X);
				glVertex3dv(_mesh->tetrahedra[i].position[k].X);
				glEnd();
				glLineWidth(1);
				glColor3d(0, 0, 0);
				glBegin(GL_LINE_STRIP);
				glVertex3dv(_mesh->tetrahedra[i].new_position[j].X);
				glVertex3dv(_mesh->tetrahedra[i].new_position[k].X);
				glEnd();
			}
		}
	}
	//要素ポリゴンの描画
	for (i = 0; i < _mesh->num_tetrahedra; i++) {
		if (_mesh->tetrahedra[i].status == ELEMENT_FAILED) {
			// 破壊した要素は黒色で表示
			color.x = 0.1;
			color.y = 0.1;
			color.z = 0.1;
			
			// 破壊した要素はワイヤーフレーム表示
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		} else {
			// 通常の要素はミーゼス応力に応じた色で表示
			calColorMap(_mesh->tetrahedra[i].mises_stress / _max_mises_stress, &color);
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		}
		for (j = 0; j < 4; j++) {
			glColor3dv(color.X);
			glBegin(GL_TRIANGLES);
			glVertex3dv(_mesh->tetrahedra[i].new_position[(j + 0) % 4].X);
			glVertex3dv(_mesh->tetrahedra[i].new_position[(j + 1) % 4].X);
			glVertex3dv(_mesh->tetrahedra[i].new_position[(j + 2) % 4].X);
			glEnd();
		}
	}
	// 描画モードを元に戻す
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

float getDepth( int _pos_window_x, int _pos_window_y )
{
	float depth;
	GLint viewport[ 4 ];
	glGetIntegerv( GL_VIEWPORT, viewport );
	//デプスバッファの取得
	glReadPixels( _pos_window_x, viewport[3]-_pos_window_y, 1, 1,
				GL_DEPTH_COMPONENT, GL_FLOAT, &depth );
	return depth;
}

void convertWorld2Window( Vec3d *_position_world, Vec3d *_position_window )
{
    GLdouble matrix_modelview[ 16 ];
	GLdouble matrix_projection[ 16 ];
	GLint viewport[ 4 ];
	glGetIntegerv( GL_VIEWPORT,viewport );
    glGetDoublev( GL_MODELVIEW_MATRIX, matrix_modelview );
    glGetDoublev( GL_PROJECTION_MATRIX, matrix_projection );
	//ワールド座標系からウィンドウ座標系へ変換
    gluProject( _position_world->x, _position_world->y, _position_world->z,
				matrix_modelview, matrix_projection, viewport,
				&_position_window->x, &_position_window->y, &_position_window->z );
}

void convertWindow2World( Vec3d *_position_window, Vec3d *_position_world)
{
	GLdouble matrix_modelview[ 16 ];
	GLdouble matrix_projection[ 16 ];
	GLint viewport[ 4 ];
    glGetIntegerv( GL_VIEWPORT, viewport );
    glGetDoublev( GL_MODELVIEW_MATRIX, matrix_modelview );
    glGetDoublev( GL_PROJECTION_MATRIX, matrix_projection );
	//ウィンドウ座���系からワールド座標系へ変換
	gluUnProject( _position_window->x, ( double )viewport[ 3 ]-_position_window->y, _position_window->z,
				matrix_modelview, matrix_projection, viewport,
				&_position_world->x, &_position_world->y, &_position_world->z );
}

void glInit( void )
{
	glEnable( GL_DEPTH_TEST );
	glEnable( GL_LINE_SMOOTH );
	glShadeModel( GL_FLAT );
	glDisable( GL_CULL_FACE );
	glCullFace( GL_FRONT );
	glFrontFace( GL_CW );
	glEnable( GL_BLEND );
	glClearColor( 1.0, 1.0, 1.0, 1.0 );
}

void setCamera( int _width, int _height )
{
	glViewport( 0, 0, _width, _height );
	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();
	glFrustum( -0.02, 0.02, -0.02*(double)_height / _width, 0.02*(double)_height / _width, 0.1, 1000 );
	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();
	gluLookAt( 0.0, 0.0, 150.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0 );
}

void renderGrid( double _scale )
{
	int i;
	glColor3d( 0, 0, 0 );
	for( i = -10; i <= 10; i++ ){
		glBegin( GL_LINE_STRIP );
		glVertex3d( -10 * _scale, 0, i * _scale );
		glVertex3d( 10 * _scale, 0, i * _scale );
		glEnd();
		glBegin( GL_LINE_STRIP );
		glVertex3d( i * _scale, 0, -10 * _scale );
		glVertex3d( i * _scale, 0, 10 * _scale );
		glEnd();
	}
}

void setMouseRotation( double _x, double _y, Matd *_dst )
{
	Matd matrix_rot_x;
	Matd matrix_rot_y;
	Matd matrix_temp;
	initMat( &matrix_rot_x );
	initMat( &matrix_rot_y );
	initMat( &matrix_temp );
	setRotationalMatrix( _x, ROT_AXIS_Y,  &matrix_rot_y);
	setRotationalMatrix( _y, ROT_AXIS_X,  &matrix_rot_x );
	multiMatandMat( &matrix_rot_y, _dst, &matrix_temp );
	multiMatandMat( &matrix_rot_x, &matrix_temp, _dst );
	releaseMat( &matrix_rot_x );
	releaseMat( &matrix_rot_y );
	releaseMat( &matrix_temp );
}

void setMouseScroll( double _s, Matd *_dst)
{
	_dst->X[0] = _dst->X[5] = _dst->X[10] = _s;
}

void renderTetrahedra(Tetrahedra *_tetrahedra) {
    // 四面体の各面を描画
    glBegin(GL_TRIANGLES);
    
    // 面1: 0-1-2
    glVertex3dv(_tetrahedra->new_position[0].X);
    glVertex3dv(_tetrahedra->new_position[1].X);
    glVertex3dv(_tetrahedra->new_position[2].X);
    
    // 面2: 0-1-3
    glVertex3dv(_tetrahedra->new_position[0].X);
    glVertex3dv(_tetrahedra->new_position[1].X);
    glVertex3dv(_tetrahedra->new_position[3].X);
    
    // 面3: 0-2-3
    glVertex3dv(_tetrahedra->new_position[0].X);
    glVertex3dv(_tetrahedra->new_position[2].X);
    glVertex3dv(_tetrahedra->new_position[3].X);
    
    // 面4: 1-2-3
    glVertex3dv(_tetrahedra->new_position[1].X);
    glVertex3dv(_tetrahedra->new_position[2].X);
    glVertex3dv(_tetrahedra->new_position[3].X);
    
    glEnd();
}
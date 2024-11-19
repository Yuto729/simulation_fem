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
	//ウィンドウ座標系からワールド座標系へ変換
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
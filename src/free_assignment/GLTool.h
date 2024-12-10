#ifndef GLTOOL_H
#define GLTOOL_H

#ifndef __APPLE__
#include <GL/glut.h>
#else
#include <GLUT/glut.h>
#endif

#include "MathTool.h"
#include "Mesh.h"

//カラーマップの計算（HからRGBへの変換）
void calColorMap(double _value, Vec3d *_color);
//ボリュームメッシュを描画する
void renderFEMMesh( Mesh *_mesh, double _max_mises_stress );
//ウィンドウ座標の点のデプス値を取得する
float getDepth( int _pos_window_x, int _pos_window_y );
//ワールド座標からウィンドウ座標へ変換する
void convertWorld2Window( Vec3d *_position_world, Vec3d *_position_window );
//ウィンドウ座標からワールド座標へ変換する
void convertWindow2World( Vec3d *_position_window, Vec3d *_position_world);
//OpenGLの初期化
void glInit( void );
//カメラの設定
void setCamera( int _width, int _height );
//格子の描画
void renderGrid( double _scale );
//回転変換行列の設定
void setMouseRotation( double _x, double _y, Matd *_dst );
//拡大縮小変換行列の設定
void setMouseScroll( double _s, Matd *_dst);

#endif

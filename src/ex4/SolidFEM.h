#ifndef SOLIDFEM_H
#define SOLIDFEM_H

#include "MathTool.h"
#include "Mesh.h"

//要素の物性値の設定
void setMaterialProperty( Mesh *_mesh, double _poisson_ratio, double _young_modulus );
//形状関数の取得
void setGeometoryMatrix( Tetrahedra *_tetrahedra, Matd *_A);
//要素ひずみの計算
void calStrain( Tetrahedra *_tetrahedra );
//要素応力の計算
void calStress( Tetrahedra *_tetrahedra );
//要素ミーゼス応力の計算
double calMisesStress( Tetrahedra *_tetrahedra );
//要素体積の計算
double calVolume( Tetrahedra *_tetrahedra );
//ひずみ変位関係式の行列[Be]の設定
void setStrainDeformationMatrix( Tetrahedra *_tetrahedra );
//応力ひずみ関係式の行列[De]の設定
void setStressStrainMatrix( Tetrahedra *_tetrahedra );
//要素剛性行列[Ke]の設定
void setStiffnessMatrix( Tetrahedra *_tetrahedra );

//全体剛性行列[K]の設定
void setTotalStiffnessMatrix( Mesh *_mesh );
//ディリクレ固定条件領域の設定
void setFixRegion( Mesh *_mesh );
//逆行列の前計算
void calPreMatrix( Mesh *_mesh );
//ディリクレ変位条件領域の設定
void setDeformRegion( Mesh *_mesh );
//ディリクレ変位条件を与える
void setDeformCondition( Mesh *_mesh, Vec3d *_deformation );
//剛性方程式を解く（全ての点の変位が求まる）
void solveStiffnessEquation( Mesh *_mesh );
//全体のミーゼス応力を計算する (戻り値は最大値)
double calTotalMisessStress( Mesh *_mesh );
//指定した座標のミーゼス応力を取得する
double getMisessStressAt( Mesh *_mesh, Vec3d _position);
//全ノードの変位と力をファイルに出力する
int saveDF( Mesh *_mesh, const char *_filename );
//変位を初期化する
void clearDeform( Mesh *_mesh );

#endif
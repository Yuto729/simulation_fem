#ifndef MESH_H
#define MESH_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "MathTool.h"

#define TOLERANCE 1E-12

enum{
	NODE_FREE,
	NODE_FIXED,
	NODE_DEFORM
};

typedef struct{
	Vec3d position;
	Vec3d new_position;
	int state;
}Node;

typedef struct{
	unsigned int node_index[4];
	Vec3d position[4];
	Vec3d new_position[4];
	VecNd deformation;
	int material;
	Matd K;
	Matd D;
	Matd B;
	Matd dN;
	double poisson_ratio;
	double young_modulus;
	double volume;
	double mises_stress;
	VecNd strain;
	VecNd stress;
}Tetrahedra;

typedef struct{
	unsigned int num_node;
	unsigned int num_tetrahedra;
	Node *node;
	Tetrahedra *tetrahedra;

	unsigned int num_S;
	unsigned int *S;
	unsigned int num_Sd;
	unsigned int *Sd;
	unsigned int num_Sn;
	unsigned int *Sn;
	int is_boundary_on;

	VecNd force;
	VecNd Fd;
	VecNd Fn;

	VecNd deformation;
	VecNd Ud;
	VecNd Un;

	Matd K;
	Matd Ks;
	Matd Ls;
	Matd Ldd;
	Matd Lnd;
}Mesh;

//要素の初期化
void initTetrahedra( Tetrahedra *_tetrahedra );
//要素の開放
void releaseTetrahedra( Tetrahedra *_tetrahedra );
//点の内外判定(初期形状を利用)
int isPointInside( Tetrahedra *_tetrahedra, Vec3d _position );

//メッシュの初期化
void initMesh( Mesh *_mesh );
//メッシュの節点，要素のメモリ確保を行う
void setMeshDim( Mesh *_mesh );
//メッシュの節点，要素のメモリを開放する
void releaseMesh( Mesh *_mesh );
//メッシュの節点，要素データを読み込む
int loadFEMMesh( Mesh *_mesh, const char *_filename );
//メッシュの節点を選択
int selectMeshNodeNear( Mesh *_mesh, Vec3d *_pos, int _state, Vec3d *_selected );
//メッシュの節点を全解除
void clearMeshNodeSelected( Mesh *_mesh );

#endif
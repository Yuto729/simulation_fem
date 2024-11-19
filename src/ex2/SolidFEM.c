#include "SolidFEM.h"

void setMaterialProperty( Mesh *_mesh, double _poisson_ratio, double _young_modulus )
{
	unsigned int i;
	//全ての要素についてポアッソン比とヤング率を設定
	for( i = 0; i < _mesh->num_tetrahedra; i ++ ){
		_mesh->tetrahedra[ i ].poisson_ratio = _poisson_ratio;
		_mesh->tetrahedra[ i ].young_modulus = _young_modulus;
	}
}

void setGeometoryMatrix( Tetrahedra *_tetrahedra, Matd *_A)
{
	unsigned int i,j;
	initMat( _A );
	setMatDim( _A, 4, 4 );
	for( i = 0; i < 4; i ++){
		for( j = 0; j < 4; j ++){
			if( j == 0 ){
				_A->X[ 4 * i + j ] = 1;
			}
			else{
				_A->X[ 4 * i + j ] = _tetrahedra->position[ i ].X[ j - 1 ];
			}
		}
	}
}

void calStrain( Tetrahedra *_tetrahedra )
{
	//[TODO2]_tetrahedra->Bと_tetrahedra->deformationから_tetrahedra->strainを計算する
	//ヒント：行列とベクトルの積には関数multiMatandVecNを使用する
}

void calStress( Tetrahedra *_tetrahedra )
{
	calStrain( _tetrahedra );
	//[TODO2]_tetrahedra->Dと_tetrahedra->strainから_tetrahedra->stressを計算する
}

double calMisesStress( Tetrahedra *_tetrahedra )
{
	calStress( _tetrahedra );
	//[TODO5]_tetrahedra->mises_stressにミーゼス応力の計算結果を格納する
	return _tetrahedra->mises_stress;
}

double calVolume( Tetrahedra *_tetrahedra)
{
	Matd A;
	initMat( &A );
	setGeometoryMatrix( _tetrahedra, &A );
	_tetrahedra->volume = fabs( detMat( &A ) ) / 6.0;
	releaseMat( &A );
	return _tetrahedra->volume;
}

void setStrainDeformationMatrix( Tetrahedra *_tetrahedra )
{
	unsigned int i,j;
	Matd A;
	Matd invA;
	initMat( &A );
	initMat( &invA );
	calVolume( _tetrahedra );
	setGeometoryMatrix( _tetrahedra, &A );
	invMat( &A, &invA );

	//[TODO2]invAから_tetrahedra->dNを設定する

	//[TODO2]_tetrahedra->dNから_tetrahedra->Bを設定する

	releaseMat( &A );
	releaseMat( &invA );
}

void setStressStrainMatrix( Tetrahedra *_tetrahedra )
{
	unsigned int i;
	double Dscale;
	Dscale=_tetrahedra->young_modulus / 
		( ( 1.0 + _tetrahedra->poisson_ratio ) * ( 1.0 - 2.0 * _tetrahedra->poisson_ratio ) );

	//[TODO2]行列_tetrahedra->Dを設定する
}

void setStiffnessMatrix( Tetrahedra *_tetrahedra )
{
	//V, B, Dを使って要素のKを設定する
	Matd trB;
	Matd temp1;
	Matd temp2;
	initMat( &trB );
	initMat( &temp1 );
	initMat( &temp2 );
	trMat( &_tetrahedra->B, &trB );

	//[TODO2]_tetrahedra->B, _tetrahedra->volume, trB, _tetrahedra->Dを用いて_tetrahedra->Kを設定する

	releaseMat( &trB );
	releaseMat( &temp1 );
	releaseMat( &temp2 );
}

void setTotalStiffnessMatrix( Mesh *_mesh )
{
	unsigned int i, j, k, l, m;
	unsigned int col, row;
	clearMat(& _mesh->K );
	for( i = 0; i<_mesh->num_tetrahedra; i ++ ){
		setStrainDeformationMatrix( &_mesh->tetrahedra[ i ] );
		setStressStrainMatrix( &_mesh->tetrahedra[ i ] );
		setStiffnessMatrix( &_mesh->tetrahedra[ i ] );
	
		//[TODO4]要素行列_mesh->tetrahedra[ i ].Kを足し込み，全体剛性行列_mesh->Kを生成する
		
	}
}

void setFixRegion( Mesh *_mesh )
{
	unsigned int i;
	_mesh->num_S = 0;
	for( i = 0; i < _mesh->num_node; i ++ ){
		if( _mesh->node[ i ].state != NODE_FIXED ){//ディリクレ拘束以外
			_mesh->S[ _mesh->num_S ] = i;//集合は全ノードの中から取得
			_mesh->num_S ++;
		}
	}
	if(_mesh->num_S > 0){
		_mesh->is_boundary_on = 1;
	}
}

void calPreMatrix( Mesh *_mesh )
{
	unsigned int i, j, k;
	int row, col;
	if(	_mesh->is_boundary_on == 1 ){
		//行列の次元を設定
		setMatDim( &_mesh->Ks,  _mesh->num_S * 3, _mesh->num_S * 3 );
		setMatDim( &_mesh->Ls,  _mesh->num_S * 3, _mesh->num_S * 3 );
		for( i = 0; i < _mesh->num_S ; i ++ ){//並び替えた後の順序
			for( j = 0; j < _mesh->num_S; j ++ ){//並び替えた後の順序
				for( k = 0; k < 3; k++ ){//自由度
					row= 3 * _mesh->S[ i ] + k ;
					col= 3 * _mesh->S[ j ];
					memcpy( &_mesh->Ks.X[ _mesh->num_S * 3 *( 3 * i + k ) + 3 * j ],
							&_mesh->K.X[ _mesh->num_node * 3 * row + col ],
							sizeof( double ) * 3 );
				}
			}
		}
		//高速化のため,予め逆行列を計算
		invMat( &_mesh->Ks, &_mesh->Ls );
	}
}

void setDeformRegion( Mesh *_mesh )
{
	unsigned int i, j, k;
	unsigned int col, row;
	_mesh->num_Sd = 0;
	_mesh->num_Sn = 0;
	for( i = 0; i < _mesh->num_S; i ++ ){
		switch( _mesh->node[ _mesh->S[ i ] ].state ){
			case NODE_DEFORM://ディリクレ変位条件
				_mesh->Sd[ _mesh->num_Sd ] = i;//集合はSの中から取得
				_mesh->num_Sd ++;
				break;
			case NODE_FREE://ノイマン条件
				_mesh->Sn[ _mesh->num_Sn ] = i;//集合はSの中から取得
				_mesh->num_Sn ++;
				break;
		}
	}
	//行列・ベクトルの次元が決まる
	if( _mesh->num_Sd > 0 && _mesh->num_Sn > 0 ){
		setMatDim( &_mesh->Ldd, _mesh->num_Sd * 3, _mesh->num_Sd * 3 );
		setMatDim( &_mesh->Lnd, _mesh->num_Sd * 3, _mesh->num_Sn * 3 );
		setVecNDim( &_mesh->Ud, _mesh->num_Sd * 3);
		setVecNDim( &_mesh->Un, _mesh->num_Sn * 3);
		setVecNDim( &_mesh->Fd, _mesh->num_Sd * 3);
		setVecNDim( &_mesh->Fn, _mesh->num_Sn * 3);

		for( i = 0; i < _mesh->num_Sd; i ++ ){//並び替えた後の順序
			for( j = 0; j < _mesh->num_Sd; j ++ ){//並び替えた後の順序
				for( k = 0; k < 3; k++ ){//自由度
					col = 3 * _mesh->Sd[ j ];
					row = 3 * _mesh->Sd[ i ] + k;
					memcpy( &_mesh->Ldd.X[ _mesh->num_Sd * 3 *( 3 * i + k ) + 3 * j ],
							&_mesh->Ls.X[ _mesh->num_S * 3 * row + col],
							sizeof( double ) * 3 );
				}
			}
		}
		for( i = 0; i < _mesh->num_Sn; i ++ ){//並び替えた後の順序
			for( j = 0; j < _mesh->num_Sd; j ++ ){//並び替えた後の順序
				for( k = 0; k < 3; k++ ){//自由度
					col = 3 * _mesh->Sd[ j ];
					row = 3 * _mesh->Sn[ i ] + k;
					memcpy( &_mesh->Lnd.X[ _mesh->num_Sd * 3 *( 3 * i + k ) + 3 * j ],
							&_mesh->Ls.X[ _mesh->num_S * 3 * row + col ],
							sizeof( double ) * 3 );
				}
			}
		}
	}
}

void setDeformCondition( Mesh *_mesh, Vec3d *_deformation )
{
	unsigned int i;
	clearVecN( &_mesh->Ud );
	for( i = 0; i < _mesh->num_Sd; i ++ ){
		if( _mesh->node[ _mesh->S[ _mesh->Sd[ i ] ] ].state == NODE_DEFORM){
			memcpy( &_mesh->Ud.X[ 3 * i ] ,_deformation->X, sizeof(double) * 3);
		}
	}
}

void solveStiffnessEquation( Mesh *_mesh )
{
	unsigned int i,j;
	Matd invLdd;
	if( _mesh->num_Sd > 0 ){
		initMat( &invLdd );
		invMat( &_mesh->Ldd, &invLdd );
		multiMatandVecN( &invLdd, &_mesh->Ud, &_mesh->Fd);
		multiMatandVecN( &_mesh->Lnd, &_mesh->Fd, &_mesh->Un );
		//境界領域集合に基づいて全体ベクトルに戻す
		for( i = 0; i < _mesh->num_Sd; i ++ ){
			memcpy( &_mesh->deformation.X[ 3 * _mesh->S[ _mesh->Sd[ i ] ] ],
					&_mesh->Ud.X[ 3 * i ], sizeof(double) * 3 );
		}
		for( i = 0; i < _mesh->num_Sn; i ++ ){
			memcpy( &_mesh->deformation.X[ 3 * _mesh->S[ _mesh->Sn[ i ] ] ],
					&_mesh->Un.X[ 3 * i ], sizeof(double) * 3 );
		}
		for( i = 0; i < _mesh->num_node; i ++ ){
			for( j = 0; j < 3; j ++ ){
				_mesh->node[ i ].new_position.X[j] = _mesh->node[ i ].position.X[j]
													+ _mesh->deformation.X[ 3 * i + j ];
			}
		}
		//要素の変位にも反映する
		for( i = 0; i < _mesh->num_tetrahedra; i ++ ){
			for( j = 0; j < 4; j ++ ){
				memcpy( &_mesh->tetrahedra[ i ].deformation.X[ 3 * j ],
					&_mesh->deformation.X[ 3 * _mesh->tetrahedra[ i ].node_index[ j ] ],
					sizeof( double ) * 3 );
				memcpy( _mesh->tetrahedra[ i ].new_position[ j ].X,
					_mesh->node[ _mesh->tetrahedra[ i ].node_index[ j ] ].new_position.X,
					sizeof( double ) * 3 );
			}	
		}	
		releaseMat( &invLdd );
	}
}

double calTotalMisessStress( Mesh *_mesh )
{
	unsigned int i,j;
	double max_mises_stress = 0;
	for( i = 0; i < _mesh->num_tetrahedra; i ++ ){
		//ミーゼス応力の計算
		calMisesStress( &_mesh->tetrahedra[ i ] );
		//値の大小判定
		if( max_mises_stress < _mesh->tetrahedra[ i ].mises_stress ){
			max_mises_stress = _mesh->tetrahedra[ i ].mises_stress;
		}
	}
	return max_mises_stress;
}

double getMisessStressAt( Mesh *_mesh, Vec3d _position )
{
	unsigned int i;
	double ms = 0;
	for( i = 0; i < _mesh->num_tetrahedra; i ++ ){
		if (isPointInside( &_mesh->tetrahedra[ i ], _position)) {
			ms = _mesh->tetrahedra[ i ].mises_stress;
			return ms;
		}
	}
	return ms;
}

int saveDF( Mesh *_mesh, const char *_filename )
{
	FILE *file;
	unsigned int i;

	if( (file = fopen( _filename, "w" ) ) == NULL || _mesh->is_boundary_on != 1){
		return -1;
	}
	fprintf( file, "index,x[mm],y[mm],z[mm],Ux[mm],Uy[mm],Uz[mm],Fx[N],Fy[N],Fz[N]\n");

	//[TODO6]
	//全体変位ベクトル_mesh->deformationと全体剛性行列_mesh->Kから全体力ベクトル_mesh->forceを計算
	//ファイルfileに節点番号，節点3次元座標，節点3次元変位，節点3次元力を出力する

	fclose( file );
	return 1;
}

void clearDeform( Mesh *_mesh )
{
	unsigned int i;
	clearVecN( &_mesh->deformation );
	clearVecN( &_mesh->force );
	for( i = 0; i <_mesh->num_node; i ++ ){
		memcpy(	&_mesh->node[ i ].new_position,
			&_mesh->node[ i ].position,
			sizeof( Vec3d ) );
	}
	for( i = 0; i <_mesh->num_tetrahedra; i ++ ){
		clearVecN( &_mesh->tetrahedra[i].deformation );
		memcpy( _mesh->tetrahedra[ i ].new_position,
				_mesh->tetrahedra[ i ].position,
				sizeof( Vec3d ) * 4 );
		_mesh->tetrahedra[ i ].mises_stress=0;
	}
}
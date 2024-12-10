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
	multiMatandVecN(&_tetrahedra->B, &_tetrahedra->deformation, &_tetrahedra->strain);
}

void calStress( Tetrahedra *_tetrahedra )
{
	calStrain( _tetrahedra );
	//[TODO2]_tetrahedra->Dと_tetrahedra->strainから_tetrahedra->stressを計算する
	multiMatandVecN(&_tetrahedra->D, &_tetrahedra->strain, &_tetrahedra->stress);
}

double calMisesStress( Tetrahedra *_tetrahedra )
{
	calStress( _tetrahedra );
	//[TODO5]_tetrahedra->mises_stressにミーゼス応力の計算結果を格納する
	double sigma_x = _tetrahedra->stress.X[0]; // σxx
	double sigma_y = _tetrahedra->stress.X[1]; // σyy
	double sigma_z = _tetrahedra->stress.X[2]; // σzz
	double tau_xy = _tetrahedra->stress.X[3]; // τxy
	double tau_xz = _tetrahedra->stress.X[4]; // τyz
	double tau_yz = _tetrahedra->stress.X[5]; // τzx
	
	_tetrahedra->mises_stress = sqrt(0.5 * ((sigma_x-sigma_y)*(sigma_x-sigma_y) + (sigma_y-sigma_z)*(sigma_y-sigma_z) + (sigma_z-sigma_x)*(sigma_z-sigma_x) + 
									6.0 * (tau_xy*tau_xy + tau_xz*tau_xz + tau_yz*tau_yz)));
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
	unsigned int i;
	Matd A;
	Matd invA;
	initMat( &A );
	initMat( &invA );
	calVolume( _tetrahedra );
	setGeometoryMatrix( _tetrahedra, &A );
	invMat( &A, &invA );
	//[TODO2]invAから_tetrahedra->dNを設定する

	//[TODO2]_tetrahedra->dNから_tetrahedra->Bを設定する

	// dNの設定 (3×4行列)
	setMatDim(&_tetrahedra->dN, 4, 3);
	clearMat(&_tetrahedra->dN);
	
	for(i = 0; i < 4; i++) {
		// 形状関数の偏微分を設定
		_tetrahedra->dN.X[4*0 + i] = invA.X[4*1 + i];  // ∂N/∂x
		_tetrahedra->dN.X[4*1 + i] = invA.X[4*2 + i];  // ∂N/∂y
		_tetrahedra->dN.X[4*2 + i] = invA.X[4*3 + i];  // ∂N/∂z
	}

	// B行列の設定 (6×12行列)
	setMatDim(&_tetrahedra->B, 12, 6);
	clearMat(&_tetrahedra->B);
	
	for(i = 0; i < 4; i++) {
		// εxx = ∂u/∂x
		_tetrahedra->B.X[12*0 + 3*i + 0] = _tetrahedra->dN.X[4*0 + i];
		
		// εyy = ∂v/∂y
		_tetrahedra->B.X[12*1 + 3*i + 1] = _tetrahedra->dN.X[4*1 + i];
		
		// εzz = ∂w/∂z
		_tetrahedra->B.X[12*2 + 3*i + 2] = _tetrahedra->dN.X[4*2 + i];
		
		// γxy = ∂u/∂y + ∂v/∂x
		_tetrahedra->B.X[12*3 + 3*i + 0] = _tetrahedra->dN.X[4*1 + i];  // ∂u/∂y
		_tetrahedra->B.X[12*3 + 3*i + 1] = _tetrahedra->dN.X[4*0 + i];  // ∂v/∂x
		
		// γyz = ∂v/∂z + ∂w/∂y
		_tetrahedra->B.X[12*4 + 3*i + 1] = _tetrahedra->dN.X[4*2 + i];  // ∂v/∂z
		_tetrahedra->B.X[12*4 + 3*i + 2] = _tetrahedra->dN.X[4*1 + i];  // ∂w/∂y
		
		// γzx = ∂w/∂x + ��u/∂z
		_tetrahedra->B.X[12*5 + 3*i + 2] = _tetrahedra->dN.X[4*0 + i];  // ∂w/∂x
		_tetrahedra->B.X[12*5 + 3*i + 0] = _tetrahedra->dN.X[4*2 + i];  // ∂u/∂z
	}

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
	setMatDim(&_tetrahedra->D, 6, 6);
	clearMat(&_tetrahedra->D);
	
	// 対角成分
	// 対角成分
	_tetrahedra->D.X[0] = Dscale * (1.0 - _tetrahedra->poisson_ratio);
	_tetrahedra->D.X[7] = Dscale * (1.0 - _tetrahedra->poisson_ratio);
	_tetrahedra->D.X[14] = Dscale * (1.0 - _tetrahedra->poisson_ratio);
	_tetrahedra->D.X[21] = Dscale * (1.0 - 2.0*_tetrahedra->poisson_ratio) / 2.0;
	_tetrahedra->D.X[28] = Dscale * (1.0 - 2.0*_tetrahedra->poisson_ratio) / 2.0;
	_tetrahedra->D.X[35] = Dscale * (1.0 - 2.0*_tetrahedra->poisson_ratio) / 2.0;
	
	// 非対角成分
	_tetrahedra->D.X[1] = Dscale * _tetrahedra->poisson_ratio;
	_tetrahedra->D.X[2] = Dscale * _tetrahedra->poisson_ratio;
	_tetrahedra->D.X[6] = Dscale * _tetrahedra->poisson_ratio;
	_tetrahedra->D.X[8] = Dscale * _tetrahedra->poisson_ratio;
	_tetrahedra->D.X[12] = Dscale * _tetrahedra->poisson_ratio;
	_tetrahedra->D.X[13] = Dscale * _tetrahedra->poisson_ratio;
}

void setStiffnessMatrix( Tetrahedra *_tetrahedra )
{
	Matd trB, temp1, temp2;
	initMat(&trB);
	initMat(&temp1);
	initMat(&temp2);

	// B^T * D * B * V の順で計算
	trMat(&_tetrahedra->B, &trB);
	multiMatandMat(&trB, &_tetrahedra->D, &temp1);
	multiMatandMat(&temp1, &_tetrahedra->B, &temp2);
	scalingMat(_tetrahedra->volume, &temp2, &_tetrahedra->K);

	releaseMat(&trB);
	releaseMat(&temp1);
	releaseMat(&temp2);
}

void setTotalStiffnessMatrix( Mesh *_mesh )
// TODO4
{
    unsigned int i, j, k, l, m;
    unsigned int col, row;
    clearMat(&_mesh->K); 

    for (i = 0; i < _mesh->num_tetrahedra; i++) {
        setStrainDeformationMatrix(&_mesh->tetrahedra[i]);
        setStressStrainMatrix(&_mesh->tetrahedra[i]);
        setStiffnessMatrix(&_mesh->tetrahedra[i]);

        for (j = 0; j < 4; j++) { 
            for (k = 0; k < 4; k++) {
                for (l = 0; l < 3; l++) { 
                    for (m = 0; m < 3; m++) {
                        row = 3 * _mesh->tetrahedra[i].node_index[j] + l;
                        col = 3 * _mesh->tetrahedra[i].node_index[k] + m;

                        _mesh->K.X[_mesh->K.ncol * row + col] +=
                            _mesh->tetrahedra[i].K.X[12 * (3 * j + l) + (3 * k + m)];
                    }
                }
            }
        }
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

	multiMatandVecN(&_mesh->K, &_mesh->deformation, &_mesh->force);
	
	for(i = 0; i < _mesh->num_node; i++) {
		fprintf(file, "%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
				i,
				_mesh->node[i].position.X[0],
				_mesh->node[i].position.X[1],
				_mesh->node[i].position.X[2],
				_mesh->deformation.X[3*i],
				_mesh->deformation.X[3*i+1],
				_mesh->deformation.X[3*i+2],
				_mesh->force.X[3*i],
				_mesh->force.X[3*i+1],
				_mesh->force.X[3*i+2]);
	}
	
	fclose(file);
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

void setFailureThreshold(Mesh *_mesh, double _threshold) {
    _mesh->failure_threshold = _threshold;
}

void updateFailureStatus(Mesh *_mesh) {
    int had_new_failure = 0;
    
    // 各要素のミーゼス応力を計算し、破壊判定を行う
    for(unsigned int i = 0; i < _mesh->num_tetrahedra; i++) {
        if(_mesh->tetrahedra[i].status == 0) {  // まだ破壊していない要素のみチェック
            calMisesStress(&_mesh->tetrahedra[i]);
            
            if(_mesh->tetrahedra[i].mises_stress > _mesh->failure_threshold) {
                _mesh->tetrahedra[i].status = 1;  // 破壊状態に設定
                _mesh->tetrahedra[i].young_modulus *= 0.000001;  // 破壊した要素の剛性を100万分の1に低下
                had_new_failure = 1;
                _mesh->num_failed_elements++;
            }
        }
    }
    
    // 破壊が発生した場合、剛性行列を再計算
    if(had_new_failure) {
        setTotalStiffnessMatrix(_mesh);
        calPreMatrix(_mesh);
    }
}

void resetFailureStatus(Mesh *_mesh, double _young_modulus) {
    for(unsigned int i = 0; i < _mesh->num_tetrahedra; i++) {
        _mesh->tetrahedra[i].status = 0;
        _mesh->tetrahedra[i].young_modulus = _young_modulus;  // 元の剛性に戻す
    }
    _mesh->num_failed_elements = 0;
    setTotalStiffnessMatrix(_mesh);
    calPreMatrix(_mesh);
}

int isSimulationComplete(Mesh *_mesh) {
    // 新しい破壊が発生しなくなったら終了
    return (_mesh->num_failed_elements == 0);
}
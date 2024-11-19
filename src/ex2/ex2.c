#include "MathTool.h"
#include "Mesh.h"
#include "SolidFEM.h"

int main ( void )
{
	Tetrahedra tetrahedra;
	initTetrahedra( &tetrahedra );
	tetrahedra.poisson_ratio = 0.49;
	tetrahedra.young_modulus = 10;
	setVec3( &tetrahedra.position[0], 0, 0, 0 );
	setVec3( &tetrahedra.position[1], 1, 0, 0 );
	setVec3( &tetrahedra.position[2], 0, 1, 0 );
	setVec3( &tetrahedra.position[3], 0, 0, 1 );

	setStrainDeformationMatrix( &tetrahedra );
	setStressStrainMatrix( &tetrahedra );	
	setStiffnessMatrix( &tetrahedra );
	printf("V=%f\n", tetrahedra.volume);
	printf("[B]=\n");
	printMat( &tetrahedra.B );
	printf("[D]=\n");
	printMat( &tetrahedra.D );
	printf("[dN]=\n");
	printMat( &tetrahedra.dN );
	printf("[K]=\n");
	printMat( &tetrahedra.K );
	releaseTetrahedra( &tetrahedra );

	return 0;
}

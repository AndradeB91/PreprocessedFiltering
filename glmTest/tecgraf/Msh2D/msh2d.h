/** msh2d.h - header files
* ------------------------------------------------------------------
*/

#ifndef MSH2D_H
#define MSH2D_H

#include "deflib2D.h"

#ifdef __cplusplus
extern "C" {
#endif


/** BIBLIOTECA COMPUTACIONAL PARA GERA��O DE MALHAS BIDIMENSIONAIS
    DE ELEMENTOS FINITOS
   @author Antonio Carlos de O. Miranda
   @version 1.0
*/

/**@name Malhas Estruturadas */
/*@{*/

/** Mapeamento transfinito bilinear para regi�es quadrilaterais com qualquer
    forma para os lados. O n�mero de segmentos em lados opostos devem ser
    iguais (Campos, 1991).
    @param bry       (in) vetor de coordenadas do contorno
    @param m         (in) n�mero de n�s na primeira dire��o
    @param n         (in) n�mero de n�s na outra dire��o
    @param elem_type (in) tipo de elemento:  T3 (3), T6 (6), Q4 (4) ou  Q8 (8)
    @param diagtype  (in) op��o de triangula��o:
                     1 - diagonal orientada para direita,
                     2 - diagonal orientada para esquerda,
                     3 - modo cruzado alternado,
                     4 - melhor diagonal
    @param nno       (out) n�mero de n�s gerados
    @param nel       (out) n�mero de elementos gerados
    @param pt        (out) coordenadas dos n�s da malha
    @param conn      (out) conectividade dos elementos
    @return          1 - sucesso, 0 - falha
*/
MSH2D_API int Msh2DBilinear (double *bry, int m, int n, int elem_type, int diagtype,
									  int *nno, int *nel, double **pt, int **conn);
MSH2D_API int Msh2DTryBilinear (double *bry, int np, int elem_type, int diagtype,
										  int *nno, int *nel, double **pt, int **conn);


/** Mapeamento transfinito bilinear para regi�es triangulares que conceitualmente
    � formado por uma regi�o quadrilateral com um dos lados colapsado para um
    ponto.  Os tr�s lados podem ter qualquer forma e os lados adjacentes ao
    lado colapsado devem ter o mesmo n�mero de segmentos (Campos, 1991).
    @param bry       (in) vetor de coordenadas do contorno
    @param m         (in) n�mero de n�s na primeira dire��o
    @param n         (in) n�mero de n�s na outra dire��o
    @param elem_type (in) tipo de elemento:  T3 (3), T6 (6), Q4 (4) ou  Q8 (8)
    @param diagtype  (in) op��o de triangula��o
                     1 - diagonal orientada para direita,
                     2 - diagonal orientada para esquerda,
                     3 - modo cruzado alternado,
                     4 - melhor diagonal
    @param nno       (out) n�mero de n�s gerados
    @param nel       (out) n�mero de elementos gerados
    @param pt        (out) coordenadas dos n�s da malha
    @param conn      (out) conectividade dos elementos
    @return          1 - sucesso, 0 - falha
*/
MSH2D_API int Msh2DCollBilinear (double *bry, int m, int n, int elem_type, int diagtype,
											int *nno, int *nel, double **pt, int **conn);


/** Mapeamento transfinito linear entre dois lados opostos de uma regi�o
    quadrilateral. Os outros dois lados devem ser linhas retas. O n�mero
    de segmentos nos lados opostos devem ser iguais (Campos, 1991).
    @param bry       (in) vetor de coordenadas do contorno
    @param m         (in) n�mero de n�s na primeira dire��o
    @param n         (in) n�mero de n�s na outra dire��o
    @param dir       (in) dire��o do lofting:
                          0 => primeira dire��o, 1 => outra dire��o
    @param weigth    (in) peso aplicav�l a dire�ao do lofting
    @param elem_type (in) tipo de elemento:  T3 (3), T6 (6), Q4 (4) ou  Q8 (8)
    @param diagtype  (in) op��o de triangula��o
                     1 - diagonal orientada para direita,
                     2 - diagonal orientada para esquerda,
                     3 - modo cruzado alternado,
                     4 - melhor diagonal
    @param nno       (out) n�mero de n�s gerados
    @param nel       (out) n�mero de elementos gerados
    @param pt        (out) coordenadas dos n�s da malha
    @param conn      (out) conectividade dos elementos
    @return          1 - sucesso, 0 - falha
*/
MSH2D_API int Msh2DLoft (double *bry, int m, int n, int dir, double weight,
								 int elem_type, int diagtype,
								 int *nno, int *nel, double **pt, int **conn);


/** Mapeamento transfinito linear entre um ponto e uma curva, que forma uma
    regi�o triangular. Os dois lados da regi�o que s�o adjacentes ao ponto
    de lofting tem de ser linhas retas e ter o mesmo n�mero de pontos
    (Campos, 1991).
    @param bry       (in) vetor de coordenadas do contorno
    @param m         (in) n�mero de n�s na primeira dire��o
    @param n         (in) n�mero de n�s na outra dire��o
    @param dir       (in) dire��o do lofting:
                          0 => primeira dire��o, 1 => outra dire��o
    @param weigth    (in) peso aplicav�l a dire�ao do lofting
    @param elem_type (in) tipo de elemento:  T3 (3), T6 (6), Q4 (4) ou  Q8 (8)
    @param diagtype  (in) op��o de triangula��o:
                     1 - diagonal orientada para direita,
                     2 - diagonal orientada para esquerda,
                     3 - modo cruzado alternado,
                     4 - melhor diagonal
    @param nno       (out) n�mero de n�s gerados
    @param nel       (out) n�mero de elementos gerados
    @param pt        (out) coordenadas dos n�s da malha
    @param conn      (out) conectividade dos elementos
    @return          1 - sucesso, 0 - falha
*/
MSH2D_API int Msh2DCollLoft (double *bry, int m, int n, double weight,
									  int elem_type, int diagtype,
									  int *nno, int *nel, double **pt, int **conn);


/** Mapeamento transfinito trilinear para uma regi�o triangular com qualquer
    forma para os lados.  O n�mero de segmentos em todos os lados devem ser
    iguais (Campos, 1991).
    @param bry       (in) vetor de coordenadas do contorno
    @param m         (in) n�mero de n�s em um lado do contorno
    @param elem_type (in) tipo de elemento:  T3 (3), T6 (6), Q4 (4) ou  Q8 (8)
    @param diagtype  (in) op��o de triangula��o:
                     1 - diagonal orientada para direita,
                     2 - diagonal orientada para esquerda,
                     3 - modo cruzado alternado,
                     4 - melhor diagonal
    @param nno       (out) n�mero de n�s gerados
    @param nel       (out) n�mero de elementos gerados
    @param pt        (out) coordenadas dos n�s da malha
    @param conn      (out) conectividade dos elementos
    @return          1 - sucesso, 0 - falha
*/
MSH2D_API int Msh2DTrilinear (double *bry, int m, int elem_type,
										int *nno, int *nel, double **pt, int **conn);
MSH2D_API int Msh2DTryTrilinear (double *bry, int np, int elem_type,
											int *nno, int *nel, double **pt, int **conn);
/*@}*/


/**@name Malhas N�o-Estruturadas */
/*@{*/

/** Triangula��o gen�rica aplic�vel a contornos fechados sem auto-interse��o.
    Os pontos do circuito externo devem estar orientados em sentido hor�rio
    e circuitos internos (furos) em sentido anti-hor�rio. A triangula��o �
    feita usando uma t�cnica de contra��o de contorno (Vianna, 1992). Os
    pontos interiores da triangula��o podem ser gerados pelo algoritmo e/ou
    fornecidos. Dois procedimentos s�o utilizados para gerar os pontos
    interiores. O primeiro procedimento utiliza-se de uma estrutura quadtree
    para gerar pontos no centro de cada c�lula e depois faz-se uma contra��o
    do contorno. O outro procedimento gera os pontos interiores em uma grade
    regular.
    @param n_loops    (in) n�mero de circuitos
    @param loop_segs  (in) vetor com o n�mero de segmentos em cada  circuito
    @param bdry_pts   (in) vetor de coordenadas do contorno
    @param gen_intpts (in) indicador para gerar pontos internos
                           ( 0 => nao gera, 1 => gera )
    @param n_add_pts  (in) n�mero de n�s adicionais ao interior do modelo
    @param qt_flag    (in) indicador de como gerar os pontos internos
                           ( 0 => grade, 1 => quadtree )
    @param type_mesh  (in) Tipo de elemento: T3 (3) ou T6 (6)
    @param add_coords (in) vetor de coordenadas do pontos adicionais ao
                           interior do modelo
    @param n_node     (out) n�mero de n�s gerados
    @param coords     (out) coordenadas dos n�s da malha
    @param n_elem     (out) n�mero de elementos gerados
    @param conn       (out) conectividade dos elementos
*/
MSH2D_API void Msh2DBoundContraction (int n_loops, int *loop_segs, double *bdry_pts,
												  int gen_intpts, int n_add_pts, int qt_flag,
												  int type_mesh, double *add_coords, int *n_node,
												  double **coords, int *n_elem, int **conn);

/** Triangula��o gen�rica aplic�vel a contornos fechados sem auto-interse��o.
    A triangula��o � feita em dois contornos ao mesmo tempo. Os dois contornos
    devem ter mesma topologia, mas podem ter geometrias diferentes.
    Os pontos do circuito externo devem estar orientados em sentido hor�rio
    e circuitos internos (furos) em sentido anti-hor�rio. A triangula��o �
    feita usando uma t�cnica de contra��o de contorno (Vianna, 1992). Os
    pontos interiores da triangula��o podem ser gerados pelo algoritmo e/ou
    fornecidos. Dois procedimentos s�o utilizados para gerar os pontos
    interiores. O primeiro procedimento utiliza-se de uma estrutura quadtree
    para gerar pontos no centro de cada c�lula e depois faz-se uma contra��o
    do contorno. O outro procedimento gera os pontos interiores em uma grade
    regular.
    @param n_loops          (in) n�mero de circuitos
    @param loop_segs        (in) vetor com o n�mero de segmentos em cada  circuito
    @param bdry_pts         (in) vetor de coordenadas do contorno
    @param warp_bdry_pts    (in) vetor de coordenadas do segundo contorno
    @param gen_intpts       (in) indicador para gerar pontos internos
                                ( 0 => nao gera, 1 => gera )
    @param gen_bdrypts      (in) indicador para gerar pontos no contorno
                                ( 0 => nao gera, 1 => gera )
    @param n_add_pts        (in) n�mero de n�s adicionais ao interior do modelo
    @param qt_flag          (in) indicador de como gerar os pontos internos
                                ( 0 => grade, 1 => quadtree )
    @param type_mesh        (in) Tipo de elemento: T3 (3) ou T6 (6)
    @param add_coords       (in) vetor de coordenadas do pontos adicionais ao
                                interior do modelo
    @param n_node           (out) n�mero de n�s gerados
    @param coords           (out) coordenadas dos n�s da malha
    @param warp_coords      (out) coordenadas dos n�s da segunda malha
    @param n_elem           (out) n�mero de elementos gerados
    @param conn             (out) conectividade dos elementos
*/
MSH2D_API void Msh2DBoundContractionWarp (int n_loops, int *loop_segs, double *bdry_pts, double *warp_bdry_pts,
												  int gen_intpts, int gen_bdrypts,
												  int n_add_pts, int qt_flag,
												  int type_mesh, double *add_coords, int *n_node,
												  double **coords, double **warp_coords, int *n_elem, int **conn);


/** Triangula��o gen�rica aplic�vel a contornos fechados sem auto-interse��o.
    Os pontos do circuito externo devem estar orientados em sentido hor�rio e
    circuitos internos (furos) em sentido anti-hor�rio. Para elementos Q4 e Q8,
    o n�mero de segmentos (lados de elementos) deve ser par em cada circuito.
    O algoritmo utiliza-se de uma estrutura quadtree e constr�i os elementos do
    interior do dom�nio utilizando padr�es de elementos para cada configura��o
    da quadtree; para a faixa do contorno gera os elementos por uma t�cnica de
    contra��o de contorno (Cavalcante Neto, 1994).
    @param n_loops    (in) n�mero de circuitos
    @param loop_segs  (in) vetor com o n�mero de segmentos em cada  circuito
    @param bdry_pts   (in) vetor de coordenadas do contorno
    @param type_mesh  (in) Tipo de elemento: T3 (3), T6 (6), Q4 (4) ou  Q8 (8)
    @param ref_quad   (in) 1 = refine quadtree, 0 == not refine
    @param n_node     (out) n�mero de n�s gerados
    @param coords     (out) coordenadas dos n�s da malha
    @param n_elem     (out) n�mero de elementos gerados
    @param conn       (out) conectividade dos elementos
*/
MSH2D_API int Msh2DQuadBound (int n_loops, int *loop_segs, double *bdry_pts, int type_mesh, int ref_quad,
										int *n_node, double **coords, int *n_elem, int **conn);


/** Triangula��o gen�rica aplic�vel a contornos fechados sem auto-interse��o.
    Os pontos do circuito externo devem estar orientados em sentido hor�rio e
    circuitos internos (furos) em sentido anti-hor�rio. O algoritmos utiliza
    uma estrutura de quadtree para dar uma grada��o de transi��o na gera��o
    dos n�s interiores. Os pontos do dom�nio s�o gerados conjuntamente com o
    avan�o da fronteira e de acordo com o gradiente dado pela �rvore quadtree.
    O algoritmo cont�m adicionalmente um procedimento de melhoria local da
    malha (Miranda, 1999).
    @param n_loops    (in) n�mero de circuitos
    @param loop_segs  (in) vetor com o n�mero de segmentos em cada  circuito
    @param bdry_pts   (in) vetor de coordenadas do contorno
    @param type_mesh  (in) Tipo de elemento: T3 (3) ou T6 (6)
    @param n_node     (out) n�mero de n�s gerados
    @param coords     (out) coordenadas dos n�s da malha
    @param n_elem     (out) n�mero de elementos gerados
    @param conn       (out) conectividade dos elementos
*/
MSH2D_API int  Msh2DShape (int n_loops, int *loop_segs, double *bdry_pts, int type_mesh,
									int *n_node, double **coords, int *n_elem, int **conn);


/** Vers�o do algoritmo Msh2DShape mas considerando restric�es internas
    @param n_pts      (in) n�mero de pontos de entrada
    @param bdry_pts   (in) vetor de coordenadas
    @bound_edge       (in) n�mero de aresta do contorno
    @inter_edge       (in) n�mero de aresta internas (restric�es)
    @edges            (in) arestas (id i, id j)
    @param type_mesh  (in) Tipo de elemento: T3 (3) ou T6 (6)
    @param n_node     (out) n�mero de n�s gerados
    @param coords     (out) coordenadas dos n�s da malha
    @param n_elem     (out) n�mero de elementos gerados
    @param conn       (out) conectividade dos elementos
*/
MSH2D_API int Msh2DEdge (int n_pts, double *bdry_pts, int bound_edge, int inter_edge,
								 int *edges, int type_mesh, int *n_node, double  **coords,
								 int *n_elem, int **Conn);

/** Adiciona parametros para Msh2DShape e Msh2DEdge
    @param _n      (in) n�mero de parametros
    @param _param  (in) vetor de parametros
                        _param[0] = 0, nao gera pontos internos
*/
MSH2D_API void Msh2DEdgeParams (int _n, double *_param);

/** Q-Morph - m�todo indireto de gera��o de malhas n�o estruturadas de
    quadrilaterais (Owen, 1999)
    @param n_loops    (in) n�mero de circuitos
    @param loop_segs  (in) vetor com o n�mero de segmentos em cada  circuito
    @param bdry_pts   (in) vetor de coordenadas do contorno
    @param type_mesh  (in) Tipo de elemento: 4 (Q4) ou 8 (Q8)
    @param n_node     (out) n�mero de n�s gerados
    @param coords     (out) coordenadas dos n�s da malha
    @param n_elem     (out) n�mero de elementos gerados
    @param conn       (out) conectividade dos elementos
*/
MSH2D_API int Msh2DQuadSeam (int n_loops, int *loop_segs, double *bdry_pts, int type_mesh,
									  int *n_node, double **coords, int *n_elem, int **conn);

MSH2D_API int Msh2DQuadSeamEdge (int n_pts, double *bdry_pts, int bound_edge, int *edges,
											int *n_node, double  **coords, int *n_elem, int **Conn);


/** aplicar e obtem o fator de refinamento aplicado as arestas do contorno
    para costru��o da QuadTree. Aplicado as funcoes Msh2DShape, Msh2DEdge e
    Msh2DQuadSeam
    @param factor    (in) fator > 0.1
*/
MSH2D_API void   Msh2DSetRefFactor (double factor);
MSH2D_API double Msh2DGetRefFactor (void);

/* register auxiliary functions */
typedef MSH2D_API void Msh2DSizeElement (double x, double y, double *size);


/* function to obtain the element size in a region. If NULL,
it uses an Quadtree based on boundary edges              */
MSH2D_API void Msh2DRegFunc (Msh2DSizeElement *msh2d_size);



/*
** ----------------------------------------------------------------------
** Begin 2D/Surface Mesh using TEMPLATES
** ----------------------------------------------------------------------
*/

/* Optional Function: MshSurfSetTargetMesh - target support mesh (the same
 *                    function above). Set the support surface. The surface
 *                    is used to project 3D points to support surface
*/

/* Mesh Algorithm
   - Create quadrilateral mesh in 2D or 3D. If the support mesh is given, all
     points is projected to this surface.
*/
MSH2D_API int Msh2DTemplate (
int     n_sides,          /* # of sides (2, 3, 4)                    (in) */
int     subdvision[4],    /* subdvision of each side                 (in) */
int     dim,              /* dimension 2 = 2D or 3 = 3D              (in) */
int     type,             /* 4 = Q4, 8 = Q8                          (in) */
int     smooth,           /* 1- smooth internal points, 0 - does not (in) */
double  *bdry_pts,        /* coordinates of all points               (in) */
int     *n_node,          /* counts # of pts for meshing            (out) */
double  **coords,         /* coordinate array used for meshing      (out) */
int     *n_elem,          /* number of elements generated           (out) */
int     **Conn            /* elem.connectivity list from meshing    (out) */
);

/* Compute only the number of elements that generates by Msh2DTemplate */
MSH2D_API int Msh2DTemplateNumElements (
int     n_sides,          /* # of sides (2, 3, 4)                   (in) */
int     subdvision[4]    /* subdvision of each side                 (in) */
);

/* Set Options to template - return value or -1.0 if fail */
/* position  --- value
**     0  -  (0 or 1) turn on/off manual settings
**     1  -  value of u1
**     2  -  value of v1
**     3  -  value of u2
**     4  -  value of v2
**
**
**                    A
**      *-----------------------------*
**       |\                           /|
**       | \                         / |
**       |  \ (u1,v1)       (u2,v2) /  |
**       |   *---------------------*   |
**       |   |                     |   |
**     N |   |                     |   | N
**       |   |                     |   |
**       |   |                     |   |
**       |   |                     |   |
**       |   |                     |   |
**       | b |                     | b |
**       *---*---------------------*---*
**                     B
**
**       Subdivision in A < B
**
**
**
**                     A
**       *-----------------------------*
**       |\                            |
**       | \                           | b
**       |  \ (u1,v1)                  |
**       |   *-------------------------*
**       |   |                         |
**     B |   |                         |
**       |   |                         |
**       |   |                         | B
**       |   |                         |
**       |   |                         |
**       | b |                         |
**       *---*-------------------------*
**                       A
**
**       Subdivision can be A != B
**
**
**/
MSH2D_API double Msh2DTemplateSetParam (
int position,
double value
);




/*@}*/


/**@name Fun��es Auxiliares */
/*@{*/

/** Libera a mem�ria do vetor de coordenadas dos pontos passada como argumento de sa�da pelos algoritmos.
    @param points (in) vetor de coordenadas dos pontos
*/
MSH2D_API void Msh2DFreeNodes (double *points);


/** Libera a mem�ria do vetor de conectividade dos elementos passada como argumento de sa�da pelos algoritmos.
    @param Conn (in) vetor de conectividade dos elementos
*/
MSH2D_API void Msh2DFreeConn(int *Conn);

/*@}*/



/**@name Refer�ncias Bibliogr�ficas */
/*@{*/

/**
(Campos, 1991)
Campos, J. A. P. - "Gera��o de Malhas de Elementos Finitos Bidimensionais
Baseada em uma Estrutura de Dados Topol�gica", Disserta��o de Mestrado,
Departamento de Engenharia Civil, 1991.

(Cavalcante Neto, 1994)
Cavalcante Neto, J. B. - "Simula��o Auto-Adaptativa Baseada em Enumera��o
Espacial Recursiva de Modelos Bidimensionais de Elementos Finitos",
Disserta��o de Mestrado, Departamento de Engenharia Civil, PUC/Rio,1994.

(Cavalcante Neto, 1998)
Cavalcante Neto, J. B. - "Gera��o de Malha e Estimativa de Erro para Modelos
Tridimensionais de Elementos Finitos com Trincas", Disserta��o de Doutorado,
Departamento de Engenharia Civil, PUC/Rio, 1998.

(Miranda, 1999)
Miranda, A. C. de O. - "Integra��o de Algoritmos de Gera��o de Malhas de
Elementos Finitos", Disserta��o de Mestrado, Departamento de Engenharia
Civil, PUC/Rio, 1999.

(Owen, 1999)
Owen, S.J., Staten, M.L., Canann, S. A Saigal, 1999. Q-Morph: Na indirect
approach to advancing front quad meshing. International Journal fo Numerical
Methods in Engineering, v. 44, pp. 1317-1340.

(Vianna, 1992)
Vianna, A. C. - "Modelagem Geom�trica Estentida para Modelos Bidimensionais
de Elementos Finitos", Disserta��o de Mestrado, Departamento de Engenharia
Civil, 1992.

*/

/*@}*/

#ifdef __cplusplus
}
#endif

#endif


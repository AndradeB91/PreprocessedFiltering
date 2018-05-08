/** MshSurf.h - header files
* ------------------------------------------------------------------
*/

#ifndef MSHSURFMAPP_H
#define MSHSURFMAPP_H

#include "deflib.h"

#ifdef __cplusplus
extern "C" {
#endif


/** GERAÇÃO DE MALHA ESTRUTURADA DE SUPERFÍCIE 
   @author Antonio Carlos de O. Miranda
   @version 1.0
   @date - april/2010
*/

/**@name Malhas Estruturadas */
/*@{*/

/** Mapeamento transfinito bilinear para regiões quadrilaterais com qualquer 
    forma para os lados. O número de segmentos em lados opostos devem ser 
    iguais (Campos, 1991).
    @param bry       (in) vetor de coordenadas do contorno
    @param m         (in) número de nós na primeira direção
    @param n         (in) número de nós na outra direção
    @param elem_type (in) tipo de elemento:  T3 (3), T6 (6), Q4 (4) ou  Q8 (8)  
    @param diagtype  (in) opção de triangulação:
                     1 - diagonal orientada para direita,
                     2 - diagonal orientada para esquerda,
                     3 - modo cruzado alternado,
                     4 - melhor diagonal
    @param nno       (out) número de nós gerados
    @param nel       (out) número de elementos gerados
    @param pt        (out) coordenadas dos nós da malha
    @param conn      (out) conectividade dos elementos
    @return          1 - sucesso, 0 - falha
*/
MSHSURF_API int MshSurfBilinear (double *bry, int m, int n, int elem_type, int diagtype, 
                     int *nno, int *nel, double **pt, int **conn);
MSHSURF_API int MshSurfTryBilinear (double *bry, int np, int elem_type, int diagtype, 
                        int *nno, int *nel, double **pt, int **conn);


/** Mapeamento transfinito bilinear para regiões triangulares que conceitualmente
    é formado por uma região quadrilateral com um dos lados colapsado para um 
    ponto.  Os três lados podem ter qualquer forma e os lados adjacentes ao 
    lado colapsado devem ter o mesmo número de segmentos (Campos, 1991).
    @param bry       (in) vetor de coordenadas do contorno
    @param m         (in) número de nós na primeira direção
    @param n         (in) número de nós na outra direção
    @param elem_type (in) tipo de elemento:  T3 (3), T6 (6), Q4 (4) ou  Q8 (8)  
    @param diagtype  (in) opção de triangulação
                     1 - diagonal orientada para direita,
                     2 - diagonal orientada para esquerda,
                     3 - modo cruzado alternado,
                     4 - melhor diagonal
    @param nno       (out) número de nós gerados
    @param nel       (out) número de elementos gerados
    @param pt        (out) coordenadas dos nós da malha
    @param conn      (out) conectividade dos elementos
    @return          1 - sucesso, 0 - falha
*/
MSHSURF_API int MshSurfCollBilinear (double *bry, int m, int n, int elem_type, int diagtype,
                        int *nno, int *nel, double **pt, int **conn);


/** Mapeamento transfinito linear entre dois lados opostos de uma região 
    quadrilateral. Os outros dois lados devem ser linhas retas. O número 
    de segmentos nos lados opostos devem ser iguais (Campos, 1991).
    @param bry       (in) vetor de coordenadas do contorno
    @param m         (in) número de nós na primeira direção
    @param n         (in) número de nós na outra direção
    @param dir       (in) direção do lofting: 
                          0 => primeira direção, 1 => outra direção
    @param weigth    (in) peso aplicavél a direçao do lofting
    @param elem_type (in) tipo de elemento:  T3 (3), T6 (6), Q4 (4) ou  Q8 (8)  
    @param diagtype  (in) opção de triangulação
                     1 - diagonal orientada para direita,
                     2 - diagonal orientada para esquerda,
                     3 - modo cruzado alternado,
                     4 - melhor diagonal
    @param nno       (out) número de nós gerados
    @param nel       (out) número de elementos gerados
    @param pt        (out) coordenadas dos nós da malha
    @param conn      (out) conectividade dos elementos
    @return          1 - sucesso, 0 - falha
*/
MSHSURF_API int MshSurfLoft (double *bry, int m, int n, int dir, double weight,
               int elem_type, int diagtype,
               int *nno, int *nel, double **pt, int **conn);


/** Mapeamento transfinito linear entre um ponto e uma curva, que forma uma 
    região triangular. Os dois lados da região que são adjacentes ao ponto 
    de lofting tem de ser linhas retas e ter o mesmo número de pontos 
    (Campos, 1991).
    @param bry       (in) vetor de coordenadas do contorno
    @param m         (in) número de nós na primeira direção
    @param n         (in) número de nós na outra direção
    @param dir       (in) direção do lofting: 
                          0 => primeira direção, 1 => outra direção
    @param weigth    (in) peso aplicavél a direçao do lofting
    @param elem_type (in) tipo de elemento:  T3 (3), T6 (6), Q4 (4) ou  Q8 (8)  
    @param diagtype  (in) opção de triangulação:
                     1 - diagonal orientada para direita,
                     2 - diagonal orientada para esquerda,
                     3 - modo cruzado alternado,
                     4 - melhor diagonal
    @param nno       (out) número de nós gerados
    @param nel       (out) número de elementos gerados
    @param pt        (out) coordenadas dos nós da malha
    @param conn      (out) conectividade dos elementos
    @return          1 - sucesso, 0 - falha
*/  
MSHSURF_API int MshSurfCollLoft (double *bry, int m, int n, double weight,
                   int elem_type, int diagtype,
                   int *nno, int *nel, double **pt, int **conn);


/** Mapeamento transfinito trilinear para uma região triangular com qualquer 
    forma para os lados.  O número de segmentos em todos os lados devem ser 
    iguais (Campos, 1991).
    @param bry       (in) vetor de coordenadas do contorno
    @param m         (in) número de nós em um lado do contorno
    @param elem_type (in) tipo de elemento:  T3 (3), T6 (6), Q4 (4) ou  Q8 (8)  
    @param diagtype  (in) opção de triangulação:
                     1 - diagonal orientada para direita,
                     2 - diagonal orientada para esquerda,
                     3 - modo cruzado alternado,
                     4 - melhor diagonal
    @param nno       (out) número de nós gerados
    @param nel       (out) número de elementos gerados
    @param pt        (out) coordenadas dos nós da malha
    @param conn      (out) conectividade dos elementos
    @return          1 - sucesso, 0 - falha
*/
MSHSURF_API int MshSurfTrilinear (double *bry, int m, int elem_type,  
                      int *nno, int *nel, double **pt, int **conn);
MSHSURF_API int MshSurfTryTrilinear (double *bry, int np, int elem_type,  
                        int *nno, int *nel, double **pt, int **conn);
/*@}*/


/** Dado um conjunto de pontos em 3D, a funcao tenta reordenar os pontos
    de modo a ser possível obter o numero de pontos em u e v em um
    mapeamento bilinear.
    @param np        (in) número de pontos no contorno
    @param pts       (in) vetor de coordenadas do contorno
    @param m         (out) numero de pontos na direcao u
    @param n         (out) numero de pontos na direcao v
    @return          1 - sucesso, 0 - falha
*/
MSHSURF_API int MshSurfAutomaticBilinearCornes (int np, double *pts, int *m, int *n);


/** Dado um conjunto de pontos em 3D, a funcao tenta reordenar os pontos
    de modo a ser possível obter o numero de pontos n em um
    mapeamento trilinear.
    @param np        (in) número de pontos no contorno
    @param pts       (in) vetor de coordenadas do contorno
    @param n         (out) numero de pontos nas direcoes
    @return          1 - sucesso, 0 - falha
*/
MSHSURF_API int MshSurfAutomaticTrilinearCornes (int np, double *pts, int *n);

#ifdef __cplusplus
}
#endif


/*@}*/

#endif


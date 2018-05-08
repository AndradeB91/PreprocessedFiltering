/** MshSurf.h - header files
* ------------------------------------------------------------------
*/

#ifndef MSHSURFMAPP_H
#define MSHSURFMAPP_H

#include "deflib.h"

#ifdef __cplusplus
extern "C" {
#endif


/** GERA��O DE MALHA ESTRUTURADA DE SUPERF�CIE 
   @author Antonio Carlos de O. Miranda
   @version 1.0
   @date - april/2010
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
MSHSURF_API int MshSurfBilinear (double *bry, int m, int n, int elem_type, int diagtype, 
                     int *nno, int *nel, double **pt, int **conn);
MSHSURF_API int MshSurfTryBilinear (double *bry, int np, int elem_type, int diagtype, 
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
MSHSURF_API int MshSurfCollBilinear (double *bry, int m, int n, int elem_type, int diagtype,
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
MSHSURF_API int MshSurfLoft (double *bry, int m, int n, int dir, double weight,
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
MSHSURF_API int MshSurfCollLoft (double *bry, int m, int n, double weight,
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
MSHSURF_API int MshSurfTrilinear (double *bry, int m, int elem_type,  
                      int *nno, int *nel, double **pt, int **conn);
MSHSURF_API int MshSurfTryTrilinear (double *bry, int np, int elem_type,  
                        int *nno, int *nel, double **pt, int **conn);
/*@}*/


/** Dado um conjunto de pontos em 3D, a funcao tenta reordenar os pontos
    de modo a ser poss�vel obter o numero de pontos em u e v em um
    mapeamento bilinear.
    @param np        (in) n�mero de pontos no contorno
    @param pts       (in) vetor de coordenadas do contorno
    @param m         (out) numero de pontos na direcao u
    @param n         (out) numero de pontos na direcao v
    @return          1 - sucesso, 0 - falha
*/
MSHSURF_API int MshSurfAutomaticBilinearCornes (int np, double *pts, int *m, int *n);


/** Dado um conjunto de pontos em 3D, a funcao tenta reordenar os pontos
    de modo a ser poss�vel obter o numero de pontos n em um
    mapeamento trilinear.
    @param np        (in) n�mero de pontos no contorno
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


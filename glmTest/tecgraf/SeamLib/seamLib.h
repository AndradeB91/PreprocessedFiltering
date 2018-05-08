#ifndef SEAMLIB_H
#define SEAMLIB_H

#ifdef __cplusplus
extern "C" {
#endif


/** MIDDLEWARE PARA A BIBLIOTECA COMPUTACIONAL PARA GERAÇÃO DE MALHAS BIDIMENSIONAIS
    DE ELEMENTOS FINITOS(msh2d.lib)

    Gerar a malha bidimensional utilizando o metodo Msh2DQuadSeam, que parte inicialmente
    de uma triangulacao inicial, quando o Msh2DQuadSeam falha, muda-se o fator de qualidade
    da malha(Msh2DSetRefFactor) e tenta-se outra vez.
*/

/** Q-Morph - método indireto de geração de malhas não estruturadas de 
    quadrilaterais (Owen, 1999)
    @param n_loops    (in) número de circuitos
    @param loop_segs  (in) vetor com o número de segmentos em cada  circuito    
    @param bdry_pts   (in) vetor de coordenadas do contorno
    @param type_mesh  (in) Tipo de elemento: 4 (Q4) ou 8 (Q8)
    @param n_node     (out) número de nós gerados
    @param coords     (out) coordenadas dos nós da malha
    @param n_elem     (out) número de elementos gerados
    @param conn       (out) conectividade dos elementos
*/
int MSH2DQuadSeam (int n_loops, int *loop_segs, double *bdry_pts, int type_mesh,
                   int *n_node, double **coords, int *n_elem, int **conn);

int MSH2DBilinear (double *bry, int m, int n, int elem_type, int diagtype, 
                   int *nno, int *nel, double **pt, int **conn);

#ifdef __cplusplus
}
#endif

#endif

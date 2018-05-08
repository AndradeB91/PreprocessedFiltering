/*
**  msh_smoothing.h
*/

#ifndef _MSHSURF_SMOOTHING_H_
#define _MSHSURF_SMOOTHING_H_

#ifdef __cplusplus
extern "C" {
#endif

void SmoothSetLaplacianFactor (double value);

void SmoothLaplacian (void *surf);

/*
** SmoothBiLaplacian
** Paper: Interactive Multiresolution Modeling on Arbitrary Meshes
** L. Kobbelt, S. Campagna, J. Vorsatz, and H-P. Seidel.
** SIGGRAPH 98 Proceedings, 105-114, 1998.
*/
void SmoothBiLaplacian (void *surf);

/*
** SmoothDenoisePoint
** Paper: Bilateral Mesh Denoising - S. Fleishman, I. Droni, and D. Cohen-Or
*/
void SmoothDenoisePoint (void *surf);

/*
** SmoothHC 
** Paper: Improved Laplacian Smoothing of Noisy Surface Meshes,
** J. Vollmer, R. Mencl, and H. Muller, Eurographics'99, 18, 3, 1999
*/
void SmoothHC (void *surf, int num_interaction);

/*
** SmoothMeanCurvFlow
** Paper: Implicit Fairing of Irregular Meshes Using Diffusion and
** Curvature Flow. Computer Graphics (SIGGRAPH 99 Proc), 317-324, 1999.
*/
void SmoothMeanCurvFlow (void *surf);

/*
** SmoothParalelogram
** Use this algorithm only to quadrilateral meshes
** Paper: A New Smoothing Algorithm for Quadrilateral and Hexahedral Meshes.
** Lecture Notes in Computer Science, Volume 3992/2006, 239-246.
*/
void SmoothParallelogram (void *surf, double f1, double f2);


/*
** Smooth2DOptim 
** Smoothing based on optimization of element shape.
** Under construction and tests. 
*/
void Smooth2D_Optim   (void *surf);


void SmoothLocalBezierSurf (void *surf);

void SmoothLocalSplineSurf (void *surf);

#ifdef __cplusplus
}
#endif

#endif


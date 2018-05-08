/* surf3d_octree.h: This function creates an internal octree structure
//                  based on points and sizes inserted in the structure.
////////////////////////////////////////////////////////////////////*/

#ifndef SURF3D_OCTREE_INC
#define SURF3D_OCTREE_INC

#ifdef __cplusplus
extern "C" {
#endif


/* input functions:
** Surf3DOctreeInit          - init the octree with boundbox parameters
** Surf3DOctreeAddVertexSize - refine the octree using a size in a specific point
** Surf3DOctreeEnd           - refine the whole octree to avoid disparity of cell size
*/
void*  Surf3DOctreeInit          (double xmax, double ymax, double zmax,
                                  double xmin, double ymin, double zmin);
void   Surf3DOctreeBound         (double xmax, double ymax, double zmax,
                                  double xmin, double ymin, double zmin);
int  Surf3DOctreeAddPointSize    (void *octree, double x, double y, double z, double size);
void Surf3DOctreeEnd             (void *octree, double max_size, int adj_ref);
void Surf3DOctreeLevelRefine     (void **tree, void *parent, int level);

/* output function:
** Surf3DOctreeSize - obtain the size of octree in a point
** Surf3DOctreeVisitLevels - visite all octree cells
*/
double Surf3DOctreeSize (void *octree, double x, double y, double z);
void   Surf3DOctreeVisitLevels (void *octree, void ((*func)(void *tree, void *parent,
                                int leaf, double min[3], double max[3])));

/* Release memory of internal octree */
void Surf3DOctreeRelease (void *octree);

void Surf3DOctreeDraw (void *octree);

#ifdef __cplusplus
}
#endif


#endif


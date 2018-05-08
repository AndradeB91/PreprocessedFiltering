//
// ArbCrackRegion2D implementation file
//
// Description -
//   This is the implementation file for the
//   ArbMshCrackRegion2D objects.
//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 2000, 2001
//   All rights reserved
//
// Author -
//   Wash Wawrzynek
//
// Revision -
//   $Revision: 1.22 $  $Date: 2002/08/22 19:26:43 $  $Author: wash $
//

#include "ArbMshCrackRegion2D.hpp"
#include "ArbMshRegion2D.hpp"
#include "ArbArray.hpp"
#include "ArbBSpline2D.hpp"
#include "ArbRmshRegion2D.hpp"
#include "ArbSet.cpp"
#include <stdio.h>
#include <math.h>
#include <assert.h>

#ifdef MEMDEBUG
#include "MemDbg.hpp"
#define new new(__FILE__,__LINE__)
#endif

#define PI      3.14159265359
#define HALF_PI 1.57079632680
#define TWO_PI  6.28318530718

#define NO_NODE 100000000

#define DEFAULT_BOUNDARY_TOLERANCE 0.25

// This factor is used when adding crack templates.  If the
// ratio of the template radius divided by the distance to the
// nearest boundary is greater than this ratio, then the template
// size is scaled down so that the template will not overlap or
// be too near a boundary

#define DEFAULT_TEMPLATE_BOUNDARY_FACTOR 0.5
//#define DEFAULT_TEMPLATE_BOUNDARY_FACTOR 0.75

// This factor is used to when setting characteristic element sizes
// for boundary segment ends.  The "insitu" characteristic sizes
// will be multiplied by these values

#define DEFAULT_SEG_END_CHAR_SIZE_FACTOR 1.0

// This value is used to specify the strategy for setting
// characteristic element sizes for boundary segment ends.

//#define DEFAULT_SEG_END_CHAR_SIZE_TYPE S_ANISOTROPIC
#define DEFAULT_SEG_END_CHAR_SIZE_TYPE S_MINIMUM
//#define DEFAULT_SEG_END_CHAR_SIZE_TYPE S_AVERAGE
//#define DEFAULT_SEG_END_CHAR_SIZE_TYPE S_VARIABLE

// This factor is used boundary segments near crack tip are
// updated because a crack is coming too close.  The ratio of
// the node spacing along the boundary divided by the distance
// from the edge of the crack template will be set to be less
// than this value.

#define DEFAULT_BOUNDARY_SPACING_FACTOR 1.1

// The quad angle bounds set the max and min acceptable angles
// in quadrialaterals.  The default behavior is to accept any
// miminum angle but for elements with large angles split the
// quad into two triangles. (angles in radians)

#define DEFAULT_MIN_QUAD_ANGLE 0.0
#define DEFAULT_MAX_QUAD_ANGLE 3.0


FILE *dfd ;
static bool DebugDumpOnly ;

// ------------------------------------------------------------

static int CompareInt(const int &e1,
                           const int &e2)
{
    if (e1 > e2) {
        return(1) ;
    } else if (e1 < e2) {
        return(-1) ;
    }
    return(0) ;
}




// %(CArbMshCrackRegion2D::CArbMshCrackRegion2D-constructor-|)
/* ++ ----------------------------------------------------------
**
**    CArbMshCrackRegion2D - constructor
**
**      CArbMshCrackRegion2D()
**
**      Description: This is the constructor for a ArbMshCrackRegion2D
**          object.
**
**
** -- */

CArbMshCrackRegion2D::CArbMshCrackRegion2D()
{
    NodeTable = new CArbHashTable<int,CArbCoord2D>() ;
    ElemTable = new CArbHashTable<int,ArbMshElement2D>() ;
    CrackTable = new CArbHashTable<int,int>() ;
//    CornerNodes = new CArbSet<int>(ArbCmpUnsigned) ;
    NoRmshMats = new CArbSet<int>(CompareInt) ;
    TemplateElems = new CArbHashTable<int,int>() ;
    RetainedNodes = new CArbHashTable<int,int>() ;

    MshTopo = new CArbMshTopo2D() ;

    ElemShapeType     = S_TRIANGLE ;
    MaxNodeId         = 0 ;
    MaxElemId         = 0 ;
    FirstElem         = true ;
    Order             = LINEAR ;
    CanResizeTemplate = true ;

    DoQuadAngleChecks = true ;
    MinQuadAngle = DEFAULT_MIN_QUAD_ANGLE ;
    MaxQuadAngle = DEFAULT_MAX_QUAD_ANGLE ;
    QuadAngleLocation = AT_NODES ;

    DefaultTipData.elem_shape_set    = true ;
    DefaultTipData.tip_elem_shape    = S_TRIANGLE ;
    DefaultTipData.elem_type_set     = true ;
    DefaultTipData.tip_elem_type     = T_TRIANGLE ;
    DefaultTipData.qrtr_pt_set       = true ;
    DefaultTipData.qrtr_pt_flag      = true ;
    DefaultTipData.const_tip_set     = true ;
    DefaultTipData.const_tip_flag    = true ;
    DefaultTipData.num_tip_set       = true ;
    DefaultTipData.num_tip_elems     = 8 ;
    DefaultTipData.default_num_tip   = true ;
    DefaultTipData.temp_radius_set   = true ;
    DefaultTipData.template_radius   = 0.0 ;
    DefaultTipData.prog_ratio_set    = true ;
    DefaultTipData.progression_ratio = 1.0 ;
    DefaultTipData.num_rings_set     = true ;
    DefaultTipData.number_of_rings   = 2 ;

    TipData  = 0 ;
    DeletedElems = 0 ;
    CreatedElems = 0 ;
    OldBdryTable = 0 ;
    NewBdryTable = 0 ;

    TriDelFactor = 0.0 ;
    QuadDelFactor = 0.0 ;

    OldBdryTable = new CArbHashTable<ArbEdgeKey,BoundaryData> ;
    NewBdryTable = new CArbHashTable<ArbEdgeKey,BoundaryData> ;

    CrackMouthNodes[0] = -1 ;
    CrackMouthNodes[1] = -1 ;

    AnisotropicRatio = 1.0 ;

    BoundaryTolerance = DEFAULT_BOUNDARY_TOLERANCE ;

    TemplateBoundaryFactor = DEFAULT_TEMPLATE_BOUNDARY_FACTOR ;
    SegEndCharSizeFactor   = DEFAULT_SEG_END_CHAR_SIZE_FACTOR ;
    SegEndCharSizeType     = DEFAULT_SEG_END_CHAR_SIZE_TYPE ;
    BoundarySpacingFactor  = DEFAULT_BOUNDARY_SPACING_FACTOR ;

    DebugDisplayFlags = 0 ;

    DebugDumpFlag = false ;

#ifdef DEBUG_DUMP
    DebugDumpFlag = true ;
    dfd = fopen("debug.log","w") ;
    DebugDumpOnly = false ;
#endif
}




// %(CArbMshCrackRegion2D::CArbMshCrackRegion2D-destructor-|~)
/* ++ ----------------------------------------------------------
**
**    CArbMshCrackRegion2D - destructor
**
**      ~CArbMshCrackRegion2D()
**
**      Description: This is the destructor for a ArbMshCrackRegion2D
**          object.
**
**
** -- */

CArbMshCrackRegion2D::~CArbMshCrackRegion2D()
{
    delete NodeTable ;
    delete ElemTable ;
    delete CrackTable ;
    delete TemplateElems ;
    delete RetainedNodes ;
    delete MshTopo ;
    if (NoRmshMats != 0)   delete NoRmshMats ;
    if (TipData != 0)      delete TipData ;
    if (DeletedElems != 0) delete DeletedElems ;
    if (CreatedElems != 0) delete CreatedElems ;
    if (OldBdryTable != 0) delete OldBdryTable ;
    if (NewBdryTable != 0) delete NewBdryTable ;
}




// %(CArbMshCrackRegion2D::AddNode-void-|-ArbMshNode-const|&)
/* ++ ----------------------------------------------------------
**
**    AddNode - add a node to the mesh
**
**      void AddNode(const ArbMshNode &inode)
**
**        inode - (in)  description of the node
**
**      Description: This function adds a new node to the mesh.
**
**      Exceptions:
**          CArbMshCDuplNode - duplicate node id's
**
** -- */

int CArbMshCrackRegion2D::AddNode(const ArbMshNode &inode)
{
    if (NodeTable->Store(inode.id,CArbCoord2D(inode.coord[0],
                                              inode.coord[1]))
                        == false) {
        return(ARB_DUPLICATE_NODE_ID) ;
    }
    if (inode.retained_flag) RetainedNodes->Store(inode.id,1) ;

    if (inode.id > MaxNodeId) MaxNodeId = inode.id ;

    if (DebugDumpFlag)
        fprintf(dfd,"node: %d %18.10g %18.10g\n",
                inode.id,inode.coord[0],inode.coord[1]) ;

    return(ARB_NORMAL_STATUS) ;
}




// %(CArbMshCrackRegion2D::AddNodeList-void-|-int-const|-ArbMshNode-const|*const)
/* ++ ----------------------------------------------------------
**
**    AddNodeList - add an array of nodes to the mesh
**
**      void AddNodeList(
**              const int        num_nodes,
**              const ArbMshNode *constnode_list)
**
**        num_nodes - (in)  number of nodes
**        node_list - (in)  list of nodes to add
**
**      Description: This function adds an array of nodes to the mesh.
**
**      Exceptions:
**          CArbMshCDuplNode - duplicate node id's
**
** -- */

int CArbMshCrackRegion2D::AddNodeList(const int num_nodes,
                         const ArbMshNode *const node_list)
{
    for (int i=0 ; i<num_nodes ; ++i) {
        int status = AddNode(node_list[i]) ;
        if (status != ARB_NORMAL_STATUS) return(status) ;
    }
    return(ARB_NORMAL_STATUS) ;
}




// %(CArbMshCrackRegion2D::AddElem-void-|-ArbMshElement2D-const|&)
/* ++ ----------------------------------------------------------
**
**    AddElem - add an element to the mesh
**
**      void AddElem(const ArbMshElement2D &elem)
**
**        elem - (in)  element to add
**
**      Description: This function adds an element to the mesh.
**
**      Exceptions:
**          CArbMshCDuplElem - duplicate element id's
**          CArbMshCInvalidElem - invalid element
**
** -- */

int CArbMshCrackRegion2D::AddElem(const ArbMshElement2D &ielem)
{
    if (FirstElem) {
        if ((ielem.num_nodes == 3) || (ielem.num_nodes == 4)) {
            Order = LINEAR ;
        } else {
            Order = QUADRATIC ;
        }
        FirstElem = false ;
    } else {
        if ((ielem.num_nodes == 3) || (ielem.num_nodes == 4)) {
            if (Order != LINEAR) return(ARB_INVALID_ELEM) ;
        } else if ((ielem.num_nodes == 6) || (ielem.num_nodes == 8)) {
            if (Order != QUADRATIC) return(ARB_INVALID_ELEM) ;
        } else {
            return(ARB_INVALID_ELEM) ;
        }
    }

    // try to insert this element into the hash table

    if (ElemTable->Store(ielem.elem_id,ielem) == false) {
        return(ARB_DUPLICATE_ELEM_ID) ;
    }

    // insert the element information into the mesh topology

//    int numn = (Order == LINEAR) ?
//                        ielem.num_nodes : ielem.num_nodes/2 ;
//    MshTopo->InsertCollapsedElement(ielem.elem_id,numn,ielem.nodes) ;

    if (Order == LINEAR) {
        MshTopo->InsertCollapsedElement(ielem.elem_id,ielem.num_nodes,
                                        ielem.nodes) ;
    } else {
        int i, nodes[8] ;
        int mid = ielem.num_nodes/2 ;
        for (i=0 ; i<mid ; ++i) nodes[i*2] = ielem.nodes[i] ;
        for (i=0 ; i<mid ; ++i) nodes[i*2+1] = ielem.nodes[mid+i] ;
        MshTopo->InsertCollapsedElement(ielem.elem_id,ielem.num_nodes,nodes) ;
    }

    if (ielem.elem_id > MaxElemId) MaxElemId = ielem.elem_id ;

    if (DebugDumpFlag) {
        fprintf(dfd,"elem: %d %d %d ",ielem.elem_id,ielem.mat_id,
                ielem.num_nodes) ;
        for (int di=0 ; di<ielem.num_nodes ; ++di)
            fprintf(dfd," %d",ielem.nodes[di]) ;
        fprintf(dfd,"\n") ;
    }
    return(ARB_NORMAL_STATUS) ;
}




// %(CArbMshCrackRegion2D::AddElemList-void-|-int-const|-ArbMshElement2D-const|*const)
/* ++ ----------------------------------------------------------
**
**    AddElemList - add an array of elements to the mesh
**
**      void AddElemList(
**              const int             num_elems,
**              const ArbMshElement2D *constelem_list)
**
**        num_elems - (in)  number of elements
**        elem_list - (in)  elements
**
**      Description: This function adds an array of elements to the
**          mesh.
**
**      Exceptions:
**          CArbMshCDuplElem - duplicate element id's
**          CArbMshCInvalidElem - invalid element
**
** -- */

int CArbMshCrackRegion2D::AddElemList(
               const int num_elems,
               const ArbMshElement2D *const elem_list)
{
    for (int i=0 ; i<num_elems ; ++i) {
        int status = AddElem(elem_list[i]) ;
        if (status != ARB_NORMAL_STATUS) return(status) ;
    }
    return(ARB_NORMAL_STATUS) ;
}




// %(CArbMshCrackRegion2D::DefineCrackTip-void-|-int-const|-int-const|)
/* ++ ----------------------------------------------------------
**
**    DefineCrackTip - define a node as a crack tip node
**
**      void DefineCrackTip(
**              const int crack_tip_id,
**              const int node_id)
**
**        crack_tip_id - (in)  crack tip id
**        node_id      - (in)  crack-tip node id
**
**      Description: This function defines a node as a crack tip and
**          associates a crack-tip id with the node.
**
**      Exceptions:
**          CArbMshCDuplCrack - duplicate crack id's
**          CArbMshCInvalidNode - invalid node id
**
** -- */

int CArbMshCrackRegion2D::DefineCrackTip(
                    const int crack_tip_id,
                    const int node_id)
{
    if (NodeTable->Fetch(node_id) == 0) return(ARB_INVALID_NODE) ;

    if (CrackTable->Fetch(crack_tip_id) != 0)
        return(ARB_DUPLICATE_CRACK_ID) ;

    CrackTable->Store(crack_tip_id,node_id) ;

    if (DebugDumpFlag)
        fprintf(dfd,"def_ct: %d %d\n",crack_tip_id,node_id) ;

    return(ARB_NORMAL_STATUS) ;
}




// %(CArbMshCrackRegion2D::SetTipTemplate-void-|-ElemShape-const|-int-const|)
/* ++ ----------------------------------------------------------
**
**    SetTipTemplate - specify crack-tip template element types
**
**      void SetTipTemplate(
**              const ElemShape ishape,
**              const int       crack_tip_id = DEFAULT_ID)
**
**        ishape       - (in)  either S_TRIANGLE or S_QUADRILATERAL
**        crack_tip_id - (in)  id of the crack-tip, or DEFAULT_ID to
**                             indicate that this data should be used
**                             for all tips not specifically identified
**                             in other calls to this routine
**
**      Description: This function allow the type of elements to be
**          used at the crack-tip to be specified.
**
**      Exceptions:
**          CArbMshCInvalidCrack - invalid crack tip id
**
** -- */

int CArbMshCrackRegion2D::SetTipTemplate(
                     const ElemShape ishape,
                     const int crack_tip_id)
{
    if (crack_tip_id == DEFAULT_ID) {
        DefaultTipData.tip_elem_shape = ishape ;
        if (DefaultTipData.default_num_tip)
            DefaultTipData.num_tip_elems =
                (ishape == S_TRIANGLE) ? 8 : 4 ;
    } else {
        CrackTipData *data = GetTipData(crack_tip_id) ;
        if (data == 0) return(ARB_INVALID_CRACK) ;
        data->elem_shape_set = true ;
        data->tip_elem_shape = ishape ;
        if (data->default_num_tip)
            data->num_tip_elems = (ishape == S_TRIANGLE) ? 8 : 4 ;
    }

    if (DebugDumpFlag)
        fprintf(dfd,"tip_temp: %s %d\n",
            (ishape == S_TRIANGLE ? "S_TRIANGLE" : "S_QUADRILATERAL"),
            crack_tip_id) ;

    return(ARB_NORMAL_STATUS) ;
}




// %(CArbMshCrackRegion2D::SetTipElemType-void-|-ElemType-const|-int-const|)
/* ++ ----------------------------------------------------------
**
**    SetTipElemType - specify crack-tip element types
**
**      void SetTipElemType(
**              const ElemType itype,
**              const int      crack_tip_id = DEFAULT_ID)
**
**        itype        - (in)  either T_TRIANGLE or T_COLLAPSED_QUAD
**        crack_tip_id - (in)  id of the crack-tip, or DEFAULT_ID to
**                             indicate that this data should be used
**                             for all tips not specifically identified
**                             in other calls to this routine.
**
**      Description: This function allow the type of elements to be
**          used at the crack-tip to be specified.
**
**      Exceptions:
**          CArbMshCInvalidCrack - invalid crack tip id
**
** -- */

int CArbMshCrackRegion2D::SetTipElemType(
                     const ElemType itype,
                     const int crack_tip_id)
{
    if (crack_tip_id == DEFAULT_ID) {
        DefaultTipData.tip_elem_type = itype ;
    } else {
        CrackTipData *data = GetTipData(crack_tip_id) ;
        if (data == 0) return(ARB_INVALID_CRACK) ;
        data->elem_type_set = true ;
        data->tip_elem_type = itype ;
    }

    if (DebugDumpFlag)
        fprintf(dfd,"tip_type: %s %d\n",
            (itype == T_TRIANGLE ? "T_TRIANGLE" : "T_COLLAPSED_QUAD"),
            crack_tip_id) ;

    return(ARB_NORMAL_STATUS) ;
}




// %(CArbMshCrackRegion2D::SetQrtrPointTip-void-|-bool-const|-int-const|)
/* ++ ----------------------------------------------------------
**
**    SetQrtrPointTip - specify quarter-point crack tip elements
**
**      void SetQrtrPointTip(
**              const bool flag = true,
**              const int  crack_tip_id = DEFAULT_ID)
**
**        flag         - (in)  true means use quarter-point elements
**        crack_tip_id - (in)  id of the crack-tip, or DEFAULT_ID to
**                             indicate that this data should be used
**                             for all tips not specifically identified
**                             in other calls to this routine.
**
**      Description: This function sets a flag that controls if
**          crack-tip elements will be quarter-point elements in
**          quadratic order meshes. The default is to use quarter-point
**          elements.
**
**      Exceptions:
**          CArbMshCInvalidCrack - invalid crack id
**
** -- */

int CArbMshCrackRegion2D::SetQrtrPointTip(
                       const bool flag,
                       const int crack_tip_id)
{
    if (crack_tip_id == DEFAULT_ID) {
        DefaultTipData.qrtr_pt_flag = flag ;
    } else {
        CrackTipData *data = GetTipData(crack_tip_id) ;
        if (data == 0) return(ARB_INVALID_CRACK) ;
        data->qrtr_pt_set = true ;
        data->qrtr_pt_flag = flag ;
    }

    if (DebugDumpFlag) fprintf(dfd,"qpt_tip: %s %d\n",
            (flag ? "TRUE" : "FALSE"),crack_tip_id) ;

    return(ARB_NORMAL_STATUS) ;
}




// %(CArbMshCrackRegion2D::SetConstrainedTip-void-|-bool-const|-int-const|)
/* ++ ----------------------------------------------------------
**
**    SetConstrainedTip - specify constrained collapsed crack-tip
**                        elements
**
**      void SetConstrainedTip(
**              const bool flag = true,
**              const int  crack_tip_id = DEFAULT_ID)
**
**        flag         - (in)  true means constrain the crack-tip
**                             elements
**        crack_tip_id - (in)  id of the crack-tip, or DEFAULT_ID to
**                             indicate that this data should be used
**                             for all tips not specifically identified
**                             in other calls to this routine.
**
**      Description: This function sets a flag that controls if
**          triangular crack-tip elements, which are modeled as
**          collapsed quadrilateral elements, have the nodes at the
**          crack-tip to be constrained to move together. The default
**          is to constrain the nodes.
**
**      Exceptions:
**          CArbMshCInvalidCrack - invalid crack id
**
** -- */

int CArbMshCrackRegion2D::SetConstrainedTip(
                       const bool flag,
                       const int crack_tip_id)
{
    if (crack_tip_id == DEFAULT_ID) {
        DefaultTipData.const_tip_flag = flag ;
    } else {
        CrackTipData *data = GetTipData(crack_tip_id) ;
        if (data == 0) return(ARB_INVALID_CRACK) ;
        data->const_tip_set = true ;
        data->const_tip_flag = flag ;
    }

    if (DebugDumpFlag) fprintf(dfd,"cons_tip: %s %d\n",
            (flag ? "TRUE" : "FALSE"),crack_tip_id) ;

    return(ARB_NORMAL_STATUS) ;
}




// %(CArbMshCrackRegion2D::SetNumTipElems-void-|-int-const|-int-const|)
/* ++ ----------------------------------------------------------
**
**    SetNumTipElems - set the number of crack-tip elements
**
**      void SetNumTipElems(
**              const int num_elem,
**              const int crack_tip_id = DEFAULT_ID)
**
**        num_elem     - (in)  number of crack tip elements
**        crack_tip_id - (in)  id of the crack-tip, or DEFAULT_ID to
**                             indicate that this data should be used
**                             for all tips not specifically identified
**                             in other calls to this routine.
**
**      Description: This function sets the number of elements to
**          insert at the crack-tip. The default is 8 for triangular
**          elements and 4 for quadrilateral elements.
**
**      Exceptions:
**          CArbMshCInvalidCrack - invalid crack id
**
** -- */

int CArbMshCrackRegion2D::SetNumTipElems(
                       const int num_elem,
                       const int crack_tip_id)
{
    if (crack_tip_id == DEFAULT_ID) {
        DefaultTipData.num_tip_elems = num_elem ;
        DefaultTipData.default_num_tip = false ;
    } else {
        CrackTipData *data = GetTipData(crack_tip_id) ;
        if (data == 0) return(ARB_INVALID_CRACK) ;
        data->num_tip_set = true ;
        data->num_tip_elems = num_elem ;
        data->default_num_tip = false ;
    }

    if (DebugDumpFlag)
        fprintf(dfd,"num_tip_elem: %d %d\n",num_elem,crack_tip_id) ;

    return(ARB_NORMAL_STATUS) ;
}




// %(CArbMshCrackRegion2D::SetTemplateRadius-void-|-double-const|-int-const|)
/* ++ ----------------------------------------------------------
**
**    SetTemplateRadius - set the crack-tip template radius
**
**      void SetTemplateRadius(
**              const double radius,
**              const int    crack_tip_id = DEFAULT_ID)
**
**        radius       - (in)  template radius
**        crack_tip_id - (in)  id of the crack-tip, or DEFAULT_ID to
**                             indicate that this data should be used
**                             for all tips not specifically identified
**                             in other calls to this routine.
**
**      Description: This function sets the radius of the crack-tip
**          element template.
**
**      Exceptions:
**          CArbMshCInvalidCrack - invalid crack id
**
** -- */

int CArbMshCrackRegion2D::SetTemplateRadius(
                    double radius,
                    const int crack_tip_id)
{
    if (crack_tip_id == DEFAULT_ID) {
        DefaultTipData.template_radius = radius ;
    } else {
        CrackTipData *data = GetTipData(crack_tip_id) ;
        if (data == 0) return(ARB_INVALID_CRACK) ;
        data->temp_radius_set = true ;
        data->template_radius = radius ;
    }

    if (DebugDumpFlag)
        fprintf(dfd,"temp_rad: %18.10f %d\n",radius,crack_tip_id) ;

    return(ARB_NORMAL_STATUS) ;
}




// %(CArbMshCrackRegion2D::SetTemplateRatio-void-|-double-const|)
/* ++ ----------------------------------------------------------
**
**    SetTemplateRatio - sets the crack-tip template progression ratio
**
**      void SetTemplateRatio(const double ratio = DEFAULT_ID)
**
**        ratio - (in)  id of the crack-tip, or DEFAULT_ID to indicate
**                      that this data should be used for all tips not
**                      specifically identified in other calls to this
**                      routine.
**
**      Description: This function sets the progression ratio for the
**          crack-tip template elements.
**
**      Exceptions:
**          CArbMshCInvalidCrack - invalid crack id
**
** -- */

int CArbMshCrackRegion2D::SetTemplateRatio(
                    double ratio,
                    const int crack_tip_id)
{
    if (crack_tip_id == DEFAULT_ID) {
        DefaultTipData.progression_ratio = ratio ;
    } else {
        CrackTipData *data = GetTipData(crack_tip_id) ;
        if (data == 0) return(ARB_INVALID_CRACK) ;
        data->prog_ratio_set = true ;
        data->progression_ratio = ratio ;
    }

    if (DebugDumpFlag)
        fprintf(dfd,"temp_ratio: %18.10f %d\n",ratio,crack_tip_id) ;

    return(ARB_NORMAL_STATUS) ;
}




// %(CArbMshCrackRegion2D::SetTemplateNumRings-void-|-int-const|-int-const|)
/* ++ ----------------------------------------------------------
**
**    SetTemplateNumRings - set the number of rings in the crack-tip
**                          template
**
**      void SetTemplateNumRings(
**              const int num_rings,
**              const int crack_tip_id = DEFAULT_ID)
**
**        num_rings    - (in)  number of element rings
**        crack_tip_id - (in)  id of the crack-tip, or DEFAULT_ID to
**                             indicate that this data should be used
**                             for all tips not specifically identified
**                             in other calls to this routine.
**
**      Description: This function sets the number of element rings
**          placed in a crack-tip template.
**
**      Exceptions:
**          CArbMshCInvalidCrack - invalid crack id
**
** -- */

int CArbMshCrackRegion2D::SetTemplateNumRings(
                    const int num_rings,
                    const int crack_tip_id)
{
    if (crack_tip_id == DEFAULT_ID) {
        DefaultTipData.number_of_rings = num_rings ;
    } else {
        CrackTipData *data = GetTipData(crack_tip_id) ;
        if (data == 0) return(ARB_INVALID_CRACK) ;
        data->num_rings_set = true ;
        data->number_of_rings = num_rings ;
    }

    if (DebugDumpFlag)
        fprintf(dfd,"num_rings: %d %d\n",num_rings,crack_tip_id) ;

    return(ARB_NORMAL_STATUS) ;
}




// %(CArbMshCrackRegion2D::SetRemeshElemShape-void-|-ElemShape-const|)
/* ++ ----------------------------------------------------------
**
**    SetRemeshElemShape - specify the type of elements used for
**                         remeshing
**
**      void SetRemeshElemShape(const ElemShape ishape)
**
**        ishape - (in)  either S_TRIANGLE or S_QUADRILATERAL
**
**      Description: This function specifies the type (triangle or
**          quadrilateral) of elements used during remeshing.
**
**
** -- */

void CArbMshCrackRegion2D::SetRemeshElemShape(const ElemShape shape)
{
    ElemShapeType = shape ;

    if (DebugDumpFlag)
        fprintf(dfd,"rmsh_elem_shape: %s\n",
            (shape == S_TRIANGLE ? "S_TRIANGLE" : "S_QUADRILATERAL")) ;
}




// %(CArbMshCrackRegion2D::SetNoRemeshMat-void-|-int-const|)
/* ++ ----------------------------------------------------------
**
**    SetNoRemeshMat - specify a "no remesh" material
**
**      void SetNoRemeshMat(const int mat_id)
**
**        mat_id - (in)  id of a material that should not be modified
**
**      Description: This function sets a flag that indicates that no
**          elements with the specified material id will be deleted or
**          modified during crack insertion or growth.
**
**
** -- */

void CArbMshCrackRegion2D::SetNoRemeshMat(const int mat_id)
{
  int mat_tmp = mat_id;
    NoRmshMats->Insert(mat_tmp) ;

    if (DebugDumpFlag)
        fprintf(dfd,"no_rmsh_mat: %d\n",mat_id) ;
}




// %(CArbMshCrackRegion2D::SetTemplateElem-void-|-int-const|-int-const|)
/* ++ ----------------------------------------------------------
**
**    SetTemplateElem - flag an element as a template element
**
**      void SetTemplateElem(
**              const int elem_id,
**              const int tip_id)
**
**        elem_id - (in)  id of an element
**        tip_id  - (in)  id of a crack tip
**
**      Description: This function flags an element as a template
**          element associated with the specified crack tip. Template
**          elements are not deleted during remeshing of other crack
**          tips.
**
**
** -- */

void CArbMshCrackRegion2D::SetTemplateElem(const int elem_id,
                                           const int tip_id)
{
    TemplateElems->Store(elem_id,tip_id) ;

    if (DebugDumpFlag)
        fprintf(dfd,"set_temp_elem: %d %d\n",elem_id,tip_id) ;
}


/* ++ ----------------------------------------------------------
**
**    SetQuadAngleChecks - turns on checks on max and min quad angles
**
**      void SetQuadElemAngleChecks(double min_angle,
**                                  double max_angle,
**                                  AngleCheckLocation location)
**
**        min_angle - (in) minimum acceptable angle (radians)
**        max_angle - (in) maximum acceptable angle (radians)
**        location  - (in) specify if angles should be checked
**                         at the nodes or at the integration points.
**
**      Description: This turns on quadrilateral angle checking and
**          sets the check parameters.  If a quadrilateral element
**          is found to have angles that are out of bounds it is
**          replaced with two triangular elements.
**
** -- */

void CArbMshCrackRegion2D::SetQuadElemAngleChecks(
                                double min_angle,
                                double max_angle,
                                AngleCheckLocation location)
{
    DoQuadAngleChecks = true ;
    MinQuadAngle = min_angle ;
    MaxQuadAngle = max_angle ;
    QuadAngleLocation = location ;
}


// %(CArbMshCrackRegion2D::GetTipTemplate-ElemShape-|-int-const|)
/* ++ ----------------------------------------------------------
**
**    GetTipTemplate - return the template element shapes
**
**      ElemShape GetTipTemplate(const int crack_tip_id = DEFAULT_ID)
**
**        crack_tip_id - (in)  crack-tip id or DEFAULT_ID
**
**      Description: This function returns the shape (quadrilateral or
**          triangle) of elements used at the crack-tip.
**
**      Return Value: either S_TRIANGLE or S_QUADRILATERAL
**
**      Exceptions:
**          CArbMshCInvalidCrack - invalid crack-tip id
**
** -- */

CArbMshCrackRegion2D::ElemShape
CArbMshCrackRegion2D::GetTipTemplate(const int crack_tip_id) const
{
    if ((crack_tip_id != DEFAULT_ID) &&
        (TipData != 0)) {
        CrackTipData *data = TipData->Fetch(crack_tip_id) ;
        if ((data != 0) && (data->elem_shape_set))
            return(data->tip_elem_shape) ;
    }
    return(DefaultTipData.tip_elem_shape) ;
}




// %(CArbMshCrackRegion2D::GetTipElemType-ElemType-|-int-const|)
/* ++ ----------------------------------------------------------
**
**    GetTipElemType - return the template element types
**
**      ElemType GetTipElemType(const int crack_tip_id = DEFAULT_ID)
**
**        crack_tip_id - (in)  crack-tip id or DEFAULT_ID
**
**      Description: This function returns the type (triangle or
**          collapsed quadrilateral) of elements used at the crack-tip.
**
**      Return Value: either T_TRIANGLE or T_COLLAPSED_QUAD
**
**      Exceptions:
**          CArbMshCInvalidCrack - invalid crack-tip id
**
** -- */

CArbMshCrackRegion2D::ElemType
CArbMshCrackRegion2D::GetTipElemType(const int crack_tip_id) const
{
    if ((crack_tip_id != DEFAULT_ID) &&
        (TipData != 0)) {
        CrackTipData *data = TipData->Fetch(crack_tip_id) ;
        if ((data != 0) && (data->elem_type_set))
            return(data->tip_elem_type) ;
    }
    return(DefaultTipData.tip_elem_type) ;
}




// %(CArbMshCrackRegion2D::GetQrtrPointTip-bool-|-int-const|)
/* ++ ----------------------------------------------------------
**
**    GetQrtrPointTip - return the quarter-point tip flag
**
**      bool GetQrtrPointTip(const int crack_tip_id = DEFAULT_ID)
**
**        crack_tip_id - (in)  crack-tip id or DEFAULT_ID
**
**      Description: This function returns the flag that indicatese if
**          quarter-point elements will be inserted at a crack tip for
**          quadratic order meshes.
**
**      Return Value: true means quarter-point elements will be used
**
**      Exceptions:
**          CArbMshCInvalidCrack - invalid crack-tip id
**
** -- */

bool CArbMshCrackRegion2D::GetQrtrPointTip(const int crack_tip_id) const
{
    if ((crack_tip_id != DEFAULT_ID) &&
        (TipData != 0)) {
        CrackTipData *data = TipData->Fetch(crack_tip_id) ;
        if ((data != 0) && (data->qrtr_pt_set))
            return(data->qrtr_pt_flag) ;
    }
    return(DefaultTipData.qrtr_pt_flag) ;
}




// %(CArbMshCrackRegion2D::GetConstrainedTip-bool-|-int-const|)
/* ++ ----------------------------------------------------------
**
**    GetConstrainedTip - return the constrained tip flag
**
**      bool GetConstrainedTip(const int crack_tip_id = DEFAULT_ID)
**
**        crack_tip_id - (in)  crack-tip id or DEFAULT_ID
**
**      Description: This function returns the flag that indicatese if
**          crack-tip nodes for collapsed quadrilateral crack-tip
**          elements will be constrained to move together.
**
**      Return Value: true means that the crack-tip elements will be
**          constrained to move together
**
**      Exceptions:
**          CArbMshCInvalidCrack - invalid crack-tip id
**
** -- */

bool CArbMshCrackRegion2D::GetConstrainedTip(const int crack_tip_id) const
{
    if ((crack_tip_id != DEFAULT_ID) &&
        (TipData != 0)) {
        CrackTipData *data = TipData->Fetch(crack_tip_id) ;
        if ((data != 0) && (data->const_tip_set))
            return(data->const_tip_flag) ;
    }
    return(DefaultTipData.const_tip_flag) ;
}




// %(CArbMshCrackRegion2D::GetNumTipElems-int-|-int-const|)
/* ++ ----------------------------------------------------------
**
**    GetNumTipElems - return the number of crack-tip elements
**
**      int GetNumTipElems(const int crack_tip_id = DEFAULT_ID)
**
**        crack_tip_id - (in)  crack-tip id or DEFAULT_ID
**
**      Description: This function returns the number of crack-tip
**          elements that will be inserted at a crack tip.
**
**      Return Value: number of elements
**
**      Exceptions:
**          CArbMshCInvalidCrack - invalid crack-tip id
**
** -- */

int CArbMshCrackRegion2D::GetNumTipElems(
                    const int crack_tip_id) const
{
    if ((crack_tip_id != DEFAULT_ID) &&
        (TipData != 0)) {
        CrackTipData *data = TipData->Fetch(crack_tip_id) ;
        if ((data != 0) && (data->num_tip_set))
            return(data->num_tip_elems) ;
    }
    return(DefaultTipData.num_tip_elems) ;
}




// %(CArbMshCrackRegion2D::GetTemplateRadius-double-|-int-const|)
/* ++ ----------------------------------------------------------
**
**    GetTemplateRadius - return the crack-tip radius
**
**      double GetTemplateRadius(const int crack_tip_id = DEFAULT_ID)
**
**        crack_tip_id - (in)  crack-tip id or DEFAULT_ID
**
**      Description: This function returns the radius of the crack-tip
**          template inserted at a crack tip.
**
**      Return Value: the template radius
**
**      Exceptions:
**          CArbMshCInvalidCrack - invalid crack-tip id
**
** -- */

double CArbMshCrackRegion2D::GetTemplateRadius(
                    const int crack_tip_id) const
{
    if ((crack_tip_id != DEFAULT_ID) &&
        (TipData != 0)) {
        CrackTipData *data = TipData->Fetch(crack_tip_id) ;
        if ((data != 0) && (data->temp_radius_set))
            return(data->template_radius) ;
    }
    return(DefaultTipData.template_radius) ;
}




// %(CArbMshCrackRegion2D::GetTemplateRatio-double-|-int-const|)
/* ++ ----------------------------------------------------------
**
**    GetTemplateRatio - return the progression ratio
**
**      double GetTemplateRatio(const int crack_tip_id = DEFAULT_ID)
**
**        crack_tip_id - (in)  crack-tip id or DEFAULT_ID
**
**      Description: This function returns the element size progression
**          ratio that will be used when inserting elements at a crack
**          tip.
**
**      Return Value: progression ratio
**
**      Exceptions:
**          CArbMshCInvalidCrack - invalid crack-tip id
**
** -- */

double CArbMshCrackRegion2D::GetTemplateRatio(
                    const int crack_tip_id) const
{
    if ((crack_tip_id != DEFAULT_ID) &&
        (TipData != 0)) {
        CrackTipData *data = TipData->Fetch(crack_tip_id) ;
        if ((data != 0) && (data->prog_ratio_set))
            return(data->progression_ratio) ;
    }
    return(DefaultTipData.progression_ratio) ;
}




// %(CArbMshCrackRegion2D::GetTemplateNumRings-int-|-int-const|)
/* ++ ----------------------------------------------------------
**
**    GetTemplateNumRings - return the number of template element rings
**
**      int GetTemplateNumRings(const int crack_tip_id = DEFAULT_ID)
**
**        crack_tip_id - (in)  crack-tip id or DEFAULT_ID
**
**      Description: This function returns the number of element rings
**          that will be used when generating a crack-tip template.
**
**      Return Value: number of element rings
**
**      Exceptions:
**          CArbMshCInvalidCrack - invalid crack-tip id
**
** -- */

int CArbMshCrackRegion2D::GetTemplateNumRings(
                    const int crack_tip_id) const
{
    if ((crack_tip_id != DEFAULT_ID) &&
        (TipData != 0)) {
        CrackTipData *data = TipData->Fetch(crack_tip_id) ;
        if ((data != 0) && (data->num_rings_set))
            return(data->number_of_rings) ;
    }
    return(DefaultTipData.number_of_rings) ;
}




// %(CArbMshCrackRegion2D::GetRemeshElemShape-ElemShape-|^const)
/* ++ ----------------------------------------------------------
**
**    GetRemeshElemShape - returns the element shape used during
**                         remeshing
**
**      ElemShape GetRemeshElemShape() const
**
**      Description: This function returns the shape of the elements
**          (triangle or quadrilateral) that will be used during
**          remeshing.
**
**      Return Value: either S_TRIANGULAR or S_QUADRILATERAL
**
**
** -- */

CArbMshCrackRegion2D::ElemShape
CArbMshCrackRegion2D::GetRemeshElemShape() const
{
    return(ElemShapeType) ;
}




// %(CArbMshCrackRegion2D::GetNoRemeshMats-int-|*^const-int-|*)
/* ++ ----------------------------------------------------------
**
**    GetNoRemeshMats - return a list of the "no remesh" materials
**
**      int *GetNoRemeshMats(int *num) const
**
**        num - (out) number of id's in the list
**
**      Description: This function returns a list of the id's of the
**          materials that will not be modified during remeshing.
**
**      Return Value: An array of element id's. Owner should of this
**          list passes to the client, which should eventially call
**          delete [].
**
**
** -- */

int *CArbMshCrackRegion2D::GetNoRemeshMats(int *num) const
{
    if (NoRmshMats == 0) {
        *num = 0 ;
        return(0) ;
    } else {
        *num = NoRmshMats->NumEntries() ;
        int *list = new int[*num] ;
        int **keys = NoRmshMats->GetKeyList() ;
        for (int i=0 ; i<(*num) ; ++i)
            list[i] = *(keys[i]) ;
        delete [] keys ;
        return(list) ;
    }
}




// %(CArbMshCrackRegion2D::GetElementOrder-ArbMshOrder-|^const)
/* ++ ----------------------------------------------------------
**
**    GetElementOrder - return the polynomial order of the mesh
**                      elements
**
**      ArbMshOrder GetElementOrder() const
**
**      Description: This function returns the polynomial order of the
**          elements used in the mesh.
**
**      Return Value: either LINEAR or QUADRATIC
**
**
** -- */

ArbMshOrder CArbMshCrackRegion2D::GetElementOrder() const
{
    return(Order) ;
}


/* ------------------------------------------------------------ */

void display_topo(CArbMshTopo2D *topo,
             CArbHashTable<int,CArbCoord2D> *NodeTable)
{
    int num_elem = topo->NumElements() ;
    int *elems = topo->GetElemList() ;

    printf("# Start disp topo\n") ;

    for (int i=0 ; i<num_elem ; ++i) {
        CArbTopoVtxOnElemIterator iter(topo,elems[i]) ;
        int first = *iter ;
        int cur = first ;

        while (iter.More()) {
            int next = *iter ;
            CArbCoord2D *nd0 = NodeTable->Fetch(cur) ;
            CArbCoord2D *nd1 = NodeTable->Fetch(next) ;
            printf("# edge %d %d\n",cur,next) ;
            printf("l %g %g %g %g\n",
                   (*nd0)[0],(*nd0)[1],(*nd1)[0],(*nd1)[1]) ;
            printf("t %g %g %d\n",(*nd0)[0],(*nd0)[1],cur) ;
            cur = next ;
            ++iter ;
        }

        CArbCoord2D *nd0 = NodeTable->Fetch(cur) ;
        CArbCoord2D *nd1 = NodeTable->Fetch(first) ;
        printf("# edge %d %d\n",cur,first) ;
        printf("l %g %g %g %g\n",
               (*nd0)[0],(*nd0)[1],(*nd1)[0],(*nd1)[1]) ;
        printf("t %g %g %d\n",(*nd0)[0],(*nd0)[1],cur) ;
    }
    printf("n\n") ;
    fflush(stdout) ;

    delete [] elems ;
}

// ----------------------------------------------------------

int CArbMshCrackRegion2D::NewInternalCrack(
                       const int start_tip_id,
                       const int end_tip_id,
                       const int num_points,
                       const CrackPt *const pts)
{
    if (DebugDumpFlag) {
        fprintf(dfd,"new_internal_crack: %d %d %d",start_tip_id,
                end_tip_id,num_points) ;
        for (int id=0 ; id<num_points ; ++id)
            fprintf(dfd," [ %18.10f %18.10f %18.10f %d ]",
                pts[id].coord[0],pts[id].coord[0],pts[id].char_elem_size,
                (pts[id].has_char_size ? 1 : 0)) ;
        fprintf(dfd,"\n") ;
        fflush(dfd) ;
        if (DebugDumpOnly) return(ARB_NORMAL_STATUS) ;
    }

    FirstNewElem = MaxElemId+1 ;
    CArbRmshRegion2D rmsh(this) ;
    rmsh.SetDebugDisplayFlags(
        CArbRmshRegion2D::DDFlags(DebugDisplayFlags)) ;
    if (TriDelFactor > 0.0)  rmsh.SetTriDeleteFactor(TriDelFactor) ;
    if (QuadDelFactor > 0.0) rmsh.SetQuadDeleteFactor(QuadDelFactor) ;
    if (DoQuadAngleChecks) {
        CArbRmshRegion2D::AngleCheckLocation loc =
            QuadAngleLocation == AT_NODES ?
            CArbRmshRegion2D::AT_NODES :
            CArbRmshRegion2D::AT_INT_PTS ;
        rmsh.SetQuadElemAngleChecks(MinQuadAngle,MaxQuadAngle,loc) ;
    }
    rmsh.SetAnisotropicRatio(AnisotropicRatio) ;

    // take the linear description of the crack and
    // transform that into a (zero area) polygonal description

    int num = num_points * 2 - 2 ;
    CrackPt *ppts = new CrackPt[num] ;

    int cur = 0 ;
    for (int i=0 ; i<num_points ; ++i) {
        ppts[cur] = pts[i] ;  ++cur ;
    }

    for (int j=num_points-2 ; j>0 ; --j) {
        ppts[cur] = pts[j] ;  ++cur ;
    }

    int tip_ids[2] = {start_tip_id,end_tip_id} ;
    int tip_indx[2] = {0,num_points-1} ;
    rmsh.AddFlaw(2,tip_ids,tip_indx,num,ppts) ;
    delete [] ppts ;

    return(ARB_NORMAL_STATUS) ;
}


// ----------------------------------------------------------

int CArbMshCrackRegion2D::NewSurfaceCrack(
                        const int tip_id,
                        const int num_points,
                        const CrackPt *const pts)
{
    if (DebugDumpFlag) {
        fprintf(dfd,"new_surface_crack: %d %d",tip_id,num_points) ;
        for (int id=0 ; id<num_points ; ++id)
            fprintf(dfd," [ %18.10f %18.10f %18.10f %d ]",
                pts[id].coord[0],pts[id].coord[0],pts[id].char_elem_size,
                (pts[id].has_char_size ? 1 : 0)) ;
        fprintf(dfd,"\n") ;
        fflush(dfd) ;
        if (DebugDumpOnly) return(ARB_NORMAL_STATUS) ;
    }

    FirstNewElem = MaxElemId+1 ;
    CArbRmshRegion2D rmsh(this) ;
    rmsh.SetDebugDisplayFlags(
        CArbRmshRegion2D::DDFlags(DebugDisplayFlags)) ;
    if (TriDelFactor > 0.0)  rmsh.SetTriDeleteFactor(TriDelFactor) ;
    if (QuadDelFactor > 0.0) rmsh.SetQuadDeleteFactor(QuadDelFactor) ;
    if (DoQuadAngleChecks) {
        CArbRmshRegion2D::AngleCheckLocation loc =
            QuadAngleLocation == AT_NODES ?
            CArbRmshRegion2D::AT_NODES :
            CArbRmshRegion2D::AT_INT_PTS ;
        rmsh.SetQuadElemAngleChecks(MinQuadAngle,MaxQuadAngle,loc) ;
    }
    rmsh.SetAnisotropicRatio(AnisotropicRatio) ;

    // take the linear description of the crack and
    // transform that into a (zero area) polygonal description

    int num = num_points * 2 - 2 ;
    CrackPt *ppts = new CrackPt[num] ;

    int cur = 0 ;
    for (int i=0 ; i<num_points ; ++i) {
        ppts[cur] = pts[i] ;  ++cur ;
    }

    for (int j=num_points-2 ; j>0 ; --j) {
        ppts[cur] = pts[j] ;  ++cur ;
    }

    int indx = num_points-1 ;
    rmsh.AddFlaw(1,&tip_id,&indx,num,ppts) ;
    delete [] ppts ;

    CArbArray<int> *MouthNodes = rmsh.GetMouthNodes() ;
    if (MouthNodes != 0) {
        if (MouthNodes->NumEntries() > 0)
            CrackMouthNodes[0] = MouthNodes->At(0) ;
        if (MouthNodes->NumEntries() > 1)
            CrackMouthNodes[1] = MouthNodes->At(1) ;
    }

    return(ARB_NORMAL_STATUS) ;
}



// ----------------------------------------------------------


void CArbMshCrackRegion2D::NewGeneralFlaw(
                         const int num_tips,
                         const int *tip_ids,
                         const int *tip_indices,
                         const int num_points,
                         const CrackPt *const pts)
{
    if (DebugDumpFlag) {
        fprintf(dfd,"new_general_flaw: %d %d [",num_tips,num_points) ;
        for (int di0=0 ; di0<num_tips ; ++di0)
            fprintf(dfd," %d",tip_ids[di0]) ;
        fprintf(dfd," ] [") ;
        for (int di1=0 ; di1<num_tips ; ++di1)
            fprintf(dfd," %d",tip_indices[di1]) ;
        fprintf(dfd," ]") ;
        for (int di2=0 ; di2<num_points ; ++di2)
            fprintf(dfd," [ %18.10f %18.10f %18.10f %d ]",
                pts[di2].coord[0],pts[di2].coord[0],pts[di2].char_elem_size,
                (pts[di2].has_char_size ? 1 : 0)) ;
        fprintf(dfd,"\n") ;
        fflush(dfd) ;
        if (DebugDumpOnly) return ;
    }

    FirstNewElem = MaxElemId+1 ;
    CArbRmshRegion2D rmsh(this) ;
    rmsh.SetDebugDisplayFlags(
        CArbRmshRegion2D::DDFlags(DebugDisplayFlags)) ;
    if (TriDelFactor > 0.0)  rmsh.SetTriDeleteFactor(TriDelFactor) ;
    if (QuadDelFactor > 0.0) rmsh.SetQuadDeleteFactor(QuadDelFactor) ;
    if (DoQuadAngleChecks) {
        CArbRmshRegion2D::AngleCheckLocation loc =
            QuadAngleLocation == AT_NODES ?
            CArbRmshRegion2D::AT_NODES :
            CArbRmshRegion2D::AT_INT_PTS ;
        rmsh.SetQuadElemAngleChecks(MinQuadAngle,MaxQuadAngle,loc) ;
    }
    rmsh.SetAnisotropicRatio(AnisotropicRatio) ;
    rmsh.AddFlaw(num_tips,tip_ids,tip_indices,num_points,pts) ;
}


// ----------------------------------------------------------

void CArbMshCrackRegion2D::CrackGrowth(
                       const int tip_id,
                       const int num_points,
                       const CrackPt *const pts)
{
    if (DebugDumpFlag) {
        fprintf(dfd,"crack_growth: %d %d",tip_id,num_points) ;
        for (int id=0 ; id<num_points ; ++id)
            fprintf(dfd," [ %18.10f %18.10f %18.10f %d ]",
                pts[id].coord[0],pts[id].coord[0],pts[id].char_elem_size,
                (pts[id].has_char_size ? 1 : 0)) ;
        fprintf(dfd,"\n") ;
        fflush(dfd) ;
        if (DebugDumpOnly) return ;
    }

    FirstNewElem = MaxElemId+1 ;
    CArbRmshRegion2D rmsh(this) ;
    rmsh.SetDebugDisplayFlags(
        CArbRmshRegion2D::DDFlags(DebugDisplayFlags)) ;
    if (TriDelFactor > 0.0)  rmsh.SetTriDeleteFactor(TriDelFactor) ;
    if (QuadDelFactor > 0.0) rmsh.SetQuadDeleteFactor(QuadDelFactor) ;
    if (DoQuadAngleChecks) {
        CArbRmshRegion2D::AngleCheckLocation loc =
            QuadAngleLocation == AT_NODES ?
            CArbRmshRegion2D::AT_NODES :
            CArbRmshRegion2D::AT_INT_PTS ;
        rmsh.SetQuadElemAngleChecks(MinQuadAngle,MaxQuadAngle,loc) ;
    }
    rmsh.SetAnisotropicRatio(AnisotropicRatio) ;
    CArbCoord2D old_tip(pts[0].coord[0],pts[0].coord[1]) ;
    rmsh.SetOldCrackTip(old_tip) ;

    // take the linear description of the crack and
    // transform that into a (zero area) polygonal description
    // adding the crack-tip coordinate

    int num = num_points * 2 - 2 ;
    CrackPt *ppts = new CrackPt[num] ;

    int cur = 0 ;
    for (int i=0 ; i<num_points ; ++i) {
        ppts[cur] = pts[i] ;  ++cur ;
    }

    for (int j=num_points-2 ; j>0 ; --j) {
        ppts[cur] = pts[j] ;  ++cur ;
    }

    int indx = num_points-1 ;
    rmsh.AddFlaw(1,&tip_id,&indx,num,ppts) ;
    delete [] ppts ;
}

// ----------------------------------------------------------


void CArbMshCrackRegion2D::GeneralFlawGrowth(
                         const int num_tips,
                         const int *tip_ids,
                         const int *tip_indices,
                         const int num_points,
                         const CrackPt *const pts)
{
    if (DebugDumpFlag) {
        fprintf(dfd,"general_flaw_growth: %d %d [",num_tips,num_points) ;
        for (int di0=0 ; di0<num_points ; ++di0)
            fprintf(dfd," %d",tip_ids[di0]) ;
        fprintf(dfd," ] [") ;
        for (int di1=0 ; di1<num_points ; ++di1)
            fprintf(dfd," %d",tip_indices[di1]) ;
        fprintf(dfd," ]") ;
        for (int di2=0 ; di2<num_points ; ++di2)
            fprintf(dfd," [ %18.10f %18.10f %18.10f %d ]",
                pts[di2].coord[0],pts[di2].coord[0],pts[di2].char_elem_size,
                (pts[di2].has_char_size ? 1 : 0)) ;
        fprintf(dfd,"\n") ;
        fflush(dfd) ;
        if (DebugDumpOnly) return ;
    }

    fprintf(stderr,"Crack Growth Stub\n") ;
    FirstNewElem = MaxElemId+1 ;
}

// %(CArbMshCrackRegion2D::NumNodes-int-|^const)
/* ++ ----------------------------------------------------------
**
**    NumNodes - return the number of nodes in the mesh
**
**      int NumNodes() const
**
**      Description: This function returns the current number of nodes
**          in a mesh (typically called after crack insertion or
**          growth).
**
**      Return Value: number of nodes
**
**
** -- */

int CArbMshCrackRegion2D::NumNodes() const
{
    return(NodeTable->NumEntries()) ;
}




// %(CArbMshCrackRegion2D::NumElements-int-|^const)
/* ++ ----------------------------------------------------------
**
**    NumElements - return the number of elements in the mesh
**
**      int NumElements() const
**
**      Description: This function returns the current number of
**          elements in a mesh (typically called after crack insertion
**          or growth).
**
**      Return Value: number of elements
**
**
** -- */

int CArbMshCrackRegion2D::NumElements() const
{
    return(ElemTable->NumEntries()) ;
}




// %(CArbMshCrackRegion2D::NumCrackTips-int-|^const)
/* ++ ----------------------------------------------------------
**
**    NumCrackTips - return the number of crack-tips in a mesh
**
**      int NumCrackTips() const
**
**      Description: This function returns the current number
**          crack-tips defined in a mesh.
**
**      Return Value: number of crack tips.
**
**
** -- */

int CArbMshCrackRegion2D::NumCrackTips() const
{
    return(CrackTable->NumEntries()) ;
}




// %(CArbMshCrackRegion2D::NumTemplateElems-int-|^const)
/* ++ ----------------------------------------------------------
**
**    NumTemplateElems - returns the number template elements
**
**      int NumTemplateElems() const
**
**      Description: This function returns the current total number of
**          elements in all crack-tip templates.
**
**      Return Value: number of crack-tip template elements
**
**
** -- */

int CArbMshCrackRegion2D::NumTemplateElems() const
{
    return(TemplateElems->NumEntries()) ;
}




// %(CArbMshCrackRegion2D::GetNodes-ArbMshNode-|*^const)
/* ++ ----------------------------------------------------------
**
**    GetNodes - return a list of all mesh nodes
**
**      ArbMshNode *GetNodes() const
**
**      Description: This function returns a list of all the nodes
**          current defined for the mesh (typically called after crack
**          insertion or growth).
**
**      Return Value: A list of all mesh nodes. Ownership of this
**          memory passes to the client, which must eventually call
**          delete [].
**
**
** -- */

ArbMshNode *CArbMshCrackRegion2D::GetNodes() const
{
    int i ;
    ArbMshNode   *local_entries ;

    CArbHashTableIterator<int,CArbCoord2D> iter(NodeTable) ;
    local_entries = new ArbMshNode[NodeTable->NumEntries()] ;

    for (i=0,iter.First() ; iter.More() ; ++i,++iter) {
        local_entries[i].id = iter.Key() ;
        local_entries[i].coord[0] = iter.Entry()->x() ;
        local_entries[i].coord[1] = iter.Entry()->y() ;
    }

    return(local_entries) ;
}




// %(CArbMshCrackRegion2D::GetElements-ArbMshElement2D-|*^const)
/* ++ ----------------------------------------------------------
**
**    GetElements - return a list of all mesh elements
**
**      ArbMshElement2D *GetElements() const
**
**      Description: This function returns a list of all the elements
**          current defined for the mesh (typically called after crack
**          insertion or growth).
**
**      Return Value: A list of all mesh elements. Ownership of this
**          memory passes to the client, which must eventually call
**          delete [].
**
**
** -- */

ArbMshElement2D *CArbMshCrackRegion2D::GetElements() const
{
    CArbHashTableIterator<int,ArbMshElement2D> iter(ElemTable) ;
    ArbMshElement2D *local_elems =
        new ArbMshElement2D[ElemTable->NumEntries()] ;

    for (int i=0 ; iter.More() ; ++i,++iter)
        local_elems[i] = *(iter.Entry()) ;

    return(local_elems) ;
}




// %(CArbMshCrackRegion2D::GetCrackTipIds-int-|*^const)
/* ++ ----------------------------------------------------------
**
**    GetCrackTipIds - return a list of all crack-tip id's
**
**      int *GetCrackTipIds() const
**
**      Description: This function returns a list of all the currently
**          defined crack-tip id's.
**
**      Return Value: A list of all crack-tip id's. Ownership of this
**          memory passes to the client, which must eventually call
**          delete [].
**
**
** -- */

int *CArbMshCrackRegion2D::GetCrackTipIds() const
{
    int num = CrackTable->NumEntries() ;
    if (num == 0) return(0) ;
    return(CrackTable->GetKeyList()) ;
}




// %(CArbMshCrackRegion2D::GetCrackTipNodes-int-|*^const)
/* ++ ----------------------------------------------------------
**
**    GetCrackTipNodes - return an array of the crack tip node id's
**
**      int *GetCrackTipNodes() const
**
**      Description: This function returns a list of crack-tip node
**          id's These will be in the same order as the array' returned
**          from GetCrackTipIds. Only one crack-tip node id is returned
**          for each crack tip, even if the crack tip nodes are not
**          constrained. Ownership of this memory passes to the client,
**          which should eventually call delete [].
**
**      Return Value: A list of crack-tip nodes. Ownership of this
**          memory passes to the client, which should eventually call
**          delete [].
**
**
** -- */

int *CArbMshCrackRegion2D::GetCrackTipNodes() const
{
    int num = CrackTable->NumEntries() ;
    if (num == 0) return(0) ;

    CArbHashTableIterator<int,int> iter(CrackTable) ;

    int *local_entries = new int[num] ;

    for (int i=0 ; iter.More() ; ++i,++iter)
        local_entries[i] = *(iter.Entry()) ;

    return(local_entries) ;
}




// %(CArbMshCrackRegion2D::GetCrackTipNodeList-int-|*-int-const|-int-|*)
/* ++ ----------------------------------------------------------
**
**    GetCrackTipNodeList - returns a list of the nodes at a crack tip
**
**      int *GetCrackTipNodeList(
**              const int crack_id,
**              int       *num_nodes)
**
**        crack_id  - (in)  crack-tip id
**        num_nodes - (out) number of node id's in the list
**
**      Description: This function returns a list of the crack-tip
**          nodes for the given crack tip. If this is a constrained
**          crack-tip node there will only be one element in the list.
**          There will be more than one node for unconstrained crack
**          tip nodes.
**
**      Return Value: A list of crack-tip nodes. Ownership of this
**          memory passes to the client, which should eventually call
**          delete [].
**
**      Exceptions:
**          CArbMshCInvalidCrack - invalid crack id
**
** -- */

int *CArbMshCrackRegion2D::GetCrackTipNodeList(
                       const int crack_tip_id,
                       int *num_nodes) const
{
    // get the crack tip

    int *tip = CrackTable->Fetch(crack_tip_id) ;
    if (tip == 0) return(0) ;

    return(FindTipNodeList(*tip,num_nodes,MshTopo)) ;
}




// %(CArbMshCrackRegion2D::GetTemplateElemList-int-|*^const)
/* ++ ----------------------------------------------------------
**
**    GetTemplateElemList - return a list of crack-tip template element
**                          id's
**
**      int *GetTemplateElemList() const
**
**      Description: This function returns a list of the element id's
**          for all elements in crack-tip templates.
**
**      Return Value: A list of crack-tip template elements. Ownership
**          of this memory passes to the client, which should
**          eventually call delete [].
**
**
** -- */

int *CArbMshCrackRegion2D::GetTemplateElemList() const
{
    int num = TemplateElems->NumEntries() ;
    if (num == 0) return(0) ;

    CArbHashTableIterator<int,int> iter(TemplateElems) ;

    int *local_entries = new int[num] ;

    for (int i=0 ; iter.More() ; ++i,++iter)
        local_entries[i] = iter.Key() ;

    return(local_entries) ;
}


void CArbMshCrackRegion2D::GetCrackMouthNodes(int *nd0,int *nd1) const
{
    *nd0 = CrackMouthNodes[0] ;
    *nd1 = CrackMouthNodes[1] ;
}


// %(CArbMshCrackRegion2D::GetTemplateElemIds-int-|*^const)
/* ++ ----------------------------------------------------------
**
**    GetTemplateElemIds - return a list of template element crack id's
**
**      int *GetTemplateElemIds() const
**
**      Description: This function returns an array of crack id's
**          associated with template elements. There is a one-to-one
**          corespondance with entries in the list returned by
**          GetTemplateElemList.
**
**      Return Value: A list of the crack id's associated with
**          crack-tip template elements. Ownership of this memory
**          passes to the client, which should eventually call delete
**          [].
**
**
** -- */

int *CArbMshCrackRegion2D::GetTemplateElemIds() const
{
    int num = TemplateElems->NumEntries() ;
    if (num == 0) return(0) ;

    CArbHashTableIterator<int,int> iter(TemplateElems) ;

    int *local_entries = new int[num] ;

    for (int i=0 ; iter.More() ; ++i,++iter)
        local_entries[i] = *(iter.Entry()) ;

    return(local_entries) ;
}




// %(CArbMshCrackRegion2D::GetDeletedElements-int-|*^const)
/* ++ ----------------------------------------------------------
**
**    GetDeletedElements - return a list of deleted elements
**
**      int *GetDeletedElements() const
**
**      Description: This function returns a list of the id's of
**          elements deleted during the most recent remeshing.
**
**      Return Value: A list of the id's of deleted elements. Ownership
**          of this memory passes to the client, which should
**          eventually call delete [].
**
**
** -- */

int *CArbMshCrackRegion2D::GetDeletedElements() const
{
    if (DeletedElems == 0) return(0) ;
    return(DeletedElems->AsVector()) ;
}




// %(CArbMshCrackRegion2D::NumDeletedElements-int-|^const)
/* ++ ----------------------------------------------------------
**
**    NumDeletedElements - return the number of elements deleted during
**                         remeshing
**
**      int NumDeletedElements() const
**
**      Description: This function returns the number of elements
**          deleted during the most recent remeshing.
**
**      Return Value: number of deleted elements
**
**
** -- */

int CArbMshCrackRegion2D::NumDeletedElements() const
{
    return(DeletedElems->NumEntries()) ;
}




// %(CArbMshCrackRegion2D::GetCreatedElements-int-|*^const)
/* ++ ----------------------------------------------------------
**
**    GetCreatedElements - return a list of created elements
**
**      int *GetCreatedElements() const
**
**      Description: This function returns a list of the id's of
**          elements created during the most recent remeshing.
**
**      Return Value: A list of the id's of created elements. Ownership
**          of this memory passes to the client, which should
**          eventually call delete [].
**
**
** -- */

int *CArbMshCrackRegion2D::GetCreatedElements() const
{
    if (CreatedElems == 0) return(0) ;
    return(CreatedElems->AsVector()) ;
}




// %(CArbMshCrackRegion2D::NumCreatedElements-int-|^const)
/* ++ ----------------------------------------------------------
**
**    NumCreatedElements - return the number of elements created during
**                         remeshing
**
**      int NumCreatedElements() const
**
**      Description: This function returns the number of elements
**          created during the most recent remeshing.
**
**      Return Value: number of created elements
**
**
** -- */

int CArbMshCrackRegion2D::NumCreatedElements() const
{
    return(CreatedElems->NumEntries()) ;
}




// %(CArbMshCrackRegion2D::GetUpdatedBoundaryData-void-|-int-|*-BoundaryData-|**-int-|*-BoundaryData-|**)
/* ++ ----------------------------------------------------------
**
**    GetUpdatedBoundaryData - return information about updated
**                             boundaries
**
**      void GetUpdatedBoundaryData(
**              int          *num_old,
**              BoundaryData **old_data,
**              int          *num_new,
**              BoundaryData **new_data)
**
**        num_old  - (out) number of items in the old boundary list
**        old_data - (out) list of boundary segments that existed in
**                         the old mesh and have been deleted or
**                         modified in the new mesh.
**        num_new  - (out) number of items in the new boundary list
**        new_data - (out) list of the boundary segments that exist in
**                         the new mesh, that do not correspond to
**                         segments in the old mesh.
**
**      Description: This function returns information about boundaries
**          that were updated during remeshing. Ownership of this
**          old_data and new_data lists passes to the client, which
**          should eventually call delete [] for both.
**
**
** -- */

void CArbMshCrackRegion2D::GetUpdatedBoundaryData(
              int *num_old,BoundaryData **old_data,
              int *num_new,BoundaryData **new_data)
{
    if (OldBdryTable != 0) {
        *num_old = OldBdryTable->NumEntries() ;
        CArbHashTableIterator<ArbEdgeKey,BoundaryData>
                iter(OldBdryTable) ;
        *old_data = new BoundaryData[*num_old] ;
        for (int i=0 ; iter.More() ; ++i,++iter)
            (*old_data)[i] = *(iter.Entry()) ;
    } else {
        *num_old = 0 ;
        *old_data = 0 ;
    }

    if (NewBdryTable != 0) {
        *num_new = NewBdryTable->NumEntries() ;
        CArbHashTableIterator<ArbEdgeKey,BoundaryData>
                iter(NewBdryTable) ;
        *new_data = new BoundaryData[*num_new] ;
        for (int i=0 ; i<*num_new ; ++i,++iter)
            (*new_data)[i] = *(iter.Entry()) ;
    } else {
        *num_new = 0 ;
        *new_data = 0 ;
    }
}




// %(CArbMshCrackRegion2D::FindTipNodeList-int-|*^const-int-const|-int-|*-CArbMshTopo2D-|*)
/* ++ ----------------------------------------------------------
**
**    FindTipNodeList - return a list of crack-tip node id's
**
**      int *FindTipNodeList(
**              const int     tip_node_id,
**              int           *num_nodes,
**              CArbMshTopo2D *topo) const
**
**        tip_node_id - (in)  crack-tip id
**        num_nodes   - (out) number of crack-tip nodes
**        topo        - (in)  mesh topology object
**
**      Description: This function returns an array of the crack tip
**          nodes for the given crack tip. If this is a constrained
**          crack-tip node there will only be one element in the list.
**          There will be more than one node for unconstrained crack
**          tip nodes.
**
**      Return Value: A list of crack-tip node id's. Ownership of the
**          list passes to the caller.
**
**
** -- */

int *CArbMshCrackRegion2D::FindTipNodeList(
                       const int tip_node_id,
                       int *num_nodes,
                       CArbMshTopo2D *topo) const
{
    int loc_tip = tip_node_id ;

    // get the crack tip coordinates

    CArbCoord2D *tipc = NodeTable->Fetch(loc_tip) ;
    int num = 1 ;

    // first loop CW from this node finding all nodes with
    // the same coordinates

    CArbTopoAdjVtxCyclicIterator citer(topo,loc_tip) ;
    while(1) {
        while (citer.CcwElem() != NO_ELEM) ++citer ;
        ++citer ;
        CArbCoord2D *next = NodeTable->Fetch(citer.AdjVtx()) ;
        if (*next == *tipc)
            loc_tip = citer.AdjVtx() ;
        else
            break ;
        citer.NewVtx(citer.AdjVtx()) ;
    }


    // Now loop CCW from this node finding all nodes with
    // the same coordinates

    CArbTopoAdjVtxIterator iter(topo,loc_tip) ;

    while(1) {
        while (iter.CcwElem() != NO_ELEM) ++iter ;
        CArbCoord2D *next = NodeTable->Fetch(iter.AdjVtx()) ;
        if (*next == *tipc)
            ++num ;
        else
            break ;
        iter.NewVtx(iter.AdjVtx()) ;
    }

    // allocate memory

    int *nodes = new int[num] ;
    nodes[0] = loc_tip ;
    int cur = 1 ;

    // Now loop CCW from this node finding all nodes with
    // the same coordinates

    iter.NewVtx(loc_tip) ;
    while(cur < num) {
        while (iter.CcwElem() != NO_ELEM) ++iter ;
        nodes[cur] = iter.AdjVtx() ;
        ++cur ;
        iter.NewVtx(iter.AdjVtx()) ;
    }
    *num_nodes = num ;
    return(nodes) ;
}




// %(CArbMshCrackRegion2D::GetTipData-CrackTipData-|*-int-|)
/* ++ ----------------------------------------------------------
**
**    GetTipData - return crack-tip data
**
**      CrackTipData *GetTipData(int crack_tip_id)
**
**        crack_tip_id - (in)  crack-tip id
**
**      Description: This function returns the current crack-tip
**          parameters for the given crack-tip id.
**
**      Return Value: The crack-tip parameters for this crack.
**          Ownership passes to the caller.
**
**      Exceptions:
**          CArbMshCInvalidCrack - invalid crack id
**
** -- */

CArbMshCrackRegion2D::CrackTipData
*CArbMshCrackRegion2D::GetTipData(int crack_tip_id)
{
    if (TipData == 0) {
        TipData = new CArbHashTable<int,CrackTipData>() ;
    }
    CrackTipData *data = TipData->Fetch(crack_tip_id) ;
    if (data == 0) {
        CrackTipData ddata = { 0, false, S_TRIANGLE,
                                  false, T_TRIANGLE,
                                  false, true,
                                  false, true,
                                  false, 8, true,
                                  false, 0.0,
                                  false, 1.0,
                                  false, 2 } ;
        ddata.tip_id = crack_tip_id ;
        TipData->Store(crack_tip_id,ddata) ;
        data = TipData->Fetch(crack_tip_id) ;
    }
    return(data) ;
}




// %(CArbMshCrackRegion2D::NewNode-ArbMshNode-|-double-|-double-|)
/* ++ ----------------------------------------------------------
**
**    NewNode - generate a new node
**
**      ArbMshNode NewNode(
**              double x,
**              double y)
**
**        x - (in)  node's x coordinate
**        y - (in)  node's y coordinate
**
**      Description: This function generates a new node and assigns it
**          the next available node id.
**
**      Return Value: the new node object
**
**
** -- */

ArbMshNode CArbMshCrackRegion2D::NewNode(double x,double y)
{
    ArbMshNode node ;
    ++MaxNodeId ;
    node.id = MaxNodeId ;
    node.coord[0] = x ;
    node.coord[1] = y ;
    node.coord[2] = 0.0 ;
    return(node) ;
}






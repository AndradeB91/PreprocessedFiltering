//
// ArbCrackRegion2D header file
//
// Description -
//   This is the header file for the ArbMshCrackRegion2D objects.
//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 2000
//   All rights reserved
//
// Author -
//   Wash Wawrzynek
//
// Revision -
//   $Revision: 1.22 $  $Date: 2002/07/16 14:19:04 $  $Author: wash $
//
 
#ifndef ArbMshCrackRegion2D_h
#define ArbMshCrackRegion2D_h

#define DEFAULT_ID 100000000

#include "ArbMsh.hpp"
#include "ArbCoord2D.hpp"
#include "ArbHashTable.hpp"
#include "ArbSet.hpp"
#include "ArbMshTopo2D.hpp"
#include "ArbMshRegion2D.hpp"
#include "ArbBSpline2D.hpp"

// Exceptions thrown by this object

#ifndef ARB_NORMAL_STATUS
#define ARB_NORMAL_STATUS      0
#endif
#ifndef ARB_DUPLICATE_NODE_ID
#define ARB_DUPLICATE_NODE_ID  1
#endif
#ifndef ARB_DUPLICATE_ELEM_ID
#define ARB_DUPLICATE_ELEM_ID  2
#endif
#ifndef ARB_DUPLICATE_CRACK_ID
#define ARB_DUPLICATE_CRACK_ID 3
#endif
#ifndef ARB_INVALID_ELEM
#define ARB_INVALID_ELEM       4
#endif
#ifndef ARB_INVALID_NODE
#define ARB_INVALID_NODE       5
#endif
#ifndef ARB_INVALID_CRACK
#define ARB_INVALID_CRACK      6
#endif
#ifndef ARB_BAD_NUMBER
#define ARB_BAD_NUMBER         7
#endif

typedef CArbHashTable<int,int> MaterialHash ;
typedef CArbSet<int> CornerNodes ;

class CArbMshCrackRegion2D {

    public:

    // object data types

        enum ElemShape {S_TRIANGLE, S_QUADRILATERAL} ;

        enum ElemType  {T_TRIANGLE, T_COLLAPSED_QUAD} ;

        enum SegEndType {S_MINIMUM, S_AVERAGE, S_ANISOTROPIC,
                         S_VARIABLE, S_USER_ANISOTROPIC} ;

        enum AngleCheckLocation { AT_NODES, AT_INT_PTS } ;

        struct CrackPt {
            double coord[2] ;
            double char_elem_size ;
            bool has_char_size ;
        } ;

        struct BoundaryData {
            int elem_id ;
            int node_id[3] ;
            bool bi_mat_flag ;
        } ;

    // constructors and destructors

        CArbMshCrackRegion2D() ;

        ~CArbMshCrackRegion2D() ;

    // routines to add nodes and elements to the mesh and 
    // to define crack tips

        int AddNode(const ArbMshNode &inode) ;

        int AddNodeList(const int num_nodes,
                        const ArbMshNode *const node_list) ;

        int AddElem(const ArbMshElement2D &elem) ;

        int AddElemList(const int num_elems,
                        const ArbMshElement2D *const elem_list) ;

        int DefineCrackTip(const int crack_tip_id,
                           const int node_id) ;

    // set and query parameters effecting crack growth

        int SetTipTemplate(
                 const ElemShape ishape,
                 const int crack_tip_id = DEFAULT_ID) ;
        int SetTipElemType(
                 const ElemType itype,
                 const int crack_tip_id = DEFAULT_ID);
        int SetQrtrPointTip(
                  const bool flag = true,
                  const int crack_tip_id = DEFAULT_ID) ;
        int SetConstrainedTip(
                  const bool flag = true,
                  const int crack_tip_id = DEFAULT_ID) ;
        int SetNumTipElems(
                 const int num_elem,
                 const int crack_tip_id = DEFAULT_ID) ;
        int SetTemplateRadius(
                  const double radius,
                  const int crack_tip_id = DEFAULT_ID) ;
        int SetTemplateRatio(
                  const double ratio,
                  const int crack_tip_id = DEFAULT_ID) ;
        int SetTemplateNumRings(
                  const int num_rings,
                  const int crack_tip_id = DEFAULT_ID) ;
        void SetRemeshElemShape(const ElemShape ishape) ;
        void SetNoRemeshMat(const int mat_id) ;
        void SetTemplateElem(const int elem_id,const int tip_id) ;

        void SetTriDeleteFactor(double f)  { TriDelFactor = f ; } ;
        void SetQuadDeleteFactor(double f) { QuadDelFactor = f ; } ;

        int SetStartNodeID(int id)
            { if (id > MaxNodeId) {
                  MaxNodeId = id-1 ;
                  return(ARB_NORMAL_STATUS) ;
              } else return(ARB_BAD_NUMBER) ;
            } ;
        int SetStartElemID(int id)
            { if (id > MaxElemId) {
                  MaxElemId = id-1 ;
                  return(ARB_NORMAL_STATUS) ;
              } else return(ARB_BAD_NUMBER) ;
            } ;

        void SetQuadElemAngleChecks(double min_angle,
                                    double max_angle,
                                    AngleCheckLocation location) ;

        ElemShape GetTipTemplate(
                 const int crack_tip_id = DEFAULT_ID) const ;
        ElemType GetTipElemType(
                 const int crack_tip_id = DEFAULT_ID) const ;
        bool GetQrtrPointTip(
                 const int crack_tip_id = DEFAULT_ID) const ;
        bool GetConstrainedTip(
                 const int crack_tip_id = DEFAULT_ID) const ;
        int GetNumTipElems(
                 const int crack_tip_id = DEFAULT_ID) const ;
        double GetTemplateRadius(
                 const int crack_tip_id = DEFAULT_ID) const ;
        double GetTemplateRatio(
                 const int crack_tip_id = DEFAULT_ID) const ;
        int GetTemplateNumRings(
                 const int crack_tip_id = DEFAULT_ID) const ;
        int *GetNoRemeshMats(int *num) const ;
        ElemShape GetRemeshElemShape() const ;
        ArbMshOrder GetElementOrder() const ;

    // define new cracks

        int NewInternalCrack(const int start_tip_id,
                             const int end_tip_id,
                             const int num_points,
                             const CrackPt *const pts) ;

        int NewSurfaceCrack(const int tip_id,
                            const int num_points,
                            const CrackPt *const pts) ;

        void NewGeneralFlaw(const int num_tips,
                            const int *tip_ids,
                            const int *tip_indices,
                            const int num_points,
                            const CrackPt *const pts) ;

    // specify crack growth

        void CrackGrowth(const int crack_tip_id,
                         const int num_points,
                         const CrackPt *const pts) ;

        void GeneralFlawGrowth(const int num_tips,
                               const int *tip_ids,
                               const int *tip_indices,
                               const int num_points,
                               const CrackPt *const pts) ;

    // routines to query info about the new mesh

        int NumNodes() const ;
        int NumElements() const ;
        int NumCrackTips() const ;
        int NumTemplateElems() const ;

        ArbMshNode *GetNodes() const ;
        ArbMshElement2D *GetElements() const ;
        int *GetCrackTipIds() const ;
        int *GetCrackTipNodes() const ;
        int *GetCrackTipNodeList(
                         const int crack_id,
                         int *num_nodes) const ;
        int *GetTemplateElemList() const ;
        int *GetTemplateElemIds() const ;

        void GetCrackMouthNodes(int *nd0,int *nd1) const ;

        int NumDeletedElements() const ;
        int *GetDeletedElements() const ;

        int NumCreatedElements() const ;
        int *GetCreatedElements() const ;

        void GetUpdatedBoundaryData(int *num_old,
                                    BoundaryData **old_data,
                                    int *num_new,
                                    BoundaryData **new_data) ;

        void SetAnisotropicRatio(double val)
            { AnisotropicRatio = val ;
              SegEndCharSizeType = S_USER_ANISOTROPIC ; } ;

        // routines to set tuning parameters

        void SetBoundaryTolerance(double tol)
            { BoundaryTolerance = tol ; } ;
        void SetTemplateBoundaryFactor(double fact)
            { TemplateBoundaryFactor = fact ; } ;
        void SetSegEndCharSizeFactor(double fact)
            { SegEndCharSizeFactor = fact ; } ;
        void SetSegEndCharSizeType(SegEndType type)
            { SegEndCharSizeType = type ; } ;
        void SetBoundarySpacingFactor(double fact)
            { BoundarySpacingFactor = fact ; } ;

        // debug display options

        enum DDFlags { CrackAll=63 } ;

        int DebugDisplayFlags ;

        void SetDebugDisplayFlags(DDFlags flags) {
            DebugDisplayFlags |= flags ; } ;

    private:

#ifdef WIN32
    public:
#endif
        struct CrackTipData {
            int tip_id ;
            bool         elem_shape_set ;
            ElemShape    tip_elem_shape ;
            bool         elem_type_set ;
            ElemType     tip_elem_type ;
            bool         qrtr_pt_set ;
            bool         qrtr_pt_flag ;
            bool         const_tip_set ;
            bool         const_tip_flag ;
            bool         num_tip_set ;
            int num_tip_elems ;
            bool         default_num_tip ;
            bool         temp_radius_set ;
            double       template_radius ;
            bool         prog_ratio_set ;
            double       progression_ratio ;
            bool         num_rings_set ;
            int number_of_rings ;
        } ;
#ifdef WIN32
    private:
#endif

        CArbHashTable<int,CArbCoord2D> *NodeTable ;
        CArbHashTable<int,ArbMshElement2D> *ElemTable ;
        CArbHashTable<int,int> *CrackTable ;
        CArbHashTable<int,int> *RetainedNodes ;

        CArbMshTopo2D *MshTopo ;

        ElemShape    ElemShapeType ;
        int MaxNodeId ;
        int MaxElemId ;
        bool         FirstElem ;
        ArbMshOrder  Order ;
        bool         CanResizeTemplate ;
        int FirstNewElem ;

        CArbSet<int> *NoRmshMats ;
        CArbHashTable<int,int> *TemplateElems ;
        
        CrackTipData DefaultTipData ;
        CArbHashTable<int,CrackTipData> *TipData ;
        CArbArray<int> *DeletedElems ;
        CArbArray<int> *CreatedElems ;

        CArbHashTable<ArbEdgeKey,BoundaryData> *OldBdryTable ;
        CArbHashTable<ArbEdgeKey,BoundaryData> *NewBdryTable ;

        double TriDelFactor ;
        double QuadDelFactor ;

        int CrackMouthNodes[2] ;

        double AnisotropicRatio ;

        bool               DoQuadAngleChecks ;
        double             MinQuadAngle ;
        double             MaxQuadAngle ;
        AngleCheckLocation QuadAngleLocation ;

        // tuning parameters

        double BoundaryTolerance ;
        double TemplateBoundaryFactor ;
        double SegEndCharSizeFactor ;
        SegEndType SegEndCharSizeType ;
        double BoundarySpacingFactor ;

        bool DebugDumpFlag ;

        int *FindTipNodeList(
                       const int tip_node_id,
                       int *num_nodes,
                       CArbMshTopo2D *topo) const ;

        CrackTipData *GetTipData(int crack_tip_id) ;

        ArbMshNode NewNode(double x,double y) ;

        int NewElemId() { return(++MaxElemId) ; } ;

        int NewNodeId() { return(++MaxNodeId) ; } ;

        friend class CArbHashTable<int,CrackTipData> ;
        friend int operator == (const CrackTipData &t1,
                                const CrackTipData &t2) ;

        friend class CArbRmshRegion2D ;
} ;

inline int operator == (const CArbMshCrackRegion2D::CrackTipData &t1,
                        const CArbMshCrackRegion2D::CrackTipData &t2)
{
    return(t1.tip_id == t2.tip_id) ;
}

inline int operator == (const CArbMshCrackRegion2D::BoundaryData &b1,
                        const CArbMshCrackRegion2D::BoundaryData &b2)
{
    return(b1 == b2) ;
}


/*
TYPEDEFS

CArbHashTable<int,int> MaterialHash - Hash table that maps region id's 
            to material id's 

CArbSet<int> CornerNodes - Set containing the id's of corner nodes 



CLASS CArbMshCrackRegion2D

  This object controls remeshing of a region for crack insertion or 
  growth. 

  The normal procedure for using this object is as follows: 

  1. Create an instance of the object. 

  2. Add nodes to the region with AddNode and/or AddNodeList. 

  3. Add elements to the region with AddElem and/or AddElemList. 

  4. If necessary, specify which nodes are crack-tip nodes with 
     DefineCrackTip. 

  5. Modify the mesh by: 

     1) Inserting a new flaw with NewZeroVolIntFlaw, 
        NewZeroVolSurfFlaw,NewFiniteVolIntFlaw, or NewFiniteVolSurfFlaw 

     2) Growing an existing flaw with ZeroVolFlawGrowth. 

  6. GetNodes and GetElements are called to retieve information about 
     the updated mesh. 


PUBLIC INTERFACE

  Public Data Structures:

    struct CrackPt

      This structure is used to describe the points that define the 
      crack. A crack is defined by an array of these points. 

      Member Variables:

        double coord[2] - the point's coordinates 

        double char_elem_size - desired characteristic element 
            size at this point 

        bool has_char_size - true if a characteristic element 
            size is specified 


    struct BoundaryData

      This data structure is used to store and communicate information 
      about a boundary segment that has been created, deleted, or 
      modified during remeshing. 

      Member Variables:

        int elem_id - id of the element adjacent to the segment 

        int node_id[3] - id's of the nodes associated with the 
            segment 

        bool bi_mat_flag - true indicates that his segment is 
            part of a bi-material interface 


    struct CrackTipData

      This data structure is used to store crack-tip parameters. If a 
      flag for a specific parameter is not set, a default value is 
      used. 

      Member Variables:

        int tip_id - the crack-tip's id 

        bool elem_shape_set - true indicates that the 
            tip_elem_shape parameter is set 

        ElemShape tip_elem_shape - either S_TRIANGLE or 
            S_QUADRILATERAL 

        bool elem_type_set - true indicates that the 
            tip_elem_type parameter is set 

        ElemType tip_elem_type - either T_TRIANGLE or 
            T_COLLAPSED_QUAD 

        bool qrtr_pt_set - true indicates that the qrtr_pt_flag 
            parameter is set 

        bool qrtr_pt_flag - true indicates that quarter-point 
            elements should be used for quadratic order 
            crack-tip elements 

        bool const_tip_set - true indicates that the 
            const_tip_flag parameter is set 

        bool const_tip_flag - true indicates that collapsed 
            quadrilateral crack-tip elements should have the 
            crack-tip nodes constrained to move together 

        bool num_tip_set - true indicates that the 
            num_tip_elems parameter is set 

        int num_tip_elems - number of elements to place at the 
            crack tip 

        bool default_num_tip - default number of elements to 
            place at the crack tip 

        bool temp_radius_set - true indicates that the 
            template_radius parameter is set 

        double template_radius - template radius to use at the 
            crack-tip 

        bool prog_ratio_set - true indicates that the 
            progression_ratio parameter is set 

        double progression_ratio - element size progresssion 
            ratio to use in the crack-tip template 

        bool num_rings_set - true indicates that the 
            number_of_rings parameter is set 

        int number_of_rings - number of element rings to use in 
            the crack-tip template 


  Public Member Functions:

    CArbMshCrackRegion2D - constructor 

      CArbMshCrackRegion2D()

      Description: This is the constructor for a ArbMshCrackRegion2D 
          object. 


    CArbMshCrackRegion2D - destructor 

      ~CArbMshCrackRegion2D()

      Description: This is the destructor for a ArbMshCrackRegion2D 
          object. 


    AddNode - add a node to the mesh 

      void AddNode(const ArbMshNode &inode)

        inode - (in)  description of the node 

      Description: This function adds a new node to the mesh. 

      Exceptions:
          CArbMshCDuplNode - duplicate node id's


    AddNodeList - add an array of nodes to the mesh 

      void AddNodeList(
              const int        num_nodes,
              const ArbMshNode *constnode_list)

        num_nodes - (in)  number of nodes 
        node_list - (in)  list of nodes to add 

      Description: This function adds an array of nodes to the mesh. 

      Exceptions:
          CArbMshCDuplNode - duplicate node id's


    AddElem - add an element to the mesh 

      void AddElem(const ArbMshElement2D &elem)

        elem - (in)  element to add 

      Description: This function adds an element to the mesh. 

      Exceptions:
          CArbMshCDuplElem - duplicate element id's
          CArbMshCInvalidElem - invalid element


    AddElemList - add an array of elements to the mesh 

      void AddElemList(
              const int             num_elems,
              const ArbMshElement2D *constelem_list)

        num_elems - (in)  number of elements 
        elem_list - (in)  elements 

      Description: This function adds an array of elements to the 
          mesh. 

      Exceptions:
          CArbMshCDuplElem - duplicate element id's
          CArbMshCInvalidElem - invalid element


    DefineCrackTip - define a node as a crack tip node 

      void DefineCrackTip(
              const int crack_tip_id,
              const int node_id)

        crack_tip_id - (in)  crack tip id 
        node_id      - (in)  crack-tip node id 

      Description: This function defines a node as a crack tip and 
          associates a crack-tip id with the node. 

      Exceptions:
          CArbMshCDuplCrack - duplicate crack id's
          CArbMshCInvalidNode - invalid node id


    SetTipTemplate - specify crack-tip template element types 

      void SetTipTemplate(
              const ElemShape ishape,
              const int       crack_tip_id = DEFAULT_ID)

        ishape       - (in)  either S_TRIANGLE or S_QUADRILATERAL 
        crack_tip_id - (in)  id of the crack-tip, or DEFAULT_ID to 
                             indicate that this data should be used 
                             for all tips not specifically identified 
                             in other calls to this routine 

      Description: This function allow the type of elements to be 
          used at the crack-tip to be specified. 

      Exceptions:
          CArbMshCInvalidCrack - invalid crack tip id


    SetTipElemType - specify crack-tip element types 

      void SetTipElemType(
              const ElemType itype,
              const int      crack_tip_id = DEFAULT_ID)

        itype        - (in)  either T_TRIANGLE or T_COLLAPSED_QUAD 
        crack_tip_id - (in)  id of the crack-tip, or DEFAULT_ID to 
                             indicate that this data should be used 
                             for all tips not specifically identified 
                             in other calls to this routine. 

      Description: This function allow the type of elements to be 
          used at the crack-tip to be specified. 

      Exceptions:
          CArbMshCInvalidCrack - invalid crack tip id


    SetQrtrPointTip - specify quarter-point crack tip elements 

      void SetQrtrPointTip(
              const bool flag = true,
              const int  crack_tip_id = DEFAULT_ID)

        flag         - (in)  true means use quarter-point elements 
        crack_tip_id - (in)  id of the crack-tip, or DEFAULT_ID to 
                             indicate that this data should be used 
                             for all tips not specifically identified 
                             in other calls to this routine. 

      Description: This function sets a flag that controls if 
          crack-tip elements will be quarter-point elements in 
          quadratic order meshes. The default is to use quarter-point 
          elements. 

      Exceptions:
          CArbMshCInvalidCrack - invalid crack id


    SetConstrainedTip - specify constrained collapsed crack-tip 
                        elements 

      void SetConstrainedTip(
              const bool flag = true,
              const int  crack_tip_id = DEFAULT_ID)

        flag         - (in)  true means constrain the crack-tip 
                             elements 
        crack_tip_id - (in)  id of the crack-tip, or DEFAULT_ID to 
                             indicate that this data should be used 
                             for all tips not specifically identified 
                             in other calls to this routine. 

      Description: This function sets a flag that controls if 
          triangular crack-tip elements, which are modeled as 
          collapsed quadrilateral elements, have the nodes at the 
          crack-tip to be constrained to move together. The default 
          is to constrain the nodes. 

      Exceptions:
          CArbMshCInvalidCrack - invalid crack id


    SetNumTipElems - set the number of crack-tip elements 

      void SetNumTipElems(
              const int num_elem,
              const int crack_tip_id = DEFAULT_ID)

        num_elem     - (in)  number of crack tip elements 
        crack_tip_id - (in)  id of the crack-tip, or DEFAULT_ID to 
                             indicate that this data should be used 
                             for all tips not specifically identified 
                             in other calls to this routine. 

      Description: This function sets the number of elements to 
          insert at the crack-tip. The default is 8 for triangular 
          elements and 4 for quadrilateral elements. 

      Exceptions:
          CArbMshCInvalidCrack - invalid crack id


    SetTemplateRadius - set the crack-tip template radius 

      void SetTemplateRadius(
              const double radius,
              const int    crack_tip_id = DEFAULT_ID)

        radius       - (in)  template radius 
        crack_tip_id - (in)  id of the crack-tip, or DEFAULT_ID to 
                             indicate that this data should be used 
                             for all tips not specifically identified 
                             in other calls to this routine. 

      Description: This function sets the radius of the crack-tip 
          element template. 

      Exceptions:
          CArbMshCInvalidCrack - invalid crack id


    SetTemplateRatio - sets the crack-tip template progression ratio 

      void SetTemplateRatio(const double ratio = DEFAULT_ID)

        ratio - (in)  id of the crack-tip, or DEFAULT_ID to indicate 
                      that this data should be used for all tips not 
                      specifically identified in other calls to this 
                      routine. 

      Description: This function sets the progression ratio for the 
          crack-tip template elements. 

      Exceptions:
          CArbMshCInvalidCrack - invalid crack id


    SetTemplateNumRings - set the number of rings in the crack-tip 
                          template 

      void SetTemplateNumRings(
              const int num_rings,
              const int crack_tip_id = DEFAULT_ID)

        num_rings    - (in)  number of element rings 
        crack_tip_id - (in)  id of the crack-tip, or DEFAULT_ID to 
                             indicate that this data should be used 
                             for all tips not specifically identified 
                             in other calls to this routine. 

      Description: This function sets the number of element rings 
          placed in a crack-tip template. 

      Exceptions:
          CArbMshCInvalidCrack - invalid crack id


    SetRemeshElemShape - specify the type of elements used for 
                         remeshing 

      void SetRemeshElemShape(const ElemShape ishape)

        ishape - (in)  either S_TRIANGLE or S_QUADRILATERAL 

      Description: This function specifies the type (triangle or 
          quadrilateral) of elements used during remeshing. 


    SetNoRemeshMat - specify a "no remesh" material 

      void SetNoRemeshMat(const int mat_id)

        mat_id - (in)  id of a material that should not be modified 

      Description: This function sets a flag that indicates that no 
          elements with the specified material id will be deleted or 
          modified during crack insertion or growth. 


    SetTemplateElem - flag an element as a template element 

      void SetTemplateElem(
              const int elem_id,
              const int tip_id)

        elem_id - (in)  id of an element 
        tip_id  - (in)  id of a crack tip 

      Description: This function flags an element as a template 
          element associated with the specified crack tip. Template 
          elements are not deleted during remeshing of other crack 
          tips. 


    SetTriDeleteFactor - set the triangular delete factor 

      void SetTriDeleteFactor(double f)

        f - (in)  delete factor 

      Description: This function sets the delete factor to use with 
          remeshing using triangular elements. The delete factor is 
          used to determine which elements should be deleted. The 
          factor is multiplied by the the characteristic sizes 
          specified for the crack points to find delete radii. All 
          elements adjacent to all node that falls in the delete are 
          removed. 


    SetQuadDeleteFactor - set the quadrilateral delete factor 

      void SetQuadDeleteFactor(double f)

        f - (in)  delete factor 

      Description: This function sets the delete factor to use with 
          remeshing using quadrilateral elements. The delete factor 
          is used to determine which elements should be deleted. The 
          factor is multiplied by the the characteristic sizes 
          specified for the crack points to find delete radii. All 
          elements adjacent to all node that falls in the delete are 
          removed. 


    GetTipTemplate - return the template element shapes 

      ElemShape GetTipTemplate(const int crack_tip_id = DEFAULT_ID)

        crack_tip_id - (in)  crack-tip id or DEFAULT_ID 

      Description: This function returns the shape (quadrilateral or 
          triangle) of elements used at the crack-tip. 

      Return Value: either S_TRIANGLE or S_QUADRILATERAL 

      Exceptions:
          CArbMshCInvalidCrack - invalid crack-tip id


    GetTipElemType - return the template element types 

      ElemType GetTipElemType(const int crack_tip_id = DEFAULT_ID)

        crack_tip_id - (in)  crack-tip id or DEFAULT_ID 

      Description: This function returns the type (triangle or 
          collapsed quadrilateral) of elements used at the crack-tip. 

      Return Value: either T_TRIANGLE or T_COLLAPSED_QUAD 

      Exceptions:
          CArbMshCInvalidCrack - invalid crack-tip id


    GetQrtrPointTip - return the quarter-point tip flag 

      bool GetQrtrPointTip(const int crack_tip_id = DEFAULT_ID)

        crack_tip_id - (in)  crack-tip id or DEFAULT_ID 

      Description: This function returns the flag that indicatese if 
          quarter-point elements will be inserted at a crack tip for 
          quadratic order meshes. 

      Return Value: true means quarter-point elements will be used 

      Exceptions:
          CArbMshCInvalidCrack - invalid crack-tip id


    GetConstrainedTip - return the constrained tip flag 

      bool GetConstrainedTip(const int crack_tip_id = DEFAULT_ID)

        crack_tip_id - (in)  crack-tip id or DEFAULT_ID 

      Description: This function returns the flag that indicatese if 
          crack-tip nodes for collapsed quadrilateral crack-tip 
          elements will be constrained to move together. 

      Return Value: true means that the crack-tip elements will be 
          constrained to move together 

      Exceptions:
          CArbMshCInvalidCrack - invalid crack-tip id


    GetNumTipElems - return the number of crack-tip elements 

      int GetNumTipElems(const int crack_tip_id = DEFAULT_ID)

        crack_tip_id - (in)  crack-tip id or DEFAULT_ID 

      Description: This function returns the number of crack-tip 
          elements that will be inserted at a crack tip. 

      Return Value: number of elements 

      Exceptions:
          CArbMshCInvalidCrack - invalid crack-tip id


    GetTemplateRadius - return the crack-tip radius 

      double GetTemplateRadius(const int crack_tip_id = DEFAULT_ID)

        crack_tip_id - (in)  crack-tip id or DEFAULT_ID 

      Description: This function returns the radius of the crack-tip 
          template inserted at a crack tip. 

      Return Value: the template radius 

      Exceptions:
          CArbMshCInvalidCrack - invalid crack-tip id


    GetTemplateRatio - return the progression ratio 

      double GetTemplateRatio(const int crack_tip_id = DEFAULT_ID)

        crack_tip_id - (in)  crack-tip id or DEFAULT_ID 

      Description: This function returns the element size progression 
          ratio that will be used when inserting elements at a crack 
          tip. 

      Return Value: progression ratio 

      Exceptions:
          CArbMshCInvalidCrack - invalid crack-tip id


    GetTemplateNumRings - return the number of template element rings 

      int GetTemplateNumRings(const int crack_tip_id = DEFAULT_ID)

        crack_tip_id - (in)  crack-tip id or DEFAULT_ID 

      Description: This function returns the number of element rings 
          that will be used when generating a crack-tip template. 

      Return Value: number of element rings 

      Exceptions:
          CArbMshCInvalidCrack - invalid crack-tip id


    GetNoRemeshMats - return a list of the "no remesh" materials 

      int *GetNoRemeshMats(int *num) const

        num - (out) number of id's in the list 

      Description: This function returns a list of the id's of the 
          materials that will not be modified during remeshing. 

      Return Value: An array of element id's. Owner should of this 
          list passes to the client, which should eventially call 
          delete []. 


    GetRemeshElemShape - returns the element shape used during 
                         remeshing 

      ElemShape GetRemeshElemShape() const

      Description: This function returns the shape of the elements 
          (triangle or quadrilateral) that will be used during 
          remeshing. 

      Return Value: either S_TRIANGULAR or S_QUADRILATERAL 


    GetElementOrder - return the polynomial order of the mesh 
                      elements 

      ArbMshOrder GetElementOrder() const

      Description: This function returns the polynomial order of the 
          elements used in the mesh. 

      Return Value: either LINEAR or QUADRATIC 


    NewZeroVolIntFlaw - add a new zero volume internal flaw 

      void NewZeroVolIntFlaw(
              const int     start_tip_id,
              const int     end_tip_id,
              const int     num_points,
              const CrackPt *constpts)

        start_tip_id - (in)  id to be assigned to the first crack tip 
                             (corresponding to the first point in the 
                             list) 
        end_tip_id   - (in)  id to be assigned to the second crack 
                             tip (corresponding to the last point in 
                             the list) 
        num_points   - (in)  number of points in the crack-point list 
        pts          - (in)  list of crack-points, which define the 
                             new crack 

      Description: This function controls the process of adding a new 
          zero volume internal flaw to a mesh. 

      Exceptions:
          CArbMshCInvalidCrack - the specified crack is invalid


    NewZeroVolSurfFlaw - add a new zero volume surface flaw 

      void NewZeroVolSurfFlaw(
              const int     tip_id,
              const int     num_points,
              const CrackPt *constpts)

        tip_id     - (in)  id to be assinged to the new crack tip 
        num_points - (in)  number of points in the crack-point list 
        pts        - (in)  list of crack-points, which define the new 
                           crack 

      Description: This function controls the process of adding a new 
          zero volume surface flaw. It is assumed that the first 
          point in the list is on the surface of the mesh and the 
          ending point is the new crack tip. 

      Exceptions:
          CArbMshCInvalidCrack - the specified crack is invalid


    NewFiniteVolIntFlaw - add a new finite volume internal flaw 

      void NewFiniteVolIntFlaw(
              const int     num_tips,
              const int     *tip_ids,
              const int     *tip_indices,
              const int     num_points,
              const CrackPt *constpts)

        num_tips    - (in)  number of crack-tips in the flaw 
        tip_ids     - (in)  list of id's to assign to crack tips 
        tip_indices - (in)  location in the crack-point array of the 
                            nodes that will be crack tips 
        num_points  - (in)  number of points in the crack-point list 
        pts         - (in)  list of crack-points, which define the 
                            new crack 

      Description: This function controls the process of adding a new 
          finite volume internal flaw. 


    NewFiniteVolSurfFlaw - add a new finite volume surface flaw 

      void NewFiniteVolSurfFlaw(
              const int     num_tips,
              const int     *tip_ids,
              const int     *tip_indices,
              const int     num_points,
              const CrackPt *constpts)

        num_tips    - (in)  number of crack-tips in the flaw 
        tip_ids     - (in)  list of id's to assign to crack tips 
        tip_indices - (in)  location in the crack-point array of the 
                            nodes that will be crack tips 
        num_points  - (in)  number of points in the crack-point list 
        pts         - (in)  list of crack-points, which define the 
                            new crack 

      Description: This function controls the process of adding a new 
          finite volume surface crack. 


    ZeroVolFlawGrowth - extend an existing crack 

      void ZeroVolFlawGrowth(
              const int     crack_tip_id,
              const int     num_points,
              const CrackPt *constpts)

        crack_tip_id - (in)  id of the crack-tip 
        num_points   - (in)  number of points in the crack-point list 
        pts          - (in)  list of crack-points, which define the 
                             crack extension 

      Description: This function controls the process of extending an 
          existing crack. 


    NumNodes - return the number of nodes in the mesh 

      int NumNodes() const

      Description: This function returns the current number of nodes 
          in a mesh (typically called after crack insertion or 
          growth). 

      Return Value: number of nodes 


    NumElements - return the number of elements in the mesh 

      int NumElements() const

      Description: This function returns the current number of 
          elements in a mesh (typically called after crack insertion 
          or growth). 

      Return Value: number of elements 


    NumCrackTips - return the number of crack-tips in a mesh 

      int NumCrackTips() const

      Description: This function returns the current number 
          crack-tips defined in a mesh. 

      Return Value: number of crack tips. 


    NumTemplateElems - returns the number template elements 

      int NumTemplateElems() const

      Description: This function returns the current total number of 
          elements in all crack-tip templates. 

      Return Value: number of crack-tip template elements 


    GetNodes - return a list of all mesh nodes 

      ArbMshNode *GetNodes() const

      Description: This function returns a list of all the nodes 
          current defined for the mesh (typically called after crack 
          insertion or growth). 

      Return Value: A list of all mesh nodes. Ownership of this 
          memory passes to the client, which must eventually call 
          delete []. 


    GetElements - return a list of all mesh elements 

      ArbMshElement2D *GetElements() const

      Description: This function returns a list of all the elements 
          current defined for the mesh (typically called after crack 
          insertion or growth). 

      Return Value: A list of all mesh elements. Ownership of this 
          memory passes to the client, which must eventually call 
          delete []. 


    GetCrackTipIds - return a list of all crack-tip id's 

      int *GetCrackTipIds() const

      Description: This function returns a list of all the currently 
          defined crack-tip id's. 

      Return Value: A list of all crack-tip id's. Ownership of this 
          memory passes to the client, which must eventually call 
          delete []. 


    GetCrackTipNodes - return an array of the crack tip node id's 

      int *GetCrackTipNodes() const

      Description: This function returns a list of crack-tip node 
          id's These will be in the same order as the array' returned 
          from GetCrackTipIds. Only one crack-tip node id is returned 
          for each crack tip, even if the crack tip nodes are not 
          constrained. Ownership of this memory passes to the client, 
          which should eventually call delete []. 

      Return Value: A list of crack-tip nodes. Ownership of this 
          memory passes to the client, which should eventually call 
          delete []. 


    GetCrackTipNodeList - returns a list of the nodes at a crack tip 

      int *GetCrackTipNodeList(
              const int crack_id,
              int       *num_nodes)

        crack_id  - (in)  crack-tip id 
        num_nodes - (out) number of node id's in the list 

      Description: This function returns a list of the crack-tip 
          nodes for the given crack tip. If this is a constrained 
          crack-tip node there will only be one element in the list. 
          There will be more than one node for unconstrained crack 
          tip nodes. 

      Return Value: A list of crack-tip nodes. Ownership of this 
          memory passes to the client, which should eventually call 
          delete []. 

      Exceptions:
          CArbMshCInvalidCrack - invalid crack id


    GetTemplateElemList - return a list of crack-tip template element 
                          id's 

      int *GetTemplateElemList() const

      Description: This function returns a list of the element id's 
          for all elements in crack-tip templates. 

      Return Value: A list of crack-tip template elements. Ownership 
          of this memory passes to the client, which should 
          eventually call delete []. 


    GetTemplateElemIds - return a list of template element crack id's 

      int *GetTemplateElemIds() const

      Description: This function returns an array of crack id's 
          associated with template elements. There is a one-to-one 
          corespondance with entries in the list returned by 
          GetTemplateElemList. 

      Return Value: A list of the crack id's associated with 
          crack-tip template elements. Ownership of this memory 
          passes to the client, which should eventually call delete 
          []. 


    NumDeletedElements - return the number of elements deleted during 
                         remeshing 

      int NumDeletedElements() const

      Description: This function returns the number of elements 
          deleted during the most recent remeshing. 

      Return Value: number of deleted elements 


    GetDeletedElements - return a list of deleted elements 

      int *GetDeletedElements() const

      Description: This function returns a list of the id's of 
          elements deleted during the most recent remeshing. 

      Return Value: A list of the id's of deleted elements. Ownership 
          of this memory passes to the client, which should 
          eventually call delete []. 


    NumCreatedElements - return the number of elements created during 
                         remeshing 

      int NumCreatedElements() const

      Description: This function returns the number of elements 
          created during the most recent remeshing. 

      Return Value: number of created elements 


    GetCreatedElements - return a list of created elements 

      int *GetCreatedElements() const

      Description: This function returns a list of the id's of 
          elements created during the most recent remeshing. 

      Return Value: A list of the id's of created elements. Ownership 
          of this memory passes to the client, which should 
          eventually call delete []. 


    GetUpdatedBoundaryData - return information about updated 
                             boundaries 

      void GetUpdatedBoundaryData(
              int          *num_old,
              BoundaryData **old_data,
              int          *num_new,
              BoundaryData **new_data)

        num_old  - (out) number of items in the old boundary list 
        old_data - (out) list of boundary segments that existed in 
                         the old mesh and have been deleted or 
                         modified in the new mesh. 
        num_new  - (out) number of items in the new boundary list 
        new_data - (out) list of the boundary segments that exist in 
                         the new mesh, that do not correspond to 
                         segments in the old mesh. 

      Description: This function returns information about boundaries 
          that were updated during remeshing. Ownership of this 
          old_data and new_data lists passes to the client, which 
          should eventually call delete [] for both. 


    SetDebugDisplayFlags - set debug display flags 

      void SetDebugDisplayFlags(DDFlags flags)

        flags - (in)  debug display flag 

      Description: This function sets the flag word that controls 
          what information is displayed for graphical debugging. 


  Public Member Variables:

        int DebugDisplayFlags - flag word for graphical debugging 
                output 


PRIVATE INTERFACE

  Private Data Structures:

    struct DeleteSegData

      This data structure is used to to store information used when 
      deleting elements for remeshing. Each instance stores information 
      for one segment of the crack definition. 

      Member Variables:

        CArbCoord2D end_coords[2] - coordinates of the end 
            points of the segment 

        double end_radii[2] - radii of the delete circle 
            centered at the end points of the segment 

        double end_radii_sqr[2] - square of the end radii size 

        double line_coefs[4][3] - Coefficient for the line 
            equations (ax + by + c = 0) for the segment's 
            delete trapezoid. The four vertices of the 
            trapezoid are defined by the end points of the two 
            segments that are mutualy tangent to the two 
            end-point delete circles 

        bool degenerate_case - true if one of the delete 
            circles is completely contained in the other. 


  Private Member Functions:

    FindTipNodeList - return a list of crack-tip node id's 

      int *FindTipNodeList(
              const int     tip_node_id,
              int           *num_nodes,
              CArbMshTopo2D *topo) const

        tip_node_id - (in)  crack-tip id 
        num_nodes   - (out) number of crack-tip nodes 
        topo        - (in)  mesh topology object 

      Description: This function returns an array of the crack tip 
          nodes for the given crack tip. If this is a constrained 
          crack-tip node there will only be one element in the list. 
          There will be more than one node for unconstrained crack 
          tip nodes. 

      Return Value: A list of crack-tip node id's. Ownership of the 
          list passes to the caller. 


    GetTipData - return crack-tip data 

      CrackTipData *GetTipData(int crack_tip_id)

        crack_tip_id - (in)  crack-tip id 

      Description: This function returns the current crack-tip 
          parameters for the given crack-tip id. 

      Return Value: The crack-tip parameters for this crack. 
          Ownership passes to the caller. 

      Exceptions:
          CArbMshCInvalidCrack - invalid crack id


    NewNode - generate a new node 

      ArbMshNode NewNode(
              double x,
              double y)

        x - (in)  node's x coordinate 
        y - (in)  node's y coordinate 

      Description: This function generates a new node and assigns it 
          the next available node id. 

      Return Value: the new node object 


    NewElemId - return a new element id 

      int NewElemId()

      Description: This function returns the next available element 
          id and increments the current max id counter. 

      Return Value: a unique element id 


    NewNodeId - return a new element id 

      int NewNodeId()

      Description: This function returns the next available node id 
          and increments the current max id counter. 

      Return Value: a unique node id 


  Private Member Variables:

    CArbHashTable<int,CArbCoord2D>* NodeTable - hash table the maps 
            node id's to node coordinates 

    CArbHashTable<int,ArbMshElement2D>* ElemTable - hash table that 
            maps element id's to element descriptions 

    CArbHashTable<int,int>* CrackTable - hash table that maps 
            crack-tip id's to crack-tip node numbers. 

    CArbMshTopo2D *MshTopo - current mesh topology object 

    ElemShape ElemShapeType - current element shape type 
            (S_TRIANGLE or S_QUADRILATERAL) 

    int MaxNodeId - maximum assigned node id 

    int MaxElemId - maximum assigned element id 

    bool FirstElem - if true there are no elements in the mesh yet 

    ArbMshOrder Order - polynomial order of mesh elements (LINEAR 
            or QUADRATIC) 

    bool CanResizeTemplate - if true the crack-tip template can be 
            resized during remeshing 

    CArbSet<int>* NoRmshMats - list of materials that will not be 
            modified during remeshing 

    CArbHashTable<int,int>* TemplateElems - hash table that maps 
            element id's to a crack-tip node id 

    CrackTipData DefaultTipData - crack-tip parameters to use if 
            not specifically overidden 

    CArbHashTable<int,CrackTipData>* TipData - hash table that maps 
            a crack-tip id to the associated crack-tip parameters 

    CArbArray<int>* DeletedElems - array of elements deleted during 
            remeshing 

    CArbArray<int>* CreatedElems - array of elements created during 
            remeshing 

    CArbHashTable<ArbEdgeKey,BoundaryData>* OldBdryTable - hash 
            table that maps boundary segment node id's to 
            associated boundary data for the old mesh 

    CArbHashTable<ArbEdgeKey,BoundaryData>* NewBdryTable - hash 
            table that maps boundary segment node id's to 
            associated boundary data for the new mesh 

    double TriDelFactor - delete factor used during remeshing with 
            triangles 

    double QuadDelFactor - delete factor used during remeshing with 
            quadrilaterals 


NON-MEMBER OPERATORS

operator == - equivalence operator for CrackTipData structures 

  inline int operator == (
          const CArbMshCrackRegion2D::CrackTipData &t1,
          const CArbMshCrackRegion2D::CrackTipData &t2)

    t1 - (in)  first operator 
    t2 - (in)  second operator 

  Description: This operator performs an equivalence check for 
      CrackTipData structures. 

  Return Value: nonzero indicates equivalence 


operator == - equivalence operator for BoundaryData structures 

  inline int operator == (
          const CArbMshCrackRegion2D::BoundaryData &b1,
          const CArbMshCrackRegion2D::BoundaryData &b2)

    b1 - (in)  first operator 
    b2 - (in)  second operator 

  Description: This operator performs an equivalence check for 
      BoundaryData structures. 

  Return Value: nonzero indicates equivalence 

*/

#endif

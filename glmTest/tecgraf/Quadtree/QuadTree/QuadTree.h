#ifndef _DATA_QUAD_TREE_H_
#define _DATA_QUAD_TREE_H_

#include "Data/Definitions.h"
#include "Data/Point2D.h"
#include "Data/BoxTree.h"
#include "Data/QuadTreeCell.h"

namespace Data
{
	class QuadTree : public Data::BoxTree
	{
	public:

		QuadTree(const Point2D &min, const Point2D &max, ULInt id = 0);
		QuadTree(Box *box, QuadTreeCell *root, ULInt id = 0);
		virtual ~QuadTree();

		virtual GraphNode *node(ULInt id) const;

		using BoxTree::in;
		virtual bool in(Real x, Real y) const;

		using BoxTree::optIn;
		virtual bool optIn(Real x, Real y) const;

		using BoxTree::on;
		virtual bool on(Real x, Real y) const;

		using BoxTree::optOn;
		virtual bool optOn(Real x, Real y) const;

		using BoxTree::out;
		virtual bool out(Real x, Real y) const;

		using BoxTree::optOut;
		virtual bool optOut(Real x, Real y) const;
	};
}

#endif //#ifndef _DATA_QUAD_TREE_H_

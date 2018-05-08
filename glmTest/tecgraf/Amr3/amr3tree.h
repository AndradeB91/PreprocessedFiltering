/*
**
** File: amr3tree.h 
**
** Author: Mauricio Riguette Mediano
**
** Date: 28/03/1997
**
*/

#ifndef amr3treeDEFINED
#define amr3treeDEFINED

#include <math.h>
#define RT_MAX_NUM_NODES 100

#ifndef ABS
#define ABS(x)  (((x)>= 0)?(x): -(x))
#endif

#ifndef MIN
#define MIN(x,y)  (((x)<(y))?(x):(y))
#endif

#ifndef MAX
#define MAX(x,y)  (((x)>(y))?(x):(y))
#endif

#ifndef double_MAX
#define double_MAX 1.6e+300 			 
#endif

#ifndef double_MIN
#define double_MIN 1.6e-300 			 
#endif

class amBox
{
  public:
    double Xmin,Ymin,Zmin,Xmax,Ymax,Zmax;
    void Dup(amBox *B);
    int Equal(amBox *);
    int Disjoint(amBox *);
    double Volume();
    void ExtendVolume(amBox *);
};

class amBoxId:public amBox
{
public:
  void *P;
};

class amNodeBoxId
{
  int  Num;
  int  Max;
  int  Min;
protected:
  amBoxId *Vet;
public:
  amNodeBoxId(int Pmin, int Pmax);
  virtual ~amNodeBoxId();
  int  Number () {return Num;};
  void Start () {Num=0;};
  void Dup(amNodeBoxId *Fromn);
  void BBox (amBox *B);
  void MakeBoxId(amBoxId *Node);
  void GetPos(amBoxId *E,int i);
  void PutPos(amBoxId *E,int i);
  void Insert(amBoxId *E,int Posit);
  void Append(amBoxId *E);
  void Remove(int i);
  void Remove(amBoxId *E,int Cont);
  void Remove(amBoxId *E);
  void ChooseBoxId(amBoxId *ebox,amBoxId *Out,int *Posit);
  void UpdateBox(amBox *B,int Posit);
  int  LocateBoxId(amBoxId *In,int *Posit);
  void SplitNodeR3tree(amNodeBoxId *Newnode);
};

class amR3tree 
{
protected:
  int  Max;
  int  Min;
  int Level;
  amNodeBoxId *Root;
  
  /* Estes atributos mantem uma pilha (log(n)) de nos da arvore, para    */
  /* retornar passo a passo (Result) cada box selecionado  por uma busca */
  /* na R3tree (Search).                                                 */

  amNodeBoxId *Vnode[RT_MAX_NUM_NODES];
  int Vpos[RT_MAX_NUM_NODES];
  int NumL;
  amBox SearchVal;

  int Remove(amNodeBoxId *Node,amBoxId *In,int Lev);
  void Delete(amNodeBoxId *Node,int Lev);
  void Number(amNodeBoxId *Node,int *num,int Lev);
  void MergeRoot();
  void Merge(amNodeBoxId *Father,amNodeBoxId *Node,short Posit);
  int Insert(amNodeBoxId *Node,amBoxId *In,amBoxId *Out,int Lev);
  void SplitNodeR3tree(amNodeBoxId *Node, amBoxId *Out);
  void Print( amNodeBoxId *Node, int ident, int *num, int Lev);

public:
  amR3tree(int Nmin,int Nmax);
  virtual ~amR3tree();
  void Insert(amBoxId *newentry);
  int Remove(amBoxId *newentry);

  void Print();
  void Create();
  void Delete();
  int  Number();
  void BBox(amBox *B);
  void Search(amBox *B);
  int  Result(amBoxId *B);
  void SearchAll();
  int  ResultAll(amBoxId *B);
};

#endif

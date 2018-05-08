/*
**
** File: martree.cpp
**
** Author: Mauricio Riguette Mediano
**
** Date: 28/03/1997
**
*/

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "amr3tree.h"

// -------------------------------------------
void amBox::Dup(amBox *B)
{
  Xmin = B->Xmin;
  Ymin = B->Ymin;
  Zmin = B->Zmin;
  Xmax = B->Xmax;
  Ymax = B->Ymax;
  Zmax = B->Zmax;
}
// -------------------------------------------
int amBox::Equal(amBox *B) 
{
  return( Xmin==B->Xmin&&Ymin==Ymin&&Zmin==Zmin&& 
          Xmax==B->Xmax&&Ymax==Ymax&&Zmax==Zmax);
}
// -------------------------------------------
int amBox::Disjoint(amBox *B) 
{
  if( Xmax<B->Xmin||Xmin>B->Xmax|| 
      Ymax<B->Ymin||Ymin>B->Ymax|| 
      Zmax<B->Zmin||Zmin>B->Zmax)
    return(1);
  return(0);
}
// -------------------------------------------
double amBox::Volume() 
{
  return((Xmax-Xmin)*(Ymax-Ymin)*(Zmax-Zmin));
}
// -------------------------------------------
void amBox::ExtendVolume(amBox *B)
{
  Xmin = MIN(Xmin,B->Xmin);
  Ymin = MIN(Ymin,B->Ymin);
  Zmin = MIN(Zmin,B->Zmin);
  Xmax = MAX(Xmax,B->Xmax);
  Ymax = MAX(Ymax,B->Ymax);
  Zmax = MAX(Zmax,B->Zmax);
}
// --------------------------------------------------------
amNodeBoxId::amNodeBoxId(int Pmin, int Pmax) 
{ 
  if(Pmin>Pmax) exit(1);
  Min = Pmin;
  Max = Pmax;
  Num=0;
  Vet = new amBoxId[Max+1];
  if(Vet==NULL) exit(1);
}
// --------------------------------------------------------
amNodeBoxId::~amNodeBoxId()
{
  delete[] Vet;
}
// --------------------------------------------------------
void amNodeBoxId::Dup(amNodeBoxId* Fromn)
{
  int i;
  Num=Fromn->Num;
  Max=Fromn->Max;
  Min=Fromn->Min;
  if (Vet!=NULL) delete[] Vet;
  if (Fromn->Vet==NULL) {
    Vet=NULL;
    return;
  };
  Vet = new amBoxId[Max+1];
  if(Vet==NULL) exit(1);
  for(i=0;i<Num;i++)
    Vet[i]=Fromn->Vet[i];
};
// --------------------------------------------------------
void amNodeBoxId::BBox(amBox *B)
{
  if (Num==0) exit(1);
  *B = Vet[0];
  for(int cont=1;cont<Num;cont++)
    B->ExtendVolume(&Vet[cont]);
};
// --------------------------------------------------------
void amNodeBoxId::MakeBoxId(amBoxId *E)
{
  BBox(E);
  E->P=this;
}
// --------------------------------------------------------
void amNodeBoxId::GetPos(amBoxId *E,int i)
{
  if (Vet==NULL) exit(1);
  if (i > Num) exit(1);
  *E = Vet[i];
};
// --------------------------------------------------------
void amNodeBoxId::PutPos(amBoxId *E,int i)
{
  if(Vet==NULL) exit(1);
  if(i > Num) exit(1);
  Vet[i]=*E;
};
// --------------------------------------------------------
void amNodeBoxId::Insert(amBoxId *E,int Posit)
{
  if(Vet==NULL) exit(1);
  if(Posit>Num) exit(1);
  for(int cont=Num;cont>Posit;cont--)
    Vet[cont]=Vet[cont-1];
  Num++;
  Vet[Posit]=*E;
};
// --------------------------------------------------------
void amNodeBoxId::Append(amBoxId *E)
{
  if(Num > Max) exit(1);
  PutPos(E,Num);
  Num++;
};
// --------------------------------------------------------
void amNodeBoxId::Remove(int i)
{
  if (Vet==NULL) exit(1);
  if (i > Num) exit(1);
  Num--;
  for(int cont=i;cont<Num;cont++)
    Vet[cont]=Vet[cont+1];
};
// --------------------------------------------------------
void amNodeBoxId::Remove(amBoxId *E,int Cont)
{
  GetPos(E,Cont);
  Remove(Cont);
};
// --------------------------------------------------------
void amNodeBoxId::Remove(amBoxId *E)
{
  GetPos(E,Num-1);
  Num--;
}
// --------------------------------------------------------
void amNodeBoxId::ChooseBoxId(amBoxId *Loc,amBoxId *Out,int *Posit)
{
  amBoxId eaux;
  double smallerarea,areabefore;

  if (Num==0) exit(1);
  *Posit=0;
  GetPos(&eaux,*Posit);
  *Out = eaux;
  areabefore = eaux.Volume()+Loc->Volume();
  eaux.ExtendVolume(Loc);
  smallerarea = eaux.Volume() - areabefore;
  for(int cont=1;cont<Number();cont++) {
    GetPos(&eaux,cont);
    areabefore = eaux.Volume()+Loc->Volume();
    eaux.ExtendVolume(Loc);
    if(eaux.Volume() - areabefore < smallerarea) {
      GetPos(Out,cont);
      *Posit = cont;
      smallerarea = eaux.Volume() - areabefore;
    }
  }
}
// --------------------------------------------------------
void amNodeBoxId::UpdateBox(amBox *B,int Posit)
{
  amBoxId E;
  GetPos(&E,Posit);
  E.Dup(B);
  PutPos(&E,Posit);
}
// --------------------------------------------------------
int amNodeBoxId::LocateBoxId(amBoxId *In,int *Posit)
{
  for(int cont=0;cont<Number();cont++)
    if(Vet[cont].P==In->P) {
      *Posit=cont;
      return(0);
    }
  return(1);
}
// --------------------------------------------------------
void amNodeBoxId::SplitNodeR3tree(amNodeBoxId *Newnode)
{
  amNodeBoxId auxnode(Min,Max);
  double extraarea;
  double areabefore;
  int imax, jmax;
  int i, j;
  amBoxId ei, ej;
  amBoxId eiaux,ejaux;

  imax=0;       //Inicializa entradas
  jmax=1;
  GetPos(&eiaux, 0); 
  GetPos(&ejaux, 1);
  areabefore = eiaux.Volume() + ejaux.Volume();
  eiaux.ExtendVolume(&ejaux);
  extraarea = eiaux.Volume() - areabefore;

  for(i=0;i<Number();i++) {  // Escolhe entradas iniciais.
    for(j=i+1;j<Number();j++) {
      GetPos(&eiaux, i);
      GetPos(&ejaux, j);
      areabefore = eiaux.Volume() + ejaux.Volume();
      eiaux.ExtendVolume(&ejaux);
      if(eiaux.Volume() - areabefore > extraarea) {
        imax = i;
        jmax = j;
        extraarea = eiaux.Volume() - areabefore;
      };
    } /* endfor */
  } /* endfor */

  if(imax>jmax) {
    Remove(&ei, imax);
    Remove(&ej, jmax);
  } else {
    Remove(&ej, jmax);
    Remove(&ei, imax);
  }

  while(Number()>0) {  // Move entradas para no' auxiliar.
    amBoxId E;
    Remove(&E);
    auxnode.Append(&E);
  } /* endfor */

  Append(&ei);
  Newnode->Append(&ej);
  for(i=0;i<Min-1;i++) {  // Insere minimo por no'.
    amBoxId E;
    int choose;
    auxnode.ChooseBoxId(&ei,&E,&choose);
    auxnode.Remove(&E, choose);
    ei.ExtendVolume(&E);
    Append(&E);
    auxnode.ChooseBoxId(&ej,&E,&choose);
    auxnode.Remove(&E, choose);
    ej.ExtendVolume(&E);
    Newnode->Append(&E);
  } /* endfor */

  GetPos(&ei,0);
  Newnode->GetPos(&ej,0);
  while(auxnode.Number()>0) {
    amBoxId E;
    amBox biaux,bjaux;
    auxnode.Remove(&E);
    biaux = ei;
    biaux.ExtendVolume(&E);
    bjaux = ej;
    bjaux.ExtendVolume(&E);
    if(biaux.Volume()-ei.Volume()<
       bjaux.Volume()-ej.Volume()){
      Append(&E);
      ei.ExtendVolume(&E);
    } else {
      Newnode->Append(&E);
      ej.ExtendVolume(&E);
    }
  } /* endfor */
}
// --------------------------------------------------------
amR3tree::amR3tree(int Nmin,int Nmax)
{
  Min=Nmin;
  Max=Nmax;
}
// --------------------------------------------------------
amR3tree::~amR3tree()
{
}
// --------------------------------------------------------
void amR3tree::MergeRoot()
{
  amBoxId enode;
  Root->GetPos(&enode,0);
  delete Root;
  Root=(amNodeBoxId *)enode.P;
}
// --------------------------------------------------------
void amR3tree::BBox(amBox *B)
{
  Root->BBox(B);
}
// --------------------------------------------------------
int amR3tree::Remove(amBoxId *newentry)
{
  if (Remove(Root,newentry,Level)==2) return(2);
  if (Level>0 && Root->Number()==1){
    MergeRoot();
    Level--;
  }
  return(0);
};
// --------------------------------------------------------
int amR3tree::Remove(amNodeBoxId *Node,amBoxId *In,int Lev)
{
  int posit;
  if (Lev==0) {
    if(Node->LocateBoxId(In,&posit)==0) {
      Node->Remove(posit);
      return(0);
    };
  } else {
    amBoxId eson;
    amNodeBoxId *son;
    for(posit=0;posit<Node->Number();posit++) {
      Node->GetPos(&eson,posit);
      if (In->Disjoint(&eson))
        continue;
      son=(amNodeBoxId *)eson.P;
      if(Remove(son,In,Lev-1)==2) {
        continue;
      } else if(son->Number()<Min) {
        Merge(Node,son,posit);
      } else {
        son->MakeBoxId(&eson);
        Node->PutPos(&eson,posit);
      };
      return(0);
    }
  };
  return(2);
}
// --------------------------------------------------------
void amR3tree::Insert(amBoxId *newentry)
{
  amBoxId out;
  if(Insert(Root,newentry,&out,Level)==2) {
    amBoxId oldentry;
    amNodeBoxId *newroot= new amNodeBoxId(Min,Max);
    if (newroot==NULL) exit(1);
    Root->MakeBoxId(&oldentry);
    newroot->Append(&oldentry);
    newroot->Append(&out);
    Root=newroot;
    Level++;
  } /* endif */
}
// --------------------------------------------------------
void amR3tree::SplitNodeR3tree(amNodeBoxId *Node, amBoxId *Out)
{
  amNodeBoxId *newnode= new amNodeBoxId(Min,Max);
  if (newnode==NULL) exit(1);
  Node->SplitNodeR3tree(newnode);
  newnode->MakeBoxId(Out);
};
// --------------------------------------------------------
int amR3tree::Insert(amNodeBoxId *Node,amBoxId *In,amBoxId *Out,int Lev)
{
  if (Lev==0) {
    Node->Append(In);
    if (Node->Number() > Max) {
      SplitNodeR3tree(Node,Out);
      return(2);
    }
  } else {
    int status,posit;
    amBoxId eson;
    amNodeBoxId *son;
    Node->ChooseBoxId(In,&eson,&posit);
    son=(amNodeBoxId *)eson.P;
    status = Insert(son,In,Out,Lev-1);
    son->MakeBoxId(&eson);
    Node->PutPos(&eson,posit);
    if (status==2) {
      Node->Append(Out);
      if(Node->Number() > Max) {
        SplitNodeR3tree(Node,Out);
        return(2);
      }
    }
  }
  return(0);
}
// --------------------------------------------------------
void amR3tree::Merge(amNodeBoxId *Father,amNodeBoxId *Node,short Posit)
{

//
// Alterar para procurar o no' mais proximo para 
// fazer o merge.
//

  int i;
  amBoxId auxentry;
  amNodeBoxId *auxnode = NULL;
  amBoxId eaux;
  short auxposit = 0;

  if(Posit + 1 < Father->Number()) auxposit = Posit + 1;
  else if(Posit > 0) auxposit = Posit - 1;

  Father->GetPos(&auxentry,auxposit);
  auxnode=(amNodeBoxId *)auxentry.P;
  if(auxnode->Number() + Node->Number() <= Max) {
    for(i=0;i<Node->Number();i++) {
      Node->GetPos(&eaux,i);
      auxnode->Append(&eaux);
    }
    auxnode->MakeBoxId(&auxentry);
    Father->PutPos(&auxentry,auxposit);
    delete(Node);
    Father->Remove(Posit);
  } else {
    while(Node->Number() < Min) {
      auxnode->Remove(&eaux);
      Node->Append(&eaux);
    }
    auxnode->MakeBoxId(&auxentry);
    Father->PutPos(&auxentry,auxposit);
    Node->MakeBoxId(&auxentry);
    Father->PutPos(&auxentry,Posit);
  }
}
// --------------------------------------------------------
void amR3tree::Create()
{
  Root=new amNodeBoxId(Min,Max);
  if (Root==NULL) exit(1);
  Level=0;
}
// --------------------------------------------------------
void amR3tree::Delete()
{
  if (Level>0)
    Delete(Root,Level);
  delete (Root);
};
// --------------------------------------------------------
void amR3tree::Delete(amNodeBoxId *Node,int Lev)
{
  if (Lev>0) {
    amBoxId e;
    amNodeBoxId *naux;
    while(Node->Number()>0) {
      Node->Remove(&e);
      naux=(amNodeBoxId *)e.P;
      Delete(naux,Lev-1);
      delete(naux);
    }
  }
}
// --------------------------------------------------------
int amR3tree::Number()
{
  int num=0;
  Number(Root,&num,Level);
  return(num);
};
// --------------------------------------------------------
void amR3tree::Number(amNodeBoxId *Node,int *Num,int Lev)
{
  amBoxId e;
  if (Lev==0) { 
    *Num += Node->Number();
  } else {
    amNodeBoxId *naux;
    while(Node->Number()>0) {
      Node->Remove(&e);
      naux=(amNodeBoxId *)e.P;
      Number(naux,Num,Lev-1);
    }
  }
}
// --------------------------------------------------------
void amR3tree::Print()
{
  int num=0;
  Print(Root,0,&num,Level);
  printf("Numero de elementos %d\n", num);
}
// --------------------------------------------------------
void amR3tree::Print( amNodeBoxId *Node, int ident, int *Num, int Lev)
{
  amBoxId e;
  int i,j;
  amBox* vet;
  vet = new amBox[Max];
  if(vet==NULL) exit(1);
  for( i=0;i<ident;i++) {    
    printf( " ");
  }
  printf( "No #%d ->\n", Node->Number());
  for( i=0;i<Node->Number();i++) { 
    Node->GetPos(& e, i);
    e.Dup(&vet[i]);
    for( j=0;j<ident+2;j++) {    
      printf( " ");
    }
    printf("#%d %g %g %g %g %g %g\n",i, vet[i].Xmin, vet[i].Ymin, vet[i].Zmin,
                                        vet[i].Xmax, vet[i].Ymax, vet[i].Zmax);
  }
  if (Lev==0) {
    *Num+=Node->Number();
  } else {
    amNodeBoxId *naux;
    for( i=0;i<Node->Number();i++) { 
      amBox B;
      Node->GetPos(& e, i);
      naux=(amNodeBoxId *)e.P;
      naux->BBox(&B);
      if(!B.Equal(&vet[i])) printf( "\nErro Box\n");
      Print(naux, ident+2, Num, Lev -1);
    }
  }
  delete[] vet;
}
// --------------------------------------------------------
void amR3tree::Search(amBox *B)
{
  Vnode[0]=Root;
  Vpos[0]=0;
  NumL=0;
  SearchVal=*B;
};
// --------------------------------------------------------
int amR3tree::Result(amBoxId *B)
{
  amBoxId currkey;
  currkey.P = NULL;

  for(;;) {
    for (;Vpos[NumL]<Vnode[NumL]->Number();Vpos[NumL]++) {
      Vnode[NumL]->GetPos(&currkey,Vpos[NumL]);
      if (!currkey.Disjoint(&SearchVal))
          break;
    }
    if (Vpos[NumL]==Vnode[NumL]->Number()) {
      NumL--;
      if (NumL<0) break;
      else Vpos[NumL]++;
    } else {
      Vnode[NumL]->GetPos(&currkey,Vpos[NumL]);
      if (NumL==Level) break;
      NumL++;
      Vnode[NumL]=(amNodeBoxId *)currkey.P;
      Vpos[NumL]=0;
    }
  }
  if (NumL==-1) return(1);
  *B=currkey;
  Vpos[NumL]++;
  return(0);
};
// --------------------------------------------------------
void amR3tree::SearchAll()
{
  Vnode[0]=Root;
  Vpos[0]=0;
  NumL=0;
};
// --------------------------------------------------------
int amR3tree::ResultAll(amBoxId *B)
{
  amBoxId currkey;
  for(;;) {
    Vnode[NumL]->GetPos(&currkey,Vpos[NumL]);
    if (Vpos[NumL]==Vnode[NumL]->Number()) {
      NumL--;
      if (NumL<0) break;
      else Vpos[NumL]++;
    } else {
      Vnode[NumL]->GetPos(&currkey,Vpos[NumL]);
      if (NumL==Level) break;
      NumL++;
      Vnode[NumL]=(amNodeBoxId *)currkey.P;
      Vpos[NumL]=0;
    };
  };
  if (NumL==-1) return(1);
  *B=currkey;
  Vpos[NumL]++;
  return(0);
};
// --------------------------------------------------------

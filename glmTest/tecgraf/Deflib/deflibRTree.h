// deflib.h : Arquivo de interface.
//
//     DESENVOLVIDO POR : Ricardo Morrot
//               e-MAIL : morrot@tecgraf.puc-rio.br (morrot@gmail.com)
//
//              PROJETO : Tecgraf/MVGEO/MTool3D
//                  ANO : 2011
//
//           BIBLIOTECA : RTREE
//            DESCRIÇÃO : A biblioteca RTREE passa a ser gerada da forma de duas formas
//                        distintas, dinâmica (DLL - Export/Import) e estática (lib).
//
//                        Utilize os símbolos de pré-processamento para:
//                        - RTREE_DLL_IMPORTS, para importar(utilizar) num aplicativo.
//                        - RTREE_DLL_EXPORTS, para gerar a DLL.
//                        - omitindo RTREE_DLL_IMPORTS e RTREE_DLL_EXPORTS é gerada a 
//                          biblioteca estática (lib).
//
//      DATA DE CRIAÇÃO : 22/07/2011
//
// ÚLTIMAS ATUALIZAÇÕES
//
//         [22/07/2011] - Morrot, arquivo gerado.
//
//////////////////////////////////////////////////////////////////////////////////////////
#ifndef __RTREE_DEFLIB_H__
#define __RTREE_DEFLIB_H__

#ifndef RTREE_DLL_EXPORTS
#  ifndef RTREE_DLL_IMPORTS
#    define RTREE_API
#  else
#    define RTREE_API __declspec(dllimport)
#  endif
#else
#  define RTREE_API __declspec(dllexport)
#endif

#endif //__RTREE_DEFLIB_H__

// deflib.h : Arquivo de interface.
//
//     DESENVOLVIDO POR : Ricardo Morrot
//               e-MAIL : morrot@tecgraf.puc-rio.br (morrot@gmail.com)
//
//              PROJETO : Tecgraf/MVGEO/MTool3D
//                  ANO : 2011
//
//           BIBLIOTECA : RTREE
//            DESCRI��O : A biblioteca RTREE passa a ser gerada da forma de duas formas
//                        distintas, din�mica (DLL - Export/Import) e est�tica (lib).
//
//                        Utilize os s�mbolos de pr�-processamento para:
//                        - RTREE_DLL_IMPORTS, para importar(utilizar) num aplicativo.
//                        - RTREE_DLL_EXPORTS, para gerar a DLL.
//                        - omitindo RTREE_DLL_IMPORTS e RTREE_DLL_EXPORTS � gerada a 
//                          biblioteca est�tica (lib).
//
//      DATA DE CRIA��O : 22/07/2011
//
// �LTIMAS ATUALIZA��ES
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

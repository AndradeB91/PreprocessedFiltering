// deflib.h : Arquivo de interface.
//
//     DESENVOLVIDO POR : Ricardo Morrot
//               e-MAIL : morrot@tecgraf.puc-rio.br (morrot@gmail.com)
//
//              PROJETO : Tecgraf/MVGEO/MTool3D
//                  ANO : 2011
//
//           BIBLIOTECA : MSHAUX
//            DESCRICAO : A biblioteca MSHAUX passa a ser gerada da forma de duas formas
//                        distintas, dinamica (DLL - Export/Import) e estatica (lib).
//
//                        Utilize os simbolos de pre-processamento para:
//                        - MSHAUX_DLL_IMPORTS, para importar(utilizar) num aplicativo.
//                        - MSHAUX_DLL_EXPORTS, para gerar a DLL.
//                        - omitindo MSHAUX_DLL_IMPORTS e MSHAUX_DLL_EXPORTS e gerada a 
//                          biblioteca estatica (lib).
//
//      DATA DE CRIACAO : 23/07/2011
//
// ULTIMAS ATUALIZACOES
//
//         [23/07/2011] - Morrot, arquivo gerado.
//
//////////////////////////////////////////////////////////////////////////////////////////
#ifndef __MSHAUX_DEFLIB_H__
#define __MSHAUX_DEFLIB_H__

#ifndef MSHAUX_DLL_EXPORTS
#  ifndef MSHAUX_DLL_IMPORTS
#    define MSHAUX_API
#  else
#    define MSHAUX_API __declspec(dllimport)
#  endif
#else
#  define MSHAUX_API __declspec(dllexport)
#endif

#endif //__MSHAUX_DEFLIB_H__


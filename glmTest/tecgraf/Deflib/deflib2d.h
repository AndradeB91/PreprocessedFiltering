// deflib.h : Arquivo de interface.
//
//     DESENVOLVIDO POR : Ricardo Morrot
//               e-MAIL : morrot@tecgraf.puc-rio.br (morrot@gmail.com)
//
//              PROJETO : Tecgraf/MVGEO/MTool3D
//                  ANO : 2011
//
//           BIBLIOTECA : MSH2D
//            DESCRICAO : A biblioteca MSH2D passa a ser gerada da forma de duas formas
//                        distintas, dinamica (DLL - Export/Import) e estatica (lib).
//
//                        Utilize os simbolos de pre-processamento para:
//                        - MSH2D_DLL_IMPORTS, para importar(utilizar) num aplicativo.
//                        - MSH2D_DLL_EXPORTS, para gerar a DLL.
//                        - omitindo MSH2D_DLL_IMPORTS e MSH2D_DLL_EXPORTS e gerada a
//                          biblioteca estatica (lib).
//
//      DATA DE CRIACAO : 23/07/2011
//
// ULTIMAS ATUALIZACOES
//
//         [23/07/2011] - Morrot, arquivo gerado.
//
//////////////////////////////////////////////////////////////////////////////////////////
#ifndef __MSH2D_DEFLIB_H__
#define __MSH2D_DEFLIB_H__

#ifndef MSH2D_DLL_EXPORTS
#  ifndef MSH2D_DLL_IMPORTS
#    define MSH2D_API
#  else
#    define MSH2D_API __declspec(dllimport)
#  endif
#else
#  define MSH2D_API __declspec(dllexport)
#endif

#endif //__MSH2D_DEFLIB_H__


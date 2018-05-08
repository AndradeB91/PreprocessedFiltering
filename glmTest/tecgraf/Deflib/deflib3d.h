// deflib.h : Arquivo de interface.
//
//     DESENVOLVIDO POR : Ricardo Morrot
//               e-MAIL : morrot@tecgraf.puc-rio.br (morrot@gmail.com)
//
//              PROJETO : Tecgraf/MVGEO/MTool3D
//                  ANO : 2011
//
//           BIBLIOTECA : MSH3D
//            DESCRICAO : A biblioteca MSH3D passa a ser gerada da forma de duas formas
//                        distintas, dinamica (DLL - Export/Import) e estatica (lib).
//
//                        Utilize os simbolos de pre-processamento para:
//                        - MSH3D_DLL_IMPORTS, para importar(utilizar) num aplicativo.
//                        - MSH3D_DLL_EXPORTS, para gerar a DLL.
//                        - omitindo MSH3D_DLL_IMPORTS e MSH3D_DLL_EXPORTS e gerada a 
//                          biblioteca estatica (lib).
//
//      DATA DE CRIACAO : 25/07/2011
//
// ULTIMAS ATUALIZACOES
//
//         [25/07/2011] - Morrot, arquivo gerado.
//
//////////////////////////////////////////////////////////////////////////////////////////
#ifndef __MSH3D_DEFLIB_H__
#define __MSH3D_DEFLIB_H__

#ifndef MSH3D_DLL_EXPORTS
#  ifndef MSH3D_DLL_IMPORTS
#    define MSH3D_API
#  else
#    define MSH3D_API __declspec(dllimport)
#  endif
#else
#  define MSH3D_API __declspec(dllexport)
#endif

#endif //__MSH3D_DEFLIB_H__


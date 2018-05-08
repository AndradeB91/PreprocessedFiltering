// deflib.h : Arquivo de interface.
//
//     DESENVOLVIDO POR : Ricardo Morrot
//               e-MAIL : morrot@tecgraf.puc-rio.br (morrot@gmail.com)
//
//              PROJETO : Tecgraf/MVGEO/MTool3D
//                  ANO : 2011
//
//           BIBLIOTECA : MSHSURF
//            DESCRICAO : A biblioteca MSHSURF passa a ser gerada de duas formas
//                        distintas, dinamica (DLL - Export/Import) e estatica (lib).
//
//                        Utilize os simbolos de pre-processamento para:
//                        - MSHSURF_DLL_IMPORTS, para importar(utilizar) num aplicativo.
//                        - MSHSURF_DLL_EXPORTS, para gerar a DLL.
//                        - omitindo MSHSURF_DLL_IMPORTS e MSHSURF_DLL_EXPORTS e gerada a
//                          biblioteca estatica (lib).
//
//      DATA DE CRIACAOO : 25/07/2011
//
// ULTIMAS ATUALIZACOES
//
//         [25/07/2011] - Morrot, arquivo gerado.
//
//////////////////////////////////////////////////////////////////////////////////////////
#ifndef __MSHSURF_DEFLIB_H__
#define __MSHSURF_DEFLIB_H__

#ifndef MSHSURF_DLL_EXPORTS
#  ifndef MSHSURF_DLL_IMPORTS
#    define MSHSURF_API
#    ifndef EXPIMP_TEMPLATE
#      define EXPIMP_TEMPLATE
#    endif
#  else
#    define MSHSURF_API __declspec(dllimport)
#    ifndef EXPIMP_TEMPLATE
#      define EXPIMP_TEMPLATE extern
#    endif
#  endif
#else
#  define MSHSURF_API __declspec(dllexport)
#  ifndef EXPIMP_TEMPLATE
#    define EXPIMP_TEMPLATE
#  endif
#endif

#endif //__MSHSURF_DEFLIB_H__


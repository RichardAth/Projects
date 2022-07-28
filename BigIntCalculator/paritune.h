#pragma once
/* \src\kernel\gmp\tune.h */
#ifdef LONG_IS_64BIT
#define __AGM_ATAN_LIMIT                 60
#define __DIVRR_GMP_LIMIT                4  /* increased from 4 */
#define __EXPNEWTON_LIMIT                66
#define __EXTGCD_HALFGCD_LIMIT           -1
#define __F2x_MUL_KARATSUBA_LIMIT        11
#define __F2x_MUL_MULII_LIMIT            41
#define __F2xqX_BARRETT_LIMIT            48
#define __F2xqX_DIVREM_BARRETT_LIMIT     97
#define __F2xqX_EXTGCD_LIMIT             97
#define __F2xqX_GCD_LIMIT                605
#define __F2xqX_HALFGCD_LIMIT            127
#define __F2xqX_INVBARRETT_LIMIT         50
#define __F2xqX_REM_BARRETT_LIMIT        101
#define __Flx_BARRETT2_LIMIT             31
#define __Flx_BARRETT_LIMIT              92
#define __Flx_DIVREM2_BARRETT_LIMIT      14
#define __Flx_DIVREM_BARRETT_LIMIT       161
#define __Flx_EXTGCD2_LIMIT              61
#define __Flx_EXTGCD_LIMIT               158
#define __Flx_GCD2_LIMIT                 1409
#define __Flx_GCD_LIMIT                  852
#define __Flx_HALFGCD2_LIMIT             36
#define __Flx_HALFGCD_LIMIT              120
#define __Flx_INVBARRETT2_LIMIT          22
#define __Flx_INVBARRETT_LIMIT           200
#define __Flx_MUL2_KARATSUBA_LIMIT       11
#define __Flx_MUL2_MULII_LIMIT           8
#define __Flx_MUL_KARATSUBA_LIMIT        33
#define __Flx_MUL_MULII_LIMIT            30
#define __Flx_REM2_BARRETT_LIMIT         89
#define __Flx_REM_BARRETT_LIMIT          159
#define __Flx_SQR2_KARATSUBA_LIMIT       15
#define __Flx_SQR2_SQRI_LIMIT            14
#define __Flx_SQR_KARATSUBA_LIMIT        93
#define __Flx_SQR_SQRI_LIMIT             37
#define __FlxqX_BARRETT_LIMIT            17
#define __FlxqX_DIVREM_BARRETT_LIMIT     46
#define __FlxqX_EXTGCD_LIMIT             44
#define __FlxqX_GCD_LIMIT                470
#define __FlxqX_HALFGCD_LIMIT            60
#define __FlxqX_INVBARRETT_LIMIT         22
#define __FlxqX_REM_BARRETT_LIMIT        48
#define __FpXQX_BARRETT_LIMIT            12
#define __FpXQX_DIVREM_BARRETT_LIMIT     30
#define __FpXQX_EXTGCD_LIMIT             28
#define __FpXQX_GCD_LIMIT                191
#define __FpXQX_HALFGCD_LIMIT            35
#define __FpXQX_INVBARRETT_LIMIT         40
#define __FpXQX_REM_BARRETT_LIMIT        30
#define __FpX_BARRETT_LIMIT              38
#define __FpX_DIVREM_BARRETT_LIMIT       113
#define __FpX_EXTGCD_LIMIT               87
#define __FpX_GCD_LIMIT                  406
#define __FpX_HALFGCD_LIMIT              58
#define __FpX_INVBARRETT_LIMIT           111
#define __FpX_REM_BARRETT_LIMIT          111
#define __Fp_POW_BARRETT_LIMIT           127
#define __Fp_POW_REDC_LIMIT              17
#define __GCD_HALFGCD_LIMIT              -1
#define __HALFGCD_LIMIT                  3
#define __INVMOD_GMP_LIMIT               3
#define __INVNEWTON_LIMIT                75
#define __LOGAGMCX_LIMIT                 22
#define __LOGAGM_LIMIT                   6
#define __MULII_FFT_LIMIT                -1
#define __MULII_KARATSUBA_LIMIT          -1
#define __MULRR_MULII_LIMIT              55
#define __RgX_MUL_LIMIT                  9
#define __RgX_SQR_LIMIT                  38
#define __SQRI_FFT_LIMIT                 -1
#define __SQRI_KARATSUBA_LIMIT           -1
#define __SQRR_SQRI_LIMIT                12
#else
#define __AGM_ATAN_LIMIT                 89
#define __DIVRR_GMP_LIMIT                4
#define __EXPNEWTON_LIMIT                197
#define __EXTGCD_HALFGCD_LIMIT           -1
#define __F2x_MUL_KARATSUBA_LIMIT        13
#define __F2x_MUL_MULII_LIMIT            774
#define __F2xqX_BARRETT_LIMIT            48
#define __F2xqX_DIVREM_BARRETT_LIMIT     127
#define __F2xqX_EXTGCD_LIMIT             127
#define __F2xqX_GCD_LIMIT                884
#define __F2xqX_HALFGCD_LIMIT            89
#define __F2xqX_INVBARRETT_LIMIT         40
#define __F2xqX_REM_BARRETT_LIMIT        127
#define __Flx_BARRETT2_LIMIT             52
#define __Flx_BARRETT_LIMIT              164
#define __Flx_DIVREM2_BARRETT_LIMIT      111
#define __Flx_DIVREM_BARRETT_LIMIT       470
#define __Flx_EXTGCD2_LIMIT              184
#define __Flx_EXTGCD_LIMIT               469
#define __Flx_GCD2_LIMIT                 1281
#define __Flx_GCD_LIMIT                  2817
#define __Flx_HALFGCD2_LIMIT             181
#define __Flx_HALFGCD_LIMIT              586
#define __Flx_INVBARRETT2_LIMIT          397
#define __Flx_INVBARRETT_LIMIT           501
#define __Flx_MUL2_KARATSUBA_LIMIT       9
#define __Flx_MUL2_MULII_LIMIT           8
#define __Flx_MUL_KARATSUBA_LIMIT        57
#define __Flx_MUL_MULII_LIMIT            146
#define __Flx_REM2_BARRETT_LIMIT         89
#define __Flx_REM_BARRETT_LIMIT          388
#define __Flx_SQR2_KARATSUBA_LIMIT       18
#define __Flx_SQR2_SQRI_LIMIT            20
#define __Flx_SQR_KARATSUBA_LIMIT        112
#define __Flx_SQR_SQRI_LIMIT             183
#define __FlxqX_BARRETT_LIMIT            17
#define __FlxqX_DIVREM_BARRETT_LIMIT     46
#define __FlxqX_EXTGCD_LIMIT             44
#define __FlxqX_GCD_LIMIT                1289
#define __FlxqX_HALFGCD_LIMIT            89
#define __FlxqX_INVBARRETT_LIMIT         22
#define __FlxqX_REM_BARRETT_LIMIT        48
#define __FpXQX_BARRETT_LIMIT            12
#define __FpXQX_DIVREM_BARRETT_LIMIT     30
#define __FpXQX_EXTGCD_LIMIT             28
#define __FpXQX_GCD_LIMIT                182
#define __FpXQX_HALFGCD_LIMIT            35
#define __FpXQX_INVBARRETT_LIMIT         40
#define __FpXQX_REM_BARRETT_LIMIT        30
#define __FpX_BARRETT_LIMIT              44
#define __FpX_DIVREM_BARRETT_LIMIT       116
#define __FpX_EXTGCD_LIMIT               81
#define __FpX_GCD_LIMIT                  414
#define __FpX_HALFGCD_LIMIT              55
#define __FpX_INVBARRETT_LIMIT           121
#define __FpX_REM_BARRETT_LIMIT          127
#define __Fp_POW_BARRETT_LIMIT           11
#define __Fp_POW_REDC_LIMIT              3
#define __GCD_HALFGCD_LIMIT              -1
#define __HALFGCD_LIMIT                  22
#define __INVMOD_GMP_LIMIT               3
#define __INVNEWTON_LIMIT                93
#define __LOGAGMCX_LIMIT                 32
#define __LOGAGM_LIMIT                   45
#define __MULII_FFT_LIMIT                -1
#define __MULII_KARATSUBA_LIMIT          -1
#define __MULRR_MULII_LIMIT              19
#define __RgX_MUL_LIMIT                  7
#define __RgX_SQR_LIMIT                  34
#define __SQRI_FFT_LIMIT                 -1
#define __SQRI_KARATSUBA_LIMIT           -1
#define __SQRR_SQRI_LIMIT                9
#endif


/*
 * MAIN FUNCTION
 * C-S-function for Matlab/Simulink
 *
 * derive - function based on algorithms bei C. Join
 *
 * Integration nach Trapezformel
 *
 * (c) Josef Zehetner, 2006-2007
 *
 */

#define S_FUNCTION_LEVEL 2
#define S_FUNCTION_NAME derive

#ifndef MATLAB_MEX_FILE
# include <rti_msg_access.h>
# include <rti_common_msg.h>
#endif

/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
#define NUM_INPUTS          1
/* Input Port  0 */
#define IN_PORT_0_NAME      u0
#define INPUT_0_WIDTH       1
#define INPUT_DIMS_0_COL    1
#define INPUT_0_DTYPE       real_T
#define INPUT_0_COMPLEX     COMPLEX_NO
#define IN_0_FRAME_BASED    FRAME_NO
#define IN_0_DIMS           1-D
#define INPUT_0_FEEDTHROUGH 1
#define IN_0_ISSIGNED        0
#define IN_0_WORDLENGTH      8
#define IN_0_FIXPOINTSCALING 1
#define IN_0_FRACTIONLENGTH  9
#define IN_0_BIAS            0
#define IN_0_SLOPE           0.125

#define NUM_OUTPUTS          1
/* Output Port  0 */
#define OUT_PORT_0_NAME      y0
#define OUTPUT_0_WIDTH       1
#define OUTPUT_DIMS_0_COL    1
#define OUTPUT_0_DTYPE       real_T
#define OUTPUT_0_COMPLEX     COMPLEX_NO
#define OUT_0_FRAME_BASED    FRAME_NO
#define OUT_0_DIMS           1-D
#define OUT_0_ISSIGNED        1
#define OUT_0_WORDLENGTH      8
#define OUT_0_FIXPOINTSCALING 1
#define OUT_0_FRACTIONLENGTH  3
#define OUT_0_BIAS            0
#define OUT_0_SLOPE           0.125

#define NPARAMS              4
/* Parameter  1 */
#define PARAMETER_0_NAME      DERIVATE
#define PARAMETER_0_DTYPE     real_T
#define PARAMETER_0_COMPLEX   COMPLEX_NO
/* Parameter  2 */
#define PARAMETER_1_NAME      FENSTERBREITE
#define PARAMETER_1_DTYPE     real_T
#define PARAMETER_1_COMPLEX   COMPLEX_NO
/* Parameter  3 */
#define PARAMETER_2_NAME      SAMPLETIME
#define PARAMETER_2_DTYPE     real_T
#define PARAMETER_2_COMPLEX   COMPLEX_NO
/* Parameter  4 */
#define PARAMETER_3_NAME      CALC_ORDER
#define PARAMETER_3_DTYPE     real_T
#define PARAMETER_3_COMPLEX   COMPLEX_NO

#define NUM_DISC_STATES      0
#define DISC_STATES_IC       [0]
#define NUM_CONT_STATES      0
#define CONT_STATES_IC       [0]

/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
#include "simstruc.h"
#define PARAM_DEF0(S) ssGetSFcnParam(S, 0)
#define PARAM_DEF1(S) ssGetSFcnParam(S, 1)
#define PARAM_DEF2(S) ssGetSFcnParam(S, 2)
#define PARAM_DEF3(S) ssGetSFcnParam(S, 3)

#define IS_PARAM_DOUBLE(pVal) (mxIsNumeric(pVal) && !mxIsLogical(pVal) &&\
!mxIsEmpty(pVal) && !mxIsSparse(pVal) && !mxIsComplex(pVal) && mxIsDouble(pVal))

// OWN DEFINES
# ifdef MATLAB_MEX_FILE
    #define MAXCOUNT 10001
#else
    #define MAXCOUNT 5001
# endif

#define AKT_POSITION        0

#define _t0 1
#define _t1 _t0*t
#define _t2 _t1*t
#define _t3 _t2*t
#define _t4 _t3*t
#define _t5 _t4*t
#define _t6 _t5*t
#define _t7 _t6*t
#define _t8 _t7*t
#define _t9 _t8*t

#define _Tt0 1
#define _Tt1 _Tt0*(T-t)
#define _Tt2 _Tt1*(T-t)
#define _Tt3 _Tt2*(T-t)
#define _Tt4 _Tt3*(T-t)
#define _Tt5 _Tt4*(T-t)
#define _Tt6 _Tt5*(T-t)
#define _Tt7 _Tt6*(T-t)
#define _Tt8 _Tt7*(T-t)
#define _Tt9 _Tt8*(T-t)

// Parametercheck ein/aus, nur f�r MATLAB_MEX_FILE
#define MDL_CHECK_PARAMETERS

# ifdef MATLAB_MEX_FILE
bool func_init(SimStruct *S)
#else
void func_init(SimStruct *S)
# endif
{

    // allgemeiner INIT Bereich
    int i = 0;

#if defined(MDL_CHECK_PARAMETERS) && defined(MATLAB_MEX_FILE)
    const real_T  *DERIVATE  = mxGetData(PARAM_DEF0(S));
    const real_T  *FENSTERBREITE  = mxGetData(PARAM_DEF1(S));
    const real_T  *SAMPLETIME  = mxGetData(PARAM_DEF2(S));
    const real_T  *CALC_ORDER  = mxGetData(PARAM_DEF3(S));

    int derivate = (int)*DERIVATE;
    double fen = (double)*FENSTERBREITE;
    double sam = (double)*SAMPLETIME;
    int order = (int)*CALC_ORDER;

    if ((fen / sam)>MAXCOUNT)
    {
        ssSetErrorStatus(S,"Fensterbreite/Sampletime �berschreitet maximale Anzahl Schleifendurchl�ufe!");
        return false;
    }
    if ((derivate)>3)
    {
        ssSetErrorStatus(S,"Nur bis zur dritten Ableitung implementiert!");
        return false;
    }
    if ((derivate)<0)
    {
        ssSetErrorStatus(S,"Negative Ableitung nicht sinnvoll!");
        return false;
    }
    switch (derivate)
    {
        case 0:
        {
            if ((order<0) || (order>4))
            {
                ssSetErrorStatus(S,"Nullte Ableitung: Ordnung 0 bis 4 implementiert");
                return false;
            }
        }break;
        case 1:
        {
            if ((order<1) || (order>5))
            {
                ssSetErrorStatus(S,"Erste Ableitung: Ordnung 1 bis 5 implementiert");
                return false;
            }
        }break;
        case 2:
        {
            if ((order<2) || (order>6))
            {
                ssSetErrorStatus(S,"Nullte Ableitung: Ordnung 2 bis 6 implementiert");
                return false;
            }
        }break;
        case 3:
        {
            if ((order<3) || (order>7))
            {
                ssSetErrorStatus(S,"Nullte Ableitung: Ordnung 3 bis 7 implementiert");
                return false;
            }
        }break;
    }
# endif

# ifdef MATLAB_MEX_FILE
    return true;
# endif
}

#if defined(MDL_CHECK_PARAMETERS) && defined(MATLAB_MEX_FILE)
static void mdlCheckParameters(SimStruct *S)
{
    int paramIndex = 0;
    bool validParam = false;
    char paramVector[] ={'1','2','3','4'};
    static char parameterErrorMsg[] ="Parameter   muss DOUBLE sein";
    {
        const mxArray *pVal0 = ssGetSFcnParam(S,0);
        if (!IS_PARAM_DOUBLE(pVal0)) {
            validParam = true;
            paramIndex = 0;
            goto EXIT_POINT;
        }
    }
    {
        const mxArray *pVal1 = ssGetSFcnParam(S,1);
        if (!IS_PARAM_DOUBLE(pVal1)) {
            validParam = true;
            paramIndex = 1;
            goto EXIT_POINT;
        }
    }
    {
        const mxArray *pVal2 = ssGetSFcnParam(S,2);
        if (!IS_PARAM_DOUBLE(pVal2)) {
            validParam = true;
            paramIndex = 2;
            goto EXIT_POINT;
        }
    }
    {
        const mxArray *pVal3 = ssGetSFcnParam(S,3);
        if (!IS_PARAM_DOUBLE(pVal3)) {
            validParam = true;
            paramIndex = 3;
            goto EXIT_POINT;
        }
    }
    EXIT_POINT:
        if (validParam) {
            parameterErrorMsg[11] = paramVector[paramIndex];
            ssSetErrorStatus(S,parameterErrorMsg);
        }
        return;
}
#endif /* MDL_CHECK_PARAMETERS */

static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumSFcnParams(S, NPARAMS);

#   ifdef MATLAB_MEX_FILE
// f�r Simulink
  if (ssGetNumSFcnParams(S) == ssGetSFcnParamsCount(S))
    {
#if defined(MDL_CHECK_PARAMETERS) && defined(MATLAB_MEX_FILE)
        mdlCheckParameters(S);
#endif /* MDL_CHECK_PARAMETERS */
        if (ssGetErrorStatus(S) != NULL)
        {
            return;
        }
  }
    else
    {
     return; /* Parameter mismatch will be reported by Simulink */
    }

    // own inits
    if (!func_init(S))
    {
        return;
    }
#   else

// f�r Realtime Systeme
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S))
    {
        rti_msg_error_set(RTI_SFUNCTION_PARAM_ERROR);
        return;
    }
    func_init(S);

#   endif

    ssSetNumContStates(              S, NUM_CONT_STATES);
    ssSetNumDiscStates(              S, NUM_DISC_STATES);
    ssSetNumInputPorts(              S, NUM_INPUTS);
    ssSetNumOutputPorts(             S, 1);
    ssSetInputPortWidth(             S, 0, INPUT_0_WIDTH);
    ssSetInputPortDataType(S, 0, SS_DOUBLE);
    ssSetInputPortComplexSignal(S, 0, INPUT_0_COMPLEX);
    ssSetInputPortDirectFeedThrough(S, 0, INPUT_0_FEEDTHROUGH);
    ssSetInputPortRequiredContiguous(S, 0, 1); /*direct input signal access*/

    ssSetOutputPortWidth(            S, 0, NUM_OUTPUTS);
    ssSetOutputPortDataType(S, 0, SS_DOUBLE);
    ssSetOutputPortComplexSignal(S, 0, OUTPUT_0_COMPLEX);

    ssSetNumSampleTimes(             S, 1);

    ssSetNumRWork(S, MAXCOUNT);
    ssSetNumIWork(S, 1);
    ssSetNumPWork(S, 0);
    ssSetNumModes(S, 0);
}

static void mdlInitializeSampleTimes(SimStruct *S)
{
    const real_T *SAMPLETIME  = mxGetData(PARAM_DEF2(S));
    ssSetSampleTime(S, 0, *SAMPLETIME);
    ssSetOffsetTime(S, 0, 0.0);
}

// C-modulo operator macht Probleme bei negativen Zahlen
int modulo(int dividend, int divisor)
{
    int res;

    res = dividend % divisor;
    if (res<0)
        res = divisor+res;

    return res;
}

double func_derive0(double runTime, double S, double windowTime, int order)
{
    double f;
    double t;
    double T;

    t = runTime;
    T = windowTime;

    switch(order)
    {
        case 0:
            f = S;
            break;
        case 1:
            f = (4*_Tt1*_t0 - 2*_Tt0*_t1)*S;
            break;
        case 2:
            f = (9*_Tt2*_t0 - 18*_Tt1*_t1 + 3*_Tt0*_t2)*S;
            break;
        case 3:
            f = (16*_Tt3*_t0 - 72*_Tt2*_t1 + 48*_Tt1*_t2 - 4*_Tt0*_t3)*S;
            break;
        case 4:
        default:
            f = (25*_Tt4*_t0 - 200*_Tt3*_t1 + 300*_Tt2*_t2 - 100*_Tt1*_t3 + 5*_Tt0*_t4)*S;
            break;
    }

    return f;
}

double func_derive1(double runTime, double S, double windowTime, int order)
{
    double f;
    double t;
    double T;

    t = runTime;
    T = windowTime;

    switch(order)
    {
        case 0:
            f = 0;
            break;
        case 1:
            f = (-6*_Tt1*_t0 + 6*_Tt0*_t1)*S;
            break;
        case 2:
            f = (-36*_Tt2*_t0 + 120*_Tt1*_t1 - 24*_Tt0*_t2)*S;
            break;
        case 3:
            f = (-120*_Tt3*_t0 + 840*_Tt2*_t1 - 660*_Tt1*_t2 + 60*_Tt0*_t3)*S;
            break;
        case 4:
            f = (-300*_Tt4*_t0 + 3600*_Tt3*_t1 - 6300*_Tt2*_t2 + 2280*_Tt1*_t3 - 120*_Tt0*_t4)*S;
            break;
        case 5:
            f = (-630*_Tt5*_t0 + 11550*_Tt4*_t1 -35700*_Tt3*_t2 + 28980*_Tt2*_t3 - 6090*_Tt1*_t4 + 210*_Tt0*_t5)*S;
            break;
        default:
            f = 0;
            break;
    }

    return f;
}

double func_derive2(double runTime, double S, double windowTime, int order)
{
    double f;
    double t;
    double T;

    t = runTime;
    T = windowTime;

    switch(order)
    {
        case 0:
        case 1:
            f = 0;
            break;
        case 2:
            f = (60*_Tt2*_t0 - 240*_Tt1*_t1 + 60*_Tt0*_t2)*S;
            break;
        case 3:
            f = (480*_Tt3*_t0 - 3960*_Tt2*_t1 + 3600*_Tt1*_t2 - 360*_Tt0*_t3)*S;
            break;
        case 4:
            f = (2100*_Tt4*_t0 - 29400*_Tt3*_t1 + 57960*_Tt2*_t2 - 22680*_Tt1*_t3 + 1260*_Tt0*_t4)*S;
            break;
        case 5:
            f = (6720*_Tt5*_t0 - 142800*_Tt4*_t1 + 490560*_Tt3*_t2 - 426720*_Tt2*_t3 + 94080*_Tt1*_t4 - 3360*_Tt0*_t5)*S;
            break;
        case 6:
            f = (17640*_Tt6*_t0 - 529200*_Tt5*_t1 + 2804760*_Tt4*_t2 - 4304160*_Tt3*_t3 + 2124360*_Tt2*_t4 - 302400*_Tt1*_t5 + 7560*_Tt0*_t6)*S;
            break;
        default:
            f = 0;
            break;
    }

    return f;
}

double func_derive3(double runTime, double S, double windowTime, int order)
{
    double f;
    double t;
    double T;

    t = runTime;
    T = windowTime;

    switch(order)
    {
        case 0:
        case 1:
        case 2:
            f = 0;
            break;
        case 3:
            f = (-840*_Tt3*_t0 + 7560*_Tt2*_t1 - 7560*_Tt1*_t2 + 840*_Tt0*_t3)*S;
            break;
        case 4:
            f = (-8400*_Tt4*_t0 + 127680*_Tt3*_t1 - 272160*_Tt2*_t2 + 114240*_Tt1*_t3 - 6720*_Tt0*_t4)*S;
            break;
        case 5:
            f = (-45360*_Tt5*_t0 + 1043280*_Tt4*_t1 - 3840480*_Tt3*_t2 + 3538080*_Tt2*_t3 - 816480*_Tt1*_t4 + 30240*_Tt0*_t5)*S;
            break;
        case 6:
            f = (-176400*_Tt6*_t0 + 5715360*_Tt5*_t1 - 32281200*_Tt4*_t2 + 52113600*_Tt3*_t3 - 26762400*_Tt2*_t4 + 3931200*_Tt1*_t5 - 100800*_Tt0*_t6)*S;
            break;
        case 7:
            f = (-554400*_Tt7*_t0 + 24060960*_Tt6*_t1 - 193263840*_Tt5*_t2 + 482882400*_Tt4*_t3 - 437698800*_Tt3*_t4 + 143866800*_Tt2*_t5 - 14691600*_Tt1*_t6 + 277200*_Tt0*_t7)*S;
            break;
        default:
            f = 0;
            break;
    }

    return f;
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
    const real_T  *u0  = (const real_T*) ssGetInputPortSignal(S,0);
    real_T        *y0  = (real_T *)ssGetOutputPortRealSignal(S,0);

    double T;
    double time;
    double d;
    int i;
    int j;
    int index1;
    int index2;
    double t1;
    double t2;
    double f1;
    double f2;
    double invT;

    const real_T  *DERIVATE  = mxGetData(PARAM_DEF0(S));
    const real_T  *FENSTERBREITE  = mxGetData(PARAM_DEF1(S));
    const real_T  *SAMPLETIME  = mxGetData(PARAM_DEF2(S));
    const real_T  *CALC_ORDER  = mxGetData(PARAM_DEF3(S));

    int     derivate    = (int)*DERIVATE;
    int     order       = (int)*CALC_ORDER;
    double  fen         = (double)*FENSTERBREITE;
    double  sam         = (double)*SAMPLETIME;

    int     pos         = ssGetIWorkValue(S,AKT_POSITION);
    double *S_vec       = ssGetRWork(S);

    int count = (int)(fen/sam)+1;
    // Fenstergr��e durch maxmiale Speicherbreite beschr�nkt
    if (count>MAXCOUNT)
    {
        count = (int)MAXCOUNT;
        fen = (double)MAXCOUNT*sam;
    }

    time = (double)ssGetT(S);

    // Ringbuffer f�llen
    pos = (++pos) % count;
    ssSetIWorkValue(S,AKT_POSITION,pos);   // aktuelle Ringbufferposition zur�ckschreiben!
    S_vec[pos] = u0[0];

// ab Zeitwert 0+sampletime Berechnung durchf�hren
// Werte zu beginn k�nnen sehr gross werden
#undef COUNT_FROM_START
//#define COUNT_FROM_START

    // Fenster f�r kleine Zeitwerte auf Zeitwert beschr�nken
    // nur sinnvoll, wenn "if (time>=sam)" aktiviert ist
#ifdef COUNT_FROM_START
    if (time>=fen)
    {
        i = count-1;
        T = fen;
    }
    else
    {
        i = (int)(time/sam);
        T = i * sam;
    }
#else
  i = count-1;
    T = fen;
#endif

    d = 0;

#ifdef COUNT_FROM_START
    // ab Zeitwert 0+sampletime Berechnung durchf�hren
    // Werte zu beginn k�nnen sehr gross werden
    if (time>=sam)
#else
    // ab Fensterbreite Berechnung durchf�hren
    // macht mehr Sinn
    if (time>=fen)
#endif
    {
        // Schleife �ber Fensterbreite
        for (j=1;j<i+1;j++)
        {
            t1 = (j-1)*sam;
            t2 = j*sam;

            // modulo-Operation f�r Ringbuffer
            index1 = modulo(pos - (j-1), count);
            index2 = modulo(pos - (j), count);

            switch(derivate)
            {
                //// derive0
                // Funktionswerte
                case 0:
                {
                    f1 = func_derive0(t1,S_vec[index1],T,order);
                    f2 = func_derive0(t2,S_vec[index2],T,order);
                }break;
                //// derive1
                // Funktionswerte
                case 1:
                {
                    f1 = -func_derive1(t1,S_vec[index1],T,order);
                    f2 = -func_derive1(t2,S_vec[index2],T,order);
                }break;
                //// derive2
                // Funktionswerte
                case 2:
                {
                    f1 = func_derive2(t1,S_vec[index1],T,order);
                    f2 = func_derive2(t2,S_vec[index2],T,order);
                }break;
                //// derive3
                // Funktionswerte
                case 3:
                {
                    f1 = -func_derive3(t1,S_vec[index1],T,order);
                    f2 = -func_derive3(t2,S_vec[index2],T,order);
                }break;
                // alle anderen
                default:
                {
                    f1 = 0;
                    f2 = 0;
                }break;
            }

            // Trapezformel
            d = d + (f1+f2)/2.0*sam;

        }

        invT = 1;
        for (i=1; i<=(derivate+order+1); i++)
        {
            invT = invT*T;
        }

        d = 1/invT*d;
    }

    // Ergebnis auf Ausgang schreiben
    y0[0] = d;
}


static void mdlTerminate(SimStruct *S)
{
}

#ifdef  MATLAB_MEX_FILE
# include "simulink.c"
#else
# include "cg_sfun.h"
#endif

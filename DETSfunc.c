/*
 * MAIN FUNCTION
 * C-S-function for Matlab/Simulink
 *
 * DETSFUNC - derivative estimation toolbox
 *
 * based on formulae by J. Reger
 *
 * exact integration
 *
 * (c) Josef Zehetner, 2006-09-28
 *
 */

#define S_FUNCTION_LEVEL 2
#define S_FUNCTION_NAME DETSfunc

#ifndef MATLAB_MEX_FILE
# include <rti_msg_access.h>
# include <rti_common_msg.h>
#endif

/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
#define NUM_INPUTS          2
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

#define IN_PORT_1_NAME      reset
#define INPUT_1_WIDTH       1
#define INPUT_DIMS_1_COL    1
#define INPUT_1_DTYPE       real_T
#define INPUT_1_COMPLEX     COMPLEX_NO
#define IN_1_FRAME_BASED    FRAME_NO
#define IN_1_DIMS           1-D
#define INPUT_1_FEEDTHROUGH 1
#define IN_1_ISSIGNED        0
#define IN_1_WORDLENGTH      8
#define IN_1_FIXPOINTSCALING 1
#define IN_1_FRACTIONLENGTH  9
#define IN_1_BIAS            0
#define IN_1_SLOPE           0.125

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

#define NPARAMS              5
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

#define PARAMETER_4_NAME      NUE
#define PARAMETER_4_DTYPE     real_T
#define PARAMETER_4_COMPLEX   COMPLEX_NO

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
#define PARAM_DEF4(S) ssGetSFcnParam(S, 4)

#define IS_PARAM_DOUBLE(pVal) (mxIsNumeric(pVal) && !mxIsLogical(pVal) &&\
!mxIsEmpty(pVal) && !mxIsSparse(pVal) && !mxIsComplex(pVal) && mxIsDouble(pVal))

// OWN DEFINES
# ifdef MATLAB_MEX_FILE
    #define MAXCOUNT 10001
#else
    #define MAXCOUNT 5001
# endif

#define AKT_POSITION        0
#define N_COUNT             1

// Parametercheck ein/aus, nur f�r MATLAB_MEX_FILE
#define MDL_CHECK_PARAMETERS

int func_init(SimStruct *S)
{

    // allgemeiner INIT Bereich
    int i = 0;

    const real_T  *CALC_ORDER  = mxGetData(PARAM_DEF3(S));
    const real_T  *NUE  = mxGetData(PARAM_DEF4(S));

    int order = (int)*CALC_ORDER;
    int nue = (int)*NUE;

#if defined(MDL_CHECK_PARAMETERS) && defined(MATLAB_MEX_FILE)
    const real_T  *DERIVATE  = mxGetData(PARAM_DEF0(S));
    const real_T  *FENSTERBREITE  = mxGetData(PARAM_DEF1(S));
    const real_T  *SAMPLETIME  = mxGetData(PARAM_DEF2(S));

    int derivate = (int)*DERIVATE;
    double fen = (double)*FENSTERBREITE;
    double sam = (double)*SAMPLETIME;

    if ((fen / sam)>MAXCOUNT)
    {
        ssSetErrorStatus(S,"Fensterbreite/Sampletime �berschreitet maximale Anzahl Schleifendurchl�ufe!");
        return 0;
    }

# endif

    return order+nue+1;
}

#if defined(MDL_CHECK_PARAMETERS) && defined(MATLAB_MEX_FILE)
static void mdlCheckParameters(SimStruct *S)
{
    int paramIndex = 0;
    bool validParam = false;
    char paramVector[] ={'1','2','3','4','5'};
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
    {
        const mxArray *pVal4 = ssGetSFcnParam(S,4);
        if (!IS_PARAM_DOUBLE(pVal4)) {
            validParam = true;
            paramIndex = 4;
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
    int ord_q;

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
    ord_q = func_init(S);
    if (!ord_q)
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
    // own inits
    ord_q = func_init(S);

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

    ssSetInputPortWidth(             S, 1, INPUT_1_WIDTH);
    ssSetInputPortDataType(S, 1, SS_DOUBLE);
    ssSetInputPortComplexSignal(S, 1, INPUT_1_COMPLEX);
    ssSetInputPortDirectFeedThrough(S, 1, INPUT_1_FEEDTHROUGH);
    ssSetInputPortRequiredContiguous(S, 1, 1); /*direct input signal access*/

    ssSetOutputPortWidth(            S, 0, NUM_OUTPUTS);
    ssSetOutputPortDataType(S, 0, SS_DOUBLE);
    ssSetOutputPortComplexSignal(S, 0, OUTPUT_0_COMPLEX);

    ssSetNumSampleTimes(             S, 1);

    ssSetNumRWork(S, 2*MAXCOUNT + ord_q);
    ssSetNumIWork(S, 2);
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

// faktorielle
int fac(int ct)
{
    int res = 1;
    int i;

    for (i=1;i<ct+1;i++)
        res = res*i;

    return res;
}

// power
double pow(double x, int p)
{
    double res = 1.0;
    int i;

    for (i=1;i<p+1;i++)
        res = res*x;

    return res;
}

// Parameter f�r Polynom Q berechnen und Speichern
void func_getQ(SimStruct *S)
{
    int ord_q;

    const real_T  *DERIVATE  = mxGetData(PARAM_DEF0(S));
    const real_T  *FENSTERBREITE  = mxGetData(PARAM_DEF1(S));
    const real_T  *SAMPLETIME  = mxGetData(PARAM_DEF2(S));
    const real_T  *CALC_ORDER  = mxGetData(PARAM_DEF3(S));
    const real_T  *NUE  = mxGetData(PARAM_DEF4(S));

    int j = (int)*DERIVATE;
    double T = (double)*FENSTERBREITE;
    double Ts = (double)*SAMPLETIME;
    int N = (int)*CALC_ORDER;
    int nue = (int)*NUE;

    int h;
    int i;
    int k1; int k2; int k;

    int n = (int)(T/Ts);

    double factor;
    double _kn;
    double _k1n;
    double kn;
    double k1n;

    double *S_vec       = ssGetRWork(S);

    ssSetIWorkValue(S,N_COUNT,n);

    ord_q = N + nue + 1;

    for (i=0;i<=ord_q;i++)
    {
        S_vec[2*MAXCOUNT + i] = 0;   // aktuelle Q Ordnung zur�ckschreiben!
    }

    h = fac(N+j+nue+1)*fac(N+1)*pow(1/T,j);    // f�r alle Werte notwendig

    for (k1=0;k1<=(N-j);k1++)
    {
        for (k2=0;k2<=j;k2++)
        {
            for (i=0;i<=(nue+k1+k2);i++)
            {
                factor = pow(-1,N+i-k1-k2)/(fac(k1)*fac(k2)*fac(i)*fac(N-j-k1)*fac(j-k2)*fac(N-k1-k2)*fac(nue-i+k1+k2)*(N-k1+1)*(N+i-k1-k2+1));
                S_vec[2*MAXCOUNT + N+i-k1-k2+1] += factor;
            }
        }
    }

    for (k=1;k<=n;k++)
    {
        factor = 0.0;
        _kn = ((double)k)/((double)n);
        _k1n = ((double)k-1.0)/((double)n);
        kn = 1;
        k1n = 1;
        for (i=0;i<=ord_q;i++)
        {
            factor += (kn-k1n)*S_vec[2*MAXCOUNT + i];
            kn = kn*_kn;
            k1n = k1n*_k1n;
        }
        S_vec[MAXCOUNT + k] = factor*h;
    }

}



double func_derive(SimStruct *S)
{
    double *S_vec       = ssGetRWork(S);

    int k;
    double d = 0.0;

    int n = ssGetIWorkValue(S,N_COUNT);
    int pos = ssGetIWorkValue(S,AKT_POSITION);
    int index;


    for (k=1;k<=n;k++)
    {
        index = modulo(pos - (k-1), n);
        d += S_vec[index] * S_vec[MAXCOUNT + k];
    }

    return d;
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
    const real_T  *u0  = (const real_T*) ssGetInputPortSignal(S,0);
    const real_T  *reset  = (const real_T*) ssGetInputPortSignal(S,1);
    real_T        *y0  = (real_T *)ssGetOutputPortRealSignal(S,0);

    double time;
    double d;

    int     pos         = ssGetIWorkValue(S,AKT_POSITION);
    double *S_vec       = ssGetRWork(S);

    int n;

    int k;

    const real_T  *FENSTERBREITE  = mxGetData(PARAM_DEF1(S));
    double T = (double)*FENSTERBREITE;
    const real_T  *SAMPLETIME  = mxGetData(PARAM_DEF2(S));
    double Ts = (double)*SAMPLETIME;

    time = (double)ssGetT(S);

    if ((time==0) || (reset[0] > 0.5))
        func_getQ(S);

    n = ssGetIWorkValue(S,N_COUNT);

    // Ringbuffer f�llen
    pos = (++pos) % n;
    ssSetIWorkValue(S,AKT_POSITION,pos);   // aktuelle Ringbufferposition zur�ckschreiben!
    S_vec[pos] = u0[0];

    if (time>T)
    {
        d = func_derive(S);
    }
    else
    {
        d = 0.0;
    }

    k = (int)(time/Ts);
    if (k>MAXCOUNT)
        k=0;

    y0[0] = d;//S_vec[MAXCOUNT + k + 1];
}


static void mdlTerminate(SimStruct *S)
{
}

#ifdef  MATLAB_MEX_FILE
# include "simulink.c"
#else
# include "cg_sfun.h"
#endif

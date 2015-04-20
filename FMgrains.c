#include "ext.h"
#include "ext_obex.h"
#include "z_dsp.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define kTableLength 16384 //Increase frames as a temporal solution for aliasing.

/* FMgrains~ */
/* A simple granular FM oscillator with 2 operators using wavetables.
Carrier Frequency and Modulator Frequency remain constant while Grain Duration is aleatoric.
Darien Brito, 2014.
 */

//-------------------------STRUCT--------------------------//

typedef struct {
    t_pxobject x_obj;
    double increment;       //  Stores increment to read the table
    double increment2;      //  Stores increment to read the modulator table
    double incrementAmp;    //  Stores increment to read the window table
    double index;           //  Index of Wavetable (NOT FM INDEX NUMBER)
    double index2;          //  Index of modWavetable (NOT FM INDEX NUMBER)
    double indexAmp;        //  Index of Window (NOT FM INDEX NUMBER)
    double freq;            //  Freq of the Carrier signal
    double modFreq;         //  Freq of the Modulator Signal
    double modIndex;        //  FM Index value for the modulator
    double sr;              //  Local Sample Rate
    double *waveTable;      //  A pointer to a wavetable
    double *waveTable2;     //  A pointer to the modulator wavetable
    double *window;         //  A pointer to the window wavetable
    double amp;             //  Variable for amplitude
    double min;             //  variable for minimum grain duration
    double max;             //  variable for maximum grain duration
    double maxFreq;         //  variable for minimum carrier frequency
    double minFreq;         //  variable for maximum carrier frequency
    long winFlag;
} fmsynth;//<---------  My struct's name

//---------------------DECLARATIONS--------------------------//

//  This will contain the class
t_class *myClass;
//  Global variables to contain message symbols
//Carrier
t_symbol *sineMess;
t_symbol *squareMess;
t_symbol *sawMess;
t_symbol *phasorMess;
t_symbol *triMess;
//Modulator

t_symbol *sineMessMod;
t_symbol *squareMessMod;
t_symbol *sawMessMod;
t_symbol *phasorMessMod;
t_symbol *triMessMod;

//  Methods
void *FMNew (double v, double a, double m, double d, double min, double max);
void FMFloat(fmsynth *x, double f);
void FMInt(fmsynth *x, long n);
void FMMod(fmsynth *x, double q);
void FMAssist(fmsynth *x, void *b, long m, long a, char *s);
void FMfree(fmsynth *x);
void FMminWin(fmsynth *x, double mn);
void FMmaxWin(fmsynth *x, double mx);
void FMminFreq(fmsynth *x, double mf);
void FMmaxFreq(fmsynth *x, double mxf);
//------------------------------------------
void FillSine (fmsynth *x);
void FillSquare (fmsynth *x);
void FillTri (fmsynth *x);
void FillPhasor (fmsynth *x);
void FillSaw (fmsynth *x);
void FillSineMod (fmsynth *x);
void FillSquareMod (fmsynth *x);
void FillTriMod (fmsynth *x);
void FillPhasorMod (fmsynth *x);
void FillSawMod (fmsynth *x);
//------------------------------------------
void FillWindow (fmsynth *x);
double alea (double min, double max);
//------------------------------------------
void FMAmp(fmsynth *, double amp);
void FMIndex(fmsynth *, double dx);
void FMDsp64(fmsynth *x, t_object *dsp64,short *count, double samplerate, long maxvectorsize, long flags);
void FMPerform64(fmsynth *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam);
void FMSet(fmsynth *x, t_symbol *s);
void FMSetMod(fmsynth *x, t_symbol *s);


//-------------------------MAX MAIN--------------------------//

int C74_EXPORT main(void) {
    t_class *c;             //Variable to store the class
    
    //OBSERVATION: from Max 6.1, the maximun number of arguments that can be declared are 4! If more required, a GIMME signature is adviced
    
    c = class_new("FMgrains~", (method) FMNew, (method) FMfree, sizeof(fmsynth), 0L, A_GIMME, 0);
    class_addmethod(c, (method) FMDsp64, "dsp64", A_CANT, 0);
    class_addmethod(c, (method) FMInt, "int", A_LONG, 0);
    class_addmethod(c, (method) FMFloat, "float", A_FLOAT, 0);
    class_addmethod(c, (method) FMAmp, "ft1", A_FLOAT, 0);
    class_addmethod(c, (method) FMMod, "ft2", A_FLOAT, 0);
    class_addmethod(c, (method) FMIndex, "ft3", A_FLOAT, 0);
    class_addmethod(c, (method) FMminWin, "ft4", A_FLOAT, 0);
    class_addmethod(c, (method) FMmaxWin, "ft5", A_FLOAT, 0);
    class_addmethod(c, (method) FMAssist, "assist", A_CANT, 0);
    class_addmethod(c, (method) FMSet, "set", A_SYM, 0);
    class_dspinit(c);
	class_register(CLASS_BOX, c);
    //Make symbols to assign to global variables for carrier
    sineMess = gensym("sine");
    squareMess = gensym("square");
    sawMess = gensym("saw");
    phasorMess = gensym("phasor");
    triMess = gensym("triangle");
    //Make symbols to assign to global variables for modulator
    sineMessMod = gensym("sinemod");
    squareMessMod = gensym("squaremod");
    sawMessMod = gensym("sawmod");
    phasorMessMod = gensym("phasormod");
    triMessMod = gensym("trianglemod");
    myClass = c;
    return 0;
}


//-----------------------ASSIST METHOD--------------------------//

// assist method (for help strings)
void FMAssist(fmsynth *x, void *b, long m, long a, char *s) {
	if (m == ASSIST_INLET){
        switch (a) {
            case 0: snprintf_zero(s, 256, "(signal/float Carrier Frequency)"); break;
            case 1: snprintf_zero(s, 256, "(float Amp)"); break;
            case 2: snprintf_zero(s, 256, "(float Modulator Frequency)"); break;
            case 3: snprintf_zero(s, 256, "(float Modulator Depth)"); break;
            case 4: snprintf_zero(s, 256, "(float Min grain duration in ms)"); break;
            case 5: snprintf_zero(s, 256, "(float Max grain duration in ms))"); break;
        }
    }
	else
		sprintf(s, "signal/float signal FM-grains");
}

//-----------------------SET METHOD--------------------------//

//CARRIER
//Compare symbol with stored symbols to determine which wavetable should be made

void FMSet(fmsynth *x, t_symbol *s) {
	if (s == sineMess)
		FillSine(x);
	else if (s == squareMess)
		FillSquare(x);
    else if (s == triMess)
        FillTri(x);
    else if (s == sawMess)
        FillSaw(x);
    else if (s == phasorMess)
        FillPhasor(x);
    //MODULATOR
    else if (s == sineMessMod)
		FillSineMod(x);
	else if (s == squareMessMod)
		FillSquareMod(x);
    else if (s == triMessMod)
        FillTriMod(x);
    else if (s == sawMessMod)
        FillSawMod(x);
    else if (s == phasorMessMod)
        FillPhasorMod(x);
}

//---------------------NEW OBJECT METHOD--------------------------//

void *FMNew (double freq, double amp, double modFreq, double dex, double minwin, double maxwin){
    fmsynth *x;
    x = (fmsynth *) object_alloc(myClass);      // Creates an object
    dsp_setup((t_pxobject *)x,0);               // Set up object and specify signal inlets

    floatin((t_object *) x, 5);                 // Inlet for max grain
    floatin((t_object *) x, 4);                 // Inlet for min grain
    floatin((t_object *) x, 3);                 // Inlet for index
    floatin((t_object *) x, 2);                 // Inlet for modFreq
    floatin((t_object *) x, 1);                 // Inlet for Amp
    
    outlet_new((t_pxobject *)x, "signal");      // Add a signal outlet
    x->freq = freq = 220.0;                     // Gimme requires initialization values
    x->amp = amp = 1.0;
    x->modFreq = modFreq = 0.;
    x->modIndex = dex= 1.;
    x->min = minwin = 0.01;
    x->max = maxwin = 0.1;
    x->winFlag = 0;
    // Allocate memory for wavetables:
    x->waveTable = (double *)sysmem_newptr(sizeof(double) * kTableLength);
    x->waveTable2 = (double *)sysmem_newptr(sizeof(double) * kTableLength);
    x->window = (double *)sysmem_newptr(sizeof(double) * kTableLength);
    
    FillSine(x);                                //Default waveform for Carrier is Sine
    FillSineMod(x);                             //Default waveform for Modulator is Sine
    FillWindow(x);                              //Default window is Hamming

    return (x);
}

//-------------------FREE ALLOCATED MEMORY--------------------------//

void FMfree(fmsynth *x) {
    if (x->waveTable) //If there is a pointer to a wavetable
        sysmem_freeptr(x->waveTable); //Free it!
    
    if (x->waveTable2) //If there is a pointer to a wavetable
        sysmem_freeptr(x->waveTable2); //Free it!
    
    if (x->window) //If there is a pointer to a wavetable
        sysmem_freeptr(x->window); //Free it!
    
    dsp_free((t_pxobject *)x); //Free the object with Max's routine
}

//--------------------INLET METHODS--------------------------//

void FMFloat(fmsynth *x, double freq) {
    freq = fabs(freq);
    x->freq = freq;
    x->increment = freq * kTableLength /x->sr;
}

void FMInt(fmsynth *x, long n) {
    n = abs(n);
    x->freq = n;
    x->increment = n * kTableLength /x->sr;
}

void FMAmp(fmsynth *x, double amp) {
    x->amp = amp;
}

void FMMod(fmsynth *x, double modFreq) {
    modFreq = fabs(modFreq);
    x->modFreq = modFreq;
    x->increment2 = modFreq * kTableLength /x->sr;
}

void FMIndex(fmsynth *x, double ix) {
    x->modIndex = ix;
}

void FMminWin(fmsynth *x, double minWin) {
    minWin = fabs(minWin)/1000.;//<---Transform it to milliseconds (MaxMSP convention)
    x->min = minWin;
}

void FMmaxWin(fmsynth *x, double maxWin) {
    maxWin = fabs(maxWin)/1000.;//<---Transform it to milliseconds (MaxMSP convention)
    x->max = maxWin;
}

void FMminFreq(fmsynth *x, double minFreq) {
    minFreq = fabs(minFreq);
    x->minFreq = minFreq;
}

void FMmaxFreq(fmsynth *x, double maxFreq) {
    maxFreq = fabs(maxFreq);
    x->maxFreq = maxFreq;
}


//------------------------64 BIT PROCESS--------------------------//

void FMDsp64(fmsynth *x, t_object *dsp64,short *count, double samplerate, long maxvectorsize, long flags) {
    x->sr = samplerate;
    x->increment2 = x->modFreq * kTableLength/samplerate;
    x->incrementAmp = (1./alea(x->min, x->max)) * kTableLength/ samplerate;
    object_method(dsp64, gensym("dsp_add64"), x, FMPerform64,0, NULL);
}

//------------------------PERFORM METHOD--------------------------//

void FMPerform64 (fmsynth *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam) {
	t_double  *out = outs[0];
	double increment,increment2,incrementAmp,index,*table, *modtable, calc, index2, mod,modindex, sampleRate, *window, indexAmp, winFlag, env,min,max, freq, amp;
    
    //Store values locally for efficiency
    increment = x->increment;
    increment2 = x->increment2;
	index = x->index;
    index2 = x->index2;
	table = x->waveTable;
    modtable = x-> waveTable2;
    freq = x->freq;
    amp = x->amp;
    modindex = x->modIndex;
    sampleRate = x->sr;
    incrementAmp = x->incrementAmp;
    indexAmp = x->indexAmp;
    window = x->window;
    winFlag = x->winFlag;
    min = x->min;
    max = x->max;
	//OUTPUT BUFFER
	while (sampleframes--) {
        
        //SIMPLE FM
        increment = fabs(freq + mod) * kTableLength/sampleRate;
        //CARRIER
        calc = (table[(int) index]); //Carrier
        env =  window[(int) indexAmp];
        //MODULATOR
        mod = (modtable[(int) index2] * fabs(freq * modindex));
        *out++ = (env * calc) * amp;
        //WAVETABLE SCAN
        index += increment; // increment index
        index2 += increment2; // increment index
        indexAmp += incrementAmp;
        //CONDITIONS
		while (index >= kTableLength) // check that increment is within bounds
			index -= kTableLength;
        
        while (index2 >= kTableLength) // check that increment is within bounds
			index2 -= kTableLength;
	
        while (indexAmp >= kTableLength) { // Envelope check
			indexAmp -= kTableLength;
            winFlag = 1;
        }
	}
    
    if (winFlag) { //Reset Envelope
        incrementAmp = (1./alea(min, max)) * kTableLength/ sampleRate;
        winFlag = 0;
    }
    //Save current values in the struct
	x->increment = increment;
    x->increment2 = increment2;
    x->incrementAmp = incrementAmp;
	x->index = index;
    x->index2 = index2;
    x->indexAmp = indexAmp;
    x->amp = amp;
    x->freq = freq;
    x->indexAmp = indexAmp;
    x->min = min;
    x->max = max;
}

//------------------------FILL WAVETABLES--------------------------//

//CARRIER

void FillSine (fmsynth *x) {
    //Fill wavetable with one period of a Sine wave
    int i;
    for(i = 0; i <kTableLength;i++)
        x->waveTable[i] = sin(TWOPI* i/ kTableLength) ;
    
}

void FillSquare (fmsynth *x) {
	//fill waveTable with one period of a Square wave
	int i;
	for (i = 0; i < kTableLength; i++)
		if (i < (kTableLength / 2))
			x->waveTable[i] = 0.99 * 0.5;
		else
			x->waveTable[i] = -0.99 * 0.5;
}

void FillSaw (fmsynth *x) {
    //fill waveTable with one period of a Sawthooth wave
    int i;
    for (i = 0; i <kTableLength; i++){
        x->waveTable[i] = ((2. * i/kTableLength) * 0.5) - 0.5; //Normalize Up
    }
}

void FillPhasor (fmsynth *x) {
    //fill waveTable with one period of a Phasor
    int i;
    for (i = 0; i <kTableLength; i++){
        x->waveTable[i] = (2. * i/kTableLength) * 0.5; //Normalize Up. Unipolar
    }
}

void FillTri (fmsynth *x) {
    //fill waveTable with one period of a Triangle wave
    int i;
    for (i = 0; i < kTableLength; i++){
        x->waveTable[i] = (fabs(((2. * i/kTableLength) * 0.5) - 0.5) * 2.0) - 0.5;
    }
}

//MODULATOR

void FillSineMod (fmsynth *x) {
    //Fill wavetable with one period of a Sine wave
    int i;
    for(i = 0; i <kTableLength;i++)
        x->waveTable2[i] = sin(TWOPI* i/ kTableLength);
    
}


void FillSquareMod (fmsynth *x) {
	//fill waveTable with one period of a Square wave
	int i;
	for (i = 0; i < kTableLength; i++)
		if (i < (kTableLength / 2))
			x->waveTable2[i] = 0.99;
		else
			x->waveTable2[i] = -0.99;
}

void FillSawMod (fmsynth *x) {
    //fill waveTable with one period of a Sawthooth wave
    int i;
    for (i = 0; i <kTableLength; i++){
        x->waveTable2[i] = (2. * i/kTableLength) - 1; //Normalize Up
    }
}

void FillPhasorMod (fmsynth *x) {
    //fill waveTable with one period of a Phasor
    int i;
    for (i = 0; i <kTableLength; i++){
        x->waveTable2[i] = (2. * i/kTableLength) * 0.5; //Normalize Up. Unipolar
    }
}

void FillTriMod (fmsynth *x) {
    //fill waveTable with one period of a Triangle wave
    int i;
    for (i = 0; i < kTableLength; i++){
        x->waveTable2[i] = abs(((2. * i/kTableLength) * 0.5) - 0.5);
    }
}



//HANNING WINDOW

void FillWindow (fmsynth *x) {
    int i;
    for (i = 0; i < kTableLength; i++) {
        x->window[i] = 0.5 * (1 - cos(2*PI*i/kTableLength-1));
    }
}

//RANDOM FUNCTION BETWEEN 2 NUMBERS
double alea (double min, double max) {
    return ((double)rand()/RAND_MAX) * (max - min) + min;
}

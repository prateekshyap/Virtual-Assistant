
/*
///////////////////////////////////////
Generated variables
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
*/

long double *** xi = NULL; //xi in problem-3 solution


long double ** alpha = NULL; //Gets calculated in forward process

long double ** beta = NULL; //Gets calculated in backward process

long double ** delta = NULL; //Gets calculated in Viterbi Algorithm

long double ** gamma = NULL; //Gets calculated in Baum Welch method


/*
///////////////////////////////////////
Model variables
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
*/

long double ** A = NULL; //Transition matrix

long double ** AComplement = NULL; //updated A

long double ** B = NULL; //Probability matrix

long double ** BComplement = NULL; //updated B


long double * Pi = NULL; //Initial probability

long double * PiComplement = NULL; //updated Pi


/*
///////////////////////////////////////
Other variables
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
*/

long double ** codebook = NULL; //codebook

long double pOfOGivenLambda = 0; //probability of an observation sequence given the model lambda

long double pStar = 0; //probability of the state sequence being helpful in modelling

long double pStarComplement = 0; //probability of the state sequence being helpful in modelling for the updated model

long double floorB = 1e-30;


int ** psi = NULL; //Gets generated in Viterbi.h


int * O = NULL; //Observation sequences

int * qStar = NULL; //State sequence

int * qStarComplement = NULL; //State sequence for the updated model

char * resultWord = NULL;


/*
//////////////////////////////////////
Sizes
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
*/

int N = 0; //number of states

int M = 0; //codebook size or number of observations

int T = 100; //size of observation sequence

int R = 0; //read it from file, it will store the number of words currently present in the universe

int duration = 0; //read it from file, it will store the duration of the recording to be considered

int p = 12; //size of Codebook vectors

int universeSize = 0; //size of the universe

/*
///////////////////////////////////////
FILE pointers
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
*/

FILE * AComplementFile = NULL; //new A will be printed

FILE * BComplementFile = NULL; //new B will be printed

FILE * PiComplementFile = NULL; //new Pi will be printed








//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\

/*
///////////////////////////////////
Rejected variables
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
*/

FILE * dataOutputFile = NULL; //output file for required output

FILE * modelOutputFile = NULL; //output file for model
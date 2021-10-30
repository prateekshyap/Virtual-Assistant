int N = 5; //number of states

int M = 32; //codebook size or number of observations

int T = 85; //size of observation sequence

long double ** A = NULL; //Transition matrix

long double ** B = NULL; //Probability matrix

long double ** delta = NULL; //Gets generated in Viterbi Algorithm

long double * Pi = NULL; //Initial probability

long double pStar = 0; //Probability of model

int ** O = NULL; //Observation sequences

int ** psi = NULL; //Gets generated in Viterbi.h

int * qStar = NULL; //State sequence
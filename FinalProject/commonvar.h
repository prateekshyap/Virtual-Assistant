long double *** xi = NULL; //xi in problem-3 solution


long double ** alpha = NULL; //Gets calculated in forward process

long double ** beta = NULL; //Gets calculated in backward process


long double ** A = NULL; //Transition matrix

long double ** AComplement = NULL; //updated A

long double ** B = NULL; //Probability matrix

long double ** BComplement = NULL; //updated B


long double ** delta = NULL; //Gets calculated in Viterbi Algorithm

long double ** gamma = NULL; //Gets calculated in Baum Welch method


long double * Pi = NULL; //Initial probability

long double * PiComplement = NULL; //updated Pi


long double pStar = 0; //Probability of model


int ** O = NULL; //Observation sequences

int ** psi = NULL; //Gets generated in Viterbi.h


int * qStar = NULL; //State sequence


int N = 5; //number of states

int M = 32; //codebook size or number of observations

int T = 85; //size of observation sequence

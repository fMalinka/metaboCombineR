#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <sstream>
#include <limits>
#include <stack>
#include <list>
#include <map>
#include <queue>
#include <string.h>
#include <stdio.h>
#include <iterator>
#include <math.h>
#include <Rcpp.h>

struct elemnt
{
    double score;
    struct elemnt *max_element;
} typedef element;

struct m_elemnt
{
    double score;    
    int index2s1;
    int index2s2;
    std::vector<int> direction;
    std::vector< struct m_elemnt *> ancestor; //in the path
};// typedef matrix_element;

typedef struct
{
    std::vector<int> experimentID;
    Rcpp::NumericMatrix experiments;
    Rcpp::NumericVector mz;
    Rcpp::NumericVector rt;
    std::vector<std::string> mzStr;
    std::vector<std::string> rtStr;
    Rcpp::NumericVector internalStandards_indexes;
    Rcpp::CharacterVector rownames;
    Rcpp::CharacterVector colnames;
    std::vector<int> profile_indexes;
    //std::string experimentName;
} experiment;


typedef struct
{
//    std::vector<struct experiment> experiments;   //1.[] alignemnt, 2.[] vector


    std::list<std::string*> profile_mz;
    std::list<std::string*>::iterator it_profile_mz;//only the same value (mz) is allowed
    std::list<unsigned long> orderPenalty;
    //std::list<double> profile_gap;    //only the same value (mz) is allowed
    //std::list<double> profile_rt;
    //Rcpp::NumericMatrix alignedExp;
    //std::vector<std::vector<double> > alignedExp;
    //std::vector<std::vector<matrix_element> > mymatrix;
    int n_alignments;   //total number of alignments
//    std::vector<int> profile_indexes;
} multipleAlingmentR;


struct mySolution
{
    multipleAlingmentR src;
    std::vector<std::vector<m_elemnt> > mymatrix;
    double score;
    unsigned long startROW;
    unsigned long endROW;

    bool operator<(const mySolution& rhs) const
    {
        return score > rhs.score;
    }
};

struct mySolutionShort
{
    multipleAlingmentR src;
    double *mymatrix;
    double score;
    double maxLocalScore;
    unsigned long startROW;
    unsigned long endROW;
    unsigned long sizeMymatrix;

    bool operator<(const mySolutionShort& rhs) const
    {
        //return score > rhs.score;
        return maxLocalScore > rhs.maxLocalScore;
    }
};

typedef struct
{
    std::vector<int> indexNonref;
    std::vector<int> indexRef;
    experiment *expNonref;
    experiment *expRef;
    int startRow;
} alingmentIndexes;


enum direction {from_none, from_top, from_left, from_oposit};
enum direction_tp {to_none, to_down, to_right, to_diagonal};

typedef std::vector<std::vector<m_elemnt> > myALIGN;

mySolutionShort generateInitionSet(std::vector<std::string> *src, std::vector<std::string> *tgt, std::map<std::string, bool> *commonSeqs, std::multimap<std::string, unsigned long> *peakPosTgz);
std::vector<std::vector<m_elemnt> > generateInitialAlign(multipleAlingmentR *initSetsrc, unsigned long dimY);
std::vector<mySolutionShort> makeCombinations(unsigned long actRow, mySolutionShort *mysol, std::map<std::string, bool> *commonSeqs, std::vector<std::string> *src_profile, std::multimap<std::string, unsigned long> *peakPosTgz);
mySolution computeSimilarityMatrix3R(mySolution *mysol, unsigned long actRow);
std::vector<std::vector<std::vector<int> > > extractAlignmentIndexesR(mySolution *mysol);
void printAlignment(mySolution *mysol, std::vector<std::string> *tgt);
void printAlignmentToFile(mySolution *mysol, std::vector<std::string> *tgt, std::string filename);
void printSemiAlignmentToFile(mySolution *mysol, std::vector<std::string> *tgt, std::string filename);
void printLocalAlignmentToFile(mySolution *mysol, std::vector<std::string> *tgt, std::string filename);
alingmentIndexes processLocalAlignments(mySolution *mysol, std::vector<std::string> *tgt);
void printAlignmentOneColumn(mySolution *mysol, std::vector<std::string> *tgt);
mySolution computeSimilarityMatrix4Complete(std::list<std::string*> src_profile_mz, std::vector<std::string> *tgt_profile);
mySolution computeSimilarityMatrix4CompleteSemiglobal(std::list<std::string *> src_profile_mz, std::vector<std::string> *tgt_profile);
mySolution computeSimilarityMatrix4CompleteLocal(std::list<std::string *> src_profile_mz, std::vector<std::string> *tgt_profile);
mySolutionShort computeSimilarityMatrix4(mySolutionShort mysol, unsigned long actRow, std::vector<std::string> *tgt_profile, double *normalizator);
std::vector<std::string> makeAlignment(std::vector<std::string> sample1, std::vector<std::string> sample2, unsigned long FILTER_LIMIT);
std::vector<std::string> mergeAlignments(mySolution *mysol, std::vector<std::string> *tgt);
void freeMatrix(mySolutionShort *tofree);

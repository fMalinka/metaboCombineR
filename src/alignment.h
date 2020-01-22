/**
    metaboCombineR, alignment.cpp
    Purpose: for alignment representation - header

    @author Frantisek Malinka
    @version 0.99.0 22/01/2020
*/

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

//structure representing the input experiment
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
alingmentIndexes processLocalAlignments(mySolution *mysol, std::vector<std::string> *tgt);
mySolution computeSimilarityMatrix4CompleteSemiglobal(std::list<std::string *> src_profile_mz, std::vector<std::string> *tgt_profile);
mySolution computeSimilarityMatrix4CompleteLocal(std::list<std::string *> src_profile_mz, std::vector<std::string> *tgt_profile);
std::vector<std::string> makeAlignment(std::vector<std::string> sample1, std::vector<std::string> sample2, unsigned long FILTER_LIMIT);

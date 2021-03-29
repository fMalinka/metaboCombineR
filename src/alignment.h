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


/**
 * @brief One element of 2-D matrix in DP
 */
struct m_elemnt
{
    double score;    /**< current DP score        */
    int index2s1;   /**< the element index in 2-D matrix of DP       */
    int index2s2;   /**< the element index in 2-D matrix of DP        */
    std::vector<int> direction; /**< a direction from we get to this element, see enum direction   */
    std::vector< struct m_elemnt *> ancestor; /**< pointer to the previous element in the path   */
};


//structure representing the input experiment
/**
 * @brief Data structure to represent the input experiment
 */
typedef struct
{
    int experimentID;   /**< unique ID        */
    Rcpp::NumericMatrix experiments;    /**< the original input matrix        */
    Rcpp::NumericVector mz; /**< orignal m/z values in numeric rcpp format        */
    Rcpp::NumericVector mz_lower_bound; /**< lower boundary of m/z values according to match mode (ppm, abs, or trunc)        */
    Rcpp::NumericVector mz_upper_bound; /**< lower boundary of m/z values according to match mode (ppm, abs, or trunc)        */
    Rcpp::NumericVector rt; /**< orignal rt values in numeric rcpp format        */
    std::vector<std::string> mzStr; /**< original m/z values in string format        */
    std::vector<std::string> rtStr; /**< orignal rt values in string forma        */
    Rcpp::NumericVector internalStandards_indexes;  /**< current DP score        */
    Rcpp::CharacterVector rownames; /**< rownames of the matrix, i.e. mz and rt        */
    Rcpp::CharacterVector colnames; /**< colnames of the matrix, i.e. sample names        */
    std::vector<int> profile_indexes;   /**< TODO TODO       */
} experiment;


//structure representing kmer of experiment
/**
 * @brief The data structure representing a kmer of experiment
 */
typedef struct
{
    Rcpp::NumericVector mz; /**< orignal m/z values in numeric rcpp format        */
    Rcpp::NumericVector mz_lower_bound; /**< lower boundary of m/z values according to match mode (ppm, abs, or trunc)        */
    Rcpp::NumericVector mz_upper_bound; /**< upper boundary of m/z values according to match mode (ppm, abs, or trunc)        */
    Rcpp::NumericVector rt; /**< orignal rt values in numeric rcpp format        */
    std::vector<std::string> mzStr; /**< original m/z values in string format        */
    std::vector<std::string> rtStr; /**< orignal rt values in string format        */
} experiment_kmer;


/**
 * @brief The element of
 */
typedef struct
{
    std::list<std::string*> profile_mz; /**< list of mz values      */
    std::list<std::string*>::iterator it_profile_mz;    /**< iterator of source multiple alignment       */
    std::list<unsigned long> orderPenalty;
    int n_alignments;   //total number of alignments
} multipleAlingmentR;


/**
 * @brief The struct stores the final multiple alignment
 */
struct mySolution
{
    multipleAlingmentR src; /**< source multiple alignment       */
    std::vector<std::vector<m_elemnt> > mymatrix;   /**< 2-d matrix that is used in DP        */
    double score;   /**< total score        */
    unsigned long startROW; /**< which row is a starting row        */
    unsigned long endROW;   /**< which row is an ending row        */

    bool operator<(const mySolution& rhs) const
    {
        return score > rhs.score;
    }
};


/**
 * @brief Indexes for an alignment
 */
typedef struct
{
    std::vector<int> indexNonref;   /**< vector of indexes for non-reference experiment        */
    std::vector<int> indexRef;  /**< vector of indexes for reference experiment       */
    experiment *expNonref;  /**< non-reference experiment        */
    experiment *expRef; /**< reference experiment        */
    int startRow;   /**< which row is a starting row        */
} alingmentIndexes;


enum direction {from_none, from_top, from_left, from_oposit};
enum direction_tp {to_none, to_down, to_right, to_diagonal};

alingmentIndexes findThebestAlignment(mySolution *mysol, std::vector<std::string> *tgt);
alingmentIndexes processAlignment(mySolution *mysol, std::vector<std::string> *tgt);

mySolution localAlignment(experiment_kmer *kmer, experiment *tgt_profile, double s_match, double s_del, double s_ins);
mySolution globalAlignment(experiment *src, experiment *tgt, double s_match, double s_del, double s_ins);

bool isMzInMatch(double mz1_1, double mz1_2, double mz2_1, double mz2_2);


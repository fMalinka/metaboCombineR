/**
    metaboCombineR, main.cpp
    Purpose: run metaboCombiner

    @author Frantisek Malinka
    @version 0.99.0 22/01/2020
*/
#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>
#include <omp.h>
#include <stdlib.h>
#include <math.h>       /* erfc */
#include <Rcpp.h>
#include "alignment.h"
#include <utility>      // std::pair, std::make_pair
#include <map>
#include <cmath>        // std::abs
#include <omp.h>
#include <limits>


#define SILENT 0

using namespace std;

/**
 * @brief The metaboCombineR class
 */
class metaboCombineR {
public:
    Rcpp::NumericMatrix runKmersAlignment(Rcpp::List allExps, double matchvalue, int windowsize, std::string matchmode, double s_match, double s_del, double s_ins, std::string mzres);
    Rcpp::NumericMatrix runRtcorrectedAlignment(Rcpp::List allExps, double matchvalue, int rtwindowsize, std::string matchmode, double s_match, double s_del, double s_ins, std::string mzres);
    std::vector<int> getAlignmentOrder() {return this->alignmentOrder;}

private:
    std::vector<experiment> experiments; /*!< Vector of the input experiments */
    Rcpp::NumericMatrix finalMatrix;    /*!< Final matrix in rcpp data structure */
    std::vector<int> alignmentOrder;    /*!< The order of alignments in which the final multiple alignment has been created */
    std::string mzres;  /*!< mz values reported in pairwise alignment, ref, nonref, or avg */

    std::vector<std::pair<int,int> > buildAlignment(std::vector<alingmentIndexes> *allAlignments);
    std::vector<std::pair<int,int> > rtCorrection(alingmentIndexes *allAlignments, int rtwindowsize);
    experiment buildMatrix(std::vector<std::pair<int,int> > *index2table, experiment *exp2Ref, experiment *expNonRef);
    experiment buildMatrixRT(std::vector<std::pair<int,int> > *index2table, experiment *exp2Ref, experiment *expNonRef);
    experiment rtcorrectedAlignment(int window_size, double s_match, double s_del, double s_ins, experiment *expNonRef, experiment *exp2Ref);
    double getSimilarity(experiment *first, experiment *second);
    std::pair<unsigned long, unsigned long> findTheMostSimilar(std::vector<experiment> *multiAlignmentExperiment);    
    std::vector<experiment_kmer> createKmers(int window_size, experiment *exp);
    experiment kmersAlignment(int window_size, double s_match, double s_del, double s_ins, experiment *expNonRef, experiment *exp2Ref);
    int updateMzBoundary(experiment *exp, double matchvalue, std::string matchmode);
    std::vector<std::vector<int> > rtCorrectionC(std::vector<std::pair<double, double> > *mzvect1, std::vector<std::pair<double, double> > *mzvect2, int orderLimit, std::map<int, bool> *mergedPeaks);
    std::vector<int> isFeatureInSubset(std::pair<double, double> *elm, std::vector<std::pair<double, double> > *subset);
};


RCPP_MODULE(metaboCombineR){
    using namespace Rcpp;
    class_<metaboCombineR>("metaboCombineR")
    // expose the default constructor
    .constructor()
    .method("runKmersAlignment", &metaboCombineR::runKmersAlignment, "run kmersAlignment algorithm")
    .method("runRtcorrectedAlignment", &metaboCombineR::runRtcorrectedAlignment, "run rtcorrectedAlignment algorithm")
    .method("getAlignmentOrder", &metaboCombineR::getAlignmentOrder, "return the order of alignments")
    ;
}

//doxygen mnual https://cecko.eu/public/doxygen

/**
 * @brief Update m/z lower and upper boundary given by the @a matchmode (ppm, abs, or trunc)
 * @param exp input experiments
 * @param matchvalue a value for mz deviation based on the @a matchmode
 * @param matchmode ppm, abs, or trunc mode
 * @return
 */
int metaboCombineR::updateMzBoundary(experiment *exp, double matchvalue, std::string matchmode)
{
    //ppm mode
    if(matchmode == "ppm")
    {
        for(int irt = 0; irt < exp->rt.size(); ++irt)
        {
            double devppm = (matchvalue * exp->mz(irt))/(1000000.0);
            exp->mz_lower_bound(irt) = exp->mz(irt) - devppm;
            exp->mz_upper_bound(irt) = exp->mz(irt) + devppm;            
        }
    }
    //absolute mode
    else if(matchmode == "abs")
    {
        for(int irt = 0; irt < exp->rt.size(); ++irt)
        {
            exp->mz_lower_bound(irt) = exp->mz(irt) - matchvalue;
            exp->mz_upper_bound(irt) = exp->mz(irt) + matchvalue;
        }
    }
    //truncation mode
    else if(matchmode == "trunc")
    {
        for(int irt = 0; irt < exp->rt.size(); ++irt)
        {     
            exp->mz_lower_bound(irt) = trunc(exp->mz(irt)*pow(10, matchvalue))/pow(10, matchvalue);
            exp->mz_upper_bound(irt) = exp->mz_lower_bound(irt);            
        }
    }
    else   //match value is not specified
    {
        return 1;
    }
    return 0;
}

//subfunction that aggregates two alignment (indexes)
/**
 * @brief Combine individually aligned k-mers into one alignment
 * @param allAlignments vector of k-mer alignments
 * @return vector of indexes pointing to the original experiments
 */
std::vector<std::pair<int,int> > metaboCombineR::buildAlignment(std::vector<alingmentIndexes> *allAlignments)
{

    std::vector<std::pair<int,int> > resultsInds;   //nonref, ref
    std::vector<std::vector<int> > commonIndexesRefCandidates((*allAlignments)[0].expRef->mz.size());
    std::vector<int> commonIndexesRef((*allAlignments)[0].expRef->mz.size(), -1);

    //indexes of matched elements
    for(unsigned long iali = 0; iali < allAlignments->size(); ++iali)
    {
        for(unsigned long irow = 0; irow < (*allAlignments)[iali].indexRef.size(); ++irow)
        {
            //are in a match
            if((*allAlignments)[iali].indexRef[irow] != -1 && (*allAlignments)[iali].indexNonref[irow] != -1)
            {
                commonIndexesRef[(*allAlignments)[iali].indexRef[irow]] = (*allAlignments)[iali].indexNonref[irow]+(*allAlignments)[iali].startRow;
                commonIndexesRefCandidates[(*allAlignments)[iali].indexRef[irow]].push_back((*allAlignments)[iali].indexNonref[irow]+(*allAlignments)[iali].startRow);     
            }

        }
    }

    //proceed unmatched (swapped) elements
    for(unsigned long iref = 0; iref < commonIndexesRefCandidates.size(); ++iref)
    {
        if(!commonIndexesRefCandidates[iref].empty())
        {
            std::map<std::string,std::vector<int> > majorityVote; //mz value, commonIndexesRef
            for(unsigned long icand = 0; icand < commonIndexesRefCandidates[iref].size(); ++icand)
            {
                majorityVote[(*allAlignments)[0].expNonref->rtStr[commonIndexesRefCandidates[iref][icand]]].push_back(commonIndexesRefCandidates[iref][icand]);
            }
            //look for votes
            std::map<std::string,std::vector<int> >::iterator it;

            std::string bestCand;
            unsigned long bestCandvotes = 0;
            double mindistance = DBL_MAX;

            //majority voting
            for (it = majorityVote.begin(); it != majorityVote.end(); it++)
            {
                if(it->second.size() > bestCandvotes)
                {
                    bestCandvotes = it->second.size();
                    bestCand = it->first;
                    mindistance = std::abs((*allAlignments)[0].expRef->rt[iref] - (*allAlignments)[0].expNonref->rt[it->second[0]]);
                }
                else
                {
                    if(it->second.size() == bestCandvotes)
                    {
                        double rtdist = std::abs((*allAlignments)[0].expRef->rt[iref] - (*allAlignments)[0].expNonref->rt[it->second[0]]); //myact2->expNonref->rt[it->second[0]]);
                        if(rtdist < mindistance)
                        {
                            bestCandvotes = it->second.size();
                            bestCand = it->first;
                            mindistance = std::abs((*allAlignments)[0].expRef->rt[iref] - (*allAlignments)[0].expNonref->rt[it->second[0]]);
                        }
                    }
                }

            }
            commonIndexesRef[iref] = majorityVote[bestCand][0];
        }
    }

    //nonreference
    std::vector<int> missingNonref(commonIndexesRef.size(), -1);
    std::map<int, bool> missingNonrefOccurance;
    for(unsigned long ires = 0; ires < commonIndexesRef.size(); ++ires)
    {
        if(commonIndexesRef[ires] != -1)
        {         
            missingNonref[commonIndexesRef[ires]] = 1;
            missingNonrefOccurance[commonIndexesRef[ires]] = true;
        }
    }
    std::vector<int> missingNonrefInd;
    for(unsigned long ic = 0; ic < (*allAlignments)[0].expNonref->mzStr.size(); ++ic)
    {
        if(missingNonrefOccurance.find(ic) == missingNonrefOccurance.end())
        {
            missingNonrefInd.push_back(ic);
        }
    }
    unsigned long mnmislastid = 0;
    for(unsigned long ires = 0; ires < commonIndexesRef.size(); ++ires)
    {
        if(commonIndexesRef[ires] != -1)
        {            
            //only on nonref
            if(!missingNonrefInd.empty())
            {
                if(mnmislastid < (missingNonrefInd.size()))
                {
                    while(missingNonrefInd[mnmislastid] < commonIndexesRef[ires])
                    {

                        resultsInds.push_back(std::pair<int,int>(missingNonrefInd[mnmislastid], -1));                        
                        mnmislastid++;
                        if(mnmislastid >= missingNonrefInd.size())
                            break;
                    }
                }
            }

            //common
            resultsInds.push_back(std::pair<int,int>(commonIndexesRef[ires], ires));
        }
        else
        {
            resultsInds.push_back(std::pair<int,int>(-1, ires));
        }        
    }
    for(unsigned long ires = mnmislastid; ires < missingNonrefInd.size(); ++ires)
    {
        resultsInds.push_back(std::pair<int,int>(missingNonrefInd[ires], -1));    
    }
    return resultsInds;
}


/**
 * @brief Reconstruct two-dimensional matrix from indexes and two experiments
 * @param index2table vector of indexes of final alignment
 * @param exp2Ref reference experiment
 * @param expNonRef non-reference experiment
 * @return final matrix
 */
experiment metaboCombineR::buildMatrix(std::vector<std::pair<int,int> > *index2table, experiment *exp2Ref, experiment *expNonRef)
{
    Rcpp::NumericMatrix mat(index2table->size(), exp2Ref->experiments.ncol() + expNonRef->experiments.ncol());
    Rcpp::CharacterVector chmat(index2table->size());
    Rcpp::CharacterVector chcolnames(exp2Ref->experiments.ncol() + expNonRef->experiments.ncol());
    //add rownames
    std::vector<std::string> newMZstr(index2table->size());
    Rcpp::NumericVector newMZ(index2table->size(), 0);

    //mz values
    for(unsigned long irow = 0; irow < index2table->size(); ++irow)
    {

        //matched
        if((*index2table)[irow].first != -1 && (*index2table)[irow].second != -1)
        {
            //matched - reported mz values
            if(this->mzres == "ref")
            {
                newMZ(irow) = exp2Ref->mz((*index2table)[irow].second);
                newMZstr[irow] = std::to_string(newMZ(irow));
            }
            else if(this->mzres == "nonref")
            {
                newMZ(irow) = expNonRef->mz((*index2table)[irow].first);
                newMZstr[irow] = std::to_string(newMZ(irow));
            }
            //average mz values
            else
            {
                newMZ(irow) = (expNonRef->mz((*index2table)[irow].first) + exp2Ref->mz((*index2table)[irow].second))/2.0;
                newMZstr[irow] = std::to_string(newMZ(irow));
            }
            chmat(irow) = std::string("M") + newMZstr[irow] + std::string("T");            
        }
        else
        {
            if((*index2table)[irow].first < (*index2table)[irow].second)
            {
                chmat(irow) = std::string("M") + exp2Ref->mzStr[(*index2table)[irow].second] + std::string("T");
                newMZstr[irow] = exp2Ref->mzStr[(*index2table)[irow].second];
                newMZ(irow) = exp2Ref->mz((*index2table)[irow].second);
            }
            else
            {
                chmat(irow) = std::string("M") + expNonRef->mzStr[(*index2table)[irow].first] + std::string("T");
                newMZstr[irow] = expNonRef->mzStr[(*index2table)[irow].first];
                newMZ(irow) = expNonRef->mz((*index2table)[irow].first);
            }
        }
    }

    Rcpp::rownames(mat) = chmat;
    //add first nonreftable
    std::vector<std::string> newRTstr(index2table->size());
    Rcpp::NumericVector newRT(index2table->size(), 0);

    //rt values
    for(unsigned long irow = 0; irow < index2table->size(); ++irow)
    {
        for(int icol = 0; icol < expNonRef->experiments.ncol(); ++icol)
        {
            if((*index2table)[irow].first == -1)
                mat(irow, icol) = NA_REAL;
            else
                mat(irow, icol) = expNonRef->experiments((*index2table)[irow].first, icol);
            chcolnames(icol) = expNonRef->colnames(icol);
        }
        if((*index2table)[irow].first != -1)
            chmat(irow) += std::string("/") + expNonRef->rtStr[(*index2table)[irow].first];

        if((*index2table)[irow].first != -1 && (*index2table)[irow].second != -1)
        {
            newRT(irow) = (expNonRef->rt((*index2table)[irow].first) + exp2Ref->rt((*index2table)[irow].second))/2.0;
        }
        else if((*index2table)[irow].first != -1)
        {
            newRT(irow) = expNonRef->rt((*index2table)[irow].first);
        }
        else
        {
            newRT(irow) = exp2Ref->rt((*index2table)[irow].second);
        }
        newRTstr[irow] = std::to_string(newRT(irow));
    }
    //add secondly reftable
    for(unsigned long irow = 0; irow < index2table->size(); ++irow)
    {
        for(int icol = expNonRef->experiments.ncol(); icol < mat.ncol(); ++icol)
        {
            if((*index2table)[irow].second == -1)
                mat(irow, icol) = NA_REAL;
            else
                mat(irow, icol) = exp2Ref->experiments((*index2table)[irow].second, icol-expNonRef->experiments.ncol());
            chcolnames(icol) = exp2Ref->colnames(icol-expNonRef->experiments.ncol());
        }
        if((*index2table)[irow].second != -1)
            chmat(irow) +=  std::string("/") + exp2Ref->rtStr[(*index2table)[irow].second];
    }

    Rcpp::colnames(mat) = chcolnames;
    experiment finalExperiment;
    finalExperiment.experiments = mat;
    finalExperiment.mz = newMZ;
    finalExperiment.mzStr = newMZstr;
    finalExperiment.rownames = chmat;
    finalExperiment.rt = newRT;
    finalExperiment.rtStr = newRTstr;
    finalExperiment.colnames = chcolnames;
    finalExperiment.experimentID = -1;	//minus one -> multiple alingnment
    finalExperiment.mz_lower_bound = Rcpp::NumericVector(newRT.size());
    finalExperiment.mz_upper_bound = Rcpp::NumericVector(newRT.size());
    return finalExperiment;
}


/**
 * @brief Checks whether a feature is in the window (subset) of rtcorrection algorithm
 * @param elm an element to be tested
 * @param subset a vector of elements to be tested for the match
 * @return a vector of indexes where the matches have been found
 */
std::vector<int> metaboCombineR::isFeatureInSubset(std::pair<double, double> *elm, std::vector<std::pair<double, double> > *subset)
{
    std::vector<int> offset; //in match
    for(int ielm = 0; ielm < subset->size(); ++ielm)
    {
        //for n/a value continue
        if(std::isnan((*subset)[ielm].first) || std::isnan((*subset)[ielm].second))
        {
            continue;
        }
        if(isMzInMatch(elm->first, elm->second, (*subset)[ielm].first, (*subset)[ielm].second))
        {
            offset.push_back(ielm);
        }        
    }    
    return offset;
}



/**
 * @brief retention time correction function
 * @param mzvect1 the first vector of m/z values
 * @param mzvect2 the second vector of m/z values
 * @param orderLimit a size of the window
 * @param mergedPeaks peaks that are merged
 * @return positions of features to be swapped/merged
 */
std::vector<std::vector<int> > metaboCombineR::rtCorrectionC(std::vector<std::pair<double, double> > *mzvect1, std::vector<std::pair<double, double> > *mzvect2, int orderLimit,  std::map<int, bool> *mergedPeaks)
{
  std::vector<int> swapIndex1;
  std::vector<int> swapIndex2;

  //from two vectors of alignment to one vector
  std::vector<std::pair<double, double> > commonFeatures(mzvect1->size());
  for(int i = 0; i < mzvect1->size(); ++i)
  {
      bool naVect1 = std::isnan((*mzvect1)[i].first) || std::isnan((*mzvect1)[i].second);
      if(naVect1)
        commonFeatures[i] = std::pair<double, double>((*mzvect2)[i].first, (*mzvect2)[i].second);
      else
        commonFeatures[i] = std::pair<double, double>((*mzvect1)[i].first, (*mzvect1)[i].second);
  }
  //no swaps are allowed
  if(orderLimit == 0)
  {
      std::vector<std::vector<int> > swapInd;
      swapInd.push_back(std::vector<int>());
      swapInd.push_back(std::vector<int>());
      return swapInd;
  }
  //the window size of one, faster version
  else if(orderLimit == 1)
  {
      //iterate over alignment
      for(int i = 0; i < mzvect1->size()-1; ++i)
      {
        bool naVect1 = std::isnan((*mzvect1)[i].first) || std::isnan((*mzvect1)[i].second);
        bool naVect2 = std::isnan((*mzvect2)[i].first) || std::isnan((*mzvect2)[i].second);
        if(!(naVect1 && naVect2))
        {
          //vect2 is non-empty
          if(naVect1)
          {
            //the window position boundary
            int lowerind = (i+1);
            std::vector<std::pair<double, double> > subset;
            subset.push_back((*mzvect1)[lowerind]);
            std::vector<int> matchedIndex  = isFeatureInSubset(&(*mzvect2)[i], &subset);
            if(matchedIndex.size() > 0)
            {
                (*mzvect1)[1 + i] = std::pair<double, double>(std::numeric_limits<double>::quiet_NaN(),
                                                                      std::numeric_limits<double>::quiet_NaN());
                swapIndex2.push_back(1 + i);
                swapIndex1.push_back(i);
                (*mergedPeaks)[1 + i] = true;
            }
          }
          else  //vect1 is empty
          {
            //the window position boundary
            int lowerind = (i+1);
            std::vector<std::pair<double, double> > subset;
            subset.push_back((*mzvect2)[lowerind]);
            std::vector<int> matchedIndex = isFeatureInSubset(&(*mzvect1)[i], &subset);
            if(matchedIndex.size() > 0)
            {
                (*mzvect2)[1 + i] = std::pair<double, double>(std::numeric_limits<double>::quiet_NaN(),
                                                                          std::numeric_limits<double>::quiet_NaN());
                swapIndex2.push_back(1 + i);
                swapIndex1.push_back(i);
                (*mergedPeaks)[1 + i] = true;
            }
          }
        }
      }
      std::vector<std::vector<int> > swapInd;
      swapInd.push_back(swapIndex1);
      swapInd.push_back(swapIndex2);
      return swapInd;
  }
  //the window size more than 1
  else
  {
      //iterate over alignment
      for(int i = 0; i < mzvect1->size()-1; ++i)
      {
        bool naVect1 = std::isnan((*mzvect1)[i].first) || std::isnan((*mzvect1)[i].second);
        bool naVect2 = std::isnan((*mzvect2)[i].first) || std::isnan((*mzvect2)[i].second);

        if(!(naVect1 && naVect2))
        {
          //vect2 is non-empty
          if(naVect1)
          {
            //the window position boundary
            int lowerind = (i+1);
            int upperind = std::min((lowerind+orderLimit), static_cast<int>(mzvect1->size()));
            //a subset of features
            std::vector<std::pair<double, double> >::const_iterator subsetLower = commonFeatures.begin() + lowerind;
            std::vector<std::pair<double, double> >::const_iterator subsetUpper = commonFeatures.begin() + upperind;
            std::vector<std::pair<double, double> > subset(subsetLower, subsetUpper);
            //is the element in the subset?
            std::vector<int> matchedIndex = isFeatureInSubset(&(*mzvect2)[i], &subset);
            bool firstOccur = false;
            //iterate over matched and merge them
            for(int imatchedIndex = 0; imatchedIndex < matchedIndex.size(); ++imatchedIndex)
            {
                (*mzvect1)[matchedIndex[imatchedIndex] + 1 + i] = std::pair<double, double>(std::numeric_limits<double>::quiet_NaN(),
                                                                  std::numeric_limits<double>::quiet_NaN());

                (*mzvect2)[matchedIndex[imatchedIndex] + 1 + i] = std::pair<double, double>(std::numeric_limits<double>::quiet_NaN(),
                                                                    std::numeric_limits<double>::quiet_NaN());
                commonFeatures[matchedIndex[imatchedIndex] + 1 + i] = std::pair<double, double>(std::numeric_limits<double>::quiet_NaN(),
                                                                             std::numeric_limits<double>::quiet_NaN());
                if(firstOccur == false)
                {
                    swapIndex2.push_back(matchedIndex[imatchedIndex] + 1 + i);
                    swapIndex1.push_back(i);
                    firstOccur = true;
                }                
                (*mergedPeaks)[matchedIndex[imatchedIndex] + 1 + i] = true;

            }
          }
          else //vect1 is empty
          {
            //the window position boundary
            int lowerind = (i+1);
            int upperind = std::min((lowerind+orderLimit), static_cast<int>(mzvect2->size()));
            //a subset of features
            std::vector<std::pair<double, double> >::const_iterator subsetLower = commonFeatures.begin() + lowerind;
            std::vector<std::pair<double, double> >::const_iterator subsetUpper = commonFeatures.begin() + upperind;
            std::vector<std::pair<double, double> > subset(subsetLower, subsetUpper);
            //is the element in the subset?
            std::vector<int> matchedIndex = isFeatureInSubset(&(*mzvect1)[i], &subset);
            bool firstOccur = false;
            //iterate over matched and merge them
            for(int imatchedIndex = 0; imatchedIndex < matchedIndex.size(); ++imatchedIndex)
            {
                (*mzvect2)[matchedIndex[imatchedIndex] + 1 + i] = std::pair<double, double>(std::numeric_limits<double>::quiet_NaN(),
                                                                          std::numeric_limits<double>::quiet_NaN());
                (*mzvect1)[matchedIndex[imatchedIndex] + 1 + i] = std::pair<double, double>(std::numeric_limits<double>::quiet_NaN(),
                                                                          std::numeric_limits<double>::quiet_NaN());
                if(firstOccur == false)
                {
                    swapIndex2.push_back(matchedIndex[imatchedIndex] + 1 + i);
                    swapIndex1.push_back(i);
                    firstOccur = true;                    
                }
                (*mergedPeaks)[matchedIndex[imatchedIndex] + 1 + i] = true;
                commonFeatures[matchedIndex[imatchedIndex] + 1 + i] = std::pair<double, double>(std::numeric_limits<double>::quiet_NaN(),
                                                                                 std::numeric_limits<double>::quiet_NaN());                
            }
          }
        }
      }
      std::vector<std::vector<int> > swapInd;
      swapInd.push_back(swapIndex1);
      swapInd.push_back(swapIndex2);
      return swapInd;
  }
  std::vector<std::vector<int> > swapInd;
  swapInd.push_back(std::vector<int>());
  swapInd.push_back(std::vector<int>());
  return swapInd;
}



//subfunction that aggregates two alignment (indexes)
/**
 * @brief rt correction of swapped features for pairwise experiments
 * @param allAlignments an alignment indexes of reference and non-reference experiment
 * @param rtwindowsize a size of the window
 * @return positions of features to be swapped
 */
std::vector<std::pair<int,int> > metaboCombineR::rtCorrection(alingmentIndexes *allAlignments, int rtwindowsize)
{

    std::vector<std::pair<int,int> > index; //index as a reference for the original table
    //mz boundary for an alignment
    std::vector<std::pair<double, double> > alingMzRef(allAlignments->indexRef.size());
    std::vector<std::pair<double, double> > alingMzNonRef(allAlignments->indexNonref.size());

    //mz boundaries
    for(unsigned long ial = 0; ial < allAlignments->indexRef.size(); ++ial)
    {
        if(allAlignments->indexRef[ial] != -1)
        {                        
            alingMzRef[ial] = std::pair<double, double>(allAlignments->expRef->mz_lower_bound[allAlignments->indexRef[ial]],
                    allAlignments->expRef->mz_upper_bound[allAlignments->indexRef[ial]]);
        }
        else
        {
            alingMzRef[ial] = std::pair<double, double>(std::numeric_limits<double>::quiet_NaN(),
                                                std::numeric_limits<double>::quiet_NaN());
        }        
    }
    for(unsigned long ial = 0; ial < allAlignments->indexNonref.size(); ++ial)
    {
        if(allAlignments->indexNonref[ial] != -1)
        {                  
            alingMzNonRef[ial] = std::pair<double, double>(allAlignments->expNonref->mz_lower_bound[allAlignments->indexNonref[ial]],
                    allAlignments->expNonref->mz_upper_bound[allAlignments->indexNonref[ial]]);
        }
        else
        {
            alingMzNonRef[ial] = std::pair<double, double>(std::numeric_limits<double>::quiet_NaN(),
                                                std::numeric_limits<double>::quiet_NaN());
        }
    }

    std::map<int, bool> mergedPeaks;
    //indicate swapped peaks
    std::vector<std::vector<int> > swapInd = rtCorrectionC(&alingMzRef, &alingMzNonRef, rtwindowsize, &mergedPeaks);
    std::vector<int> mz1swappedInd = swapInd[0];
    std::vector<int> mz2swappedInd = swapInd[1];

    //common elements
    for(unsigned long icom = 0; icom < allAlignments->indexRef.size(); ++icom)
    {
        //avoid already merged features
        if (mergedPeaks.find(icom) != mergedPeaks.end())
        {
            std::pair<int,int> mypair;
            mypair.second = -1;
            mypair.first = -1;
            index.push_back(mypair);
        }
        else
        {
            //in match
            if(allAlignments->indexRef[icom] != -1 && allAlignments->indexNonref[icom] != -1)
            {
                std::pair<int,int> mypair;
                mypair.second = allAlignments->indexRef[icom];
                mypair.first = allAlignments->indexNonref[icom];
                index.push_back(mypair);
            }
            //for non-reference
            else if(allAlignments->indexRef[icom] == -1)
            {
                std::pair<int,int> mypair;
                mypair.second = -1;
                mypair.first = allAlignments->indexNonref[icom];
                index.push_back(mypair);
            }
            //for reference
            else if(allAlignments->indexNonref[icom] == -1)
            {
                std::pair<int,int> mypair;
                mypair.second = allAlignments->indexRef[icom];
                mypair.first = -1;
                index.push_back(mypair);
            }
            else
            {
                std::pair<int,int> mypair;
                mypair.second = -1;
                mypair.first = -1;
                index.push_back(mypair);
            }
        }

    }

    //change indexes for swapped features
    for(int iswap = 0; iswap < mz1swappedInd.size(); ++iswap)
    {

            if(index[mz1swappedInd[iswap]].first == -1)
            {
                index[mz1swappedInd[iswap]].first = allAlignments->indexNonref[mz2swappedInd[iswap]];
            }
            else
            {
                index[mz1swappedInd[iswap]].second = allAlignments->indexRef[mz2swappedInd[iswap]];
            }

    }

    //the final vector format
    std::vector<std::pair<int,int> > finalindex;
    for(unsigned long ii = 0; ii < index.size(); ++ii)
    {
        finalindex.push_back(index[ii]);
    }
    return finalindex;

}


/**
 * @brief create a vector of k-mers from an input experiment
 * @param window_size a k parameter of k-mer
 * @param exp the input experiment
 * @return a vector of k-mers in the experiment data structure
 */
std::vector<experiment_kmer> metaboCombineR::createKmers(int window_size, experiment *exp)
{
    std::vector<experiment_kmer> kmerExperiments;
    std::vector<std::vector<std::string > > mysamples_windows;
    for(int istart = window_size-1; istart < exp->mz.size(); ++istart)
    {
        std::vector<std::string > actSample;
        actSample.resize(window_size);

        Rcpp::NumericVector mz(window_size);
        Rcpp::NumericVector mz_lower_bound(window_size);
        Rcpp::NumericVector mz_upper_bound(window_size);
        std::vector<std::string> mzStr(window_size);

        for(int iwin = 0; iwin < window_size; ++iwin)
        {
            actSample[(window_size-1)-iwin] = exp->mzStr[istart-iwin];
            mz((window_size-1)-iwin) =  exp->mz(istart-iwin);
            mz_lower_bound((window_size-1)-iwin) =  exp->mz_lower_bound(istart-iwin);
            mz_upper_bound((window_size-1)-iwin) =  exp->mz_upper_bound(istart-iwin);
            mzStr[(window_size-1)-iwin] =  exp->mzStr[istart-iwin];
        }
        experiment_kmer kmer;
        kmer.mz = mz;
        kmer.mzStr = mzStr;
        kmer.mz_lower_bound = mz_lower_bound;
        kmer.mz_upper_bound = mz_upper_bound;
        kmerExperiments.push_back(kmer);
        mysamples_windows.push_back(actSample);
    }
    return kmerExperiments;
}

/**
 * @brief Reconstruct two-dimensional matrix from indexes and two experiments for rtCorrectedAlignment algorithm
 * @param index2table vector of indexes of final alignment
 * @param exp2Ref reference experiment
 * @param expNonRef non-reference experiment
 * @return final matrix
 */
experiment metaboCombineR::buildMatrixRT(std::vector<std::pair<int,int> > *index2table, experiment *exp2Ref, experiment *expNonRef)
{
    Rcpp::NumericMatrix mat(index2table->size(), exp2Ref->experiments.ncol() + expNonRef->experiments.ncol());
    Rcpp::CharacterVector chmat(index2table->size());
    Rcpp::CharacterVector chcolnames(exp2Ref->experiments.ncol() + expNonRef->experiments.ncol());
    //add rownames
    std::vector<std::string> newMZstr(index2table->size());
    Rcpp::NumericVector newMZ(index2table->size(), 0);

    for(unsigned long irow = 0; irow < index2table->size(); ++irow)
    {
        if((*index2table)[irow].first == -1 && (*index2table)[irow].second == -1)
            continue;
        if((*index2table)[irow].first != -1 && (*index2table)[irow].second != -1)
        {
            //matched, reported mz values
            if(this->mzres == "ref")
            {
                newMZ(irow) = (exp2Ref->mz((*index2table)[irow].first));
                newMZstr[irow] = exp2Ref->mzStr[(*index2table)[irow].first];
            }
            else if(this->mzres == "nonref")
            {
                newMZ(irow) = (expNonRef->mz((*index2table)[irow].first));
                newMZstr[irow] = expNonRef->mzStr[(*index2table)[irow].first];
            }
            //average mz values
            else
            {
                newMZ(irow) = (expNonRef->mz((*index2table)[irow].first) + exp2Ref->mz((*index2table)[irow].second))/2.0;
                newMZstr[irow] = std::to_string(newMZ(irow));
            }
            chmat(irow) = std::string("M") + newMZstr[irow] + std::string("T");            
        }
        else
        {
            //reference experiment
            if((*index2table)[irow].first < (*index2table)[irow].second)
            {
                chmat(irow) = std::string("M") + exp2Ref->mzStr[(*index2table)[irow].second] + std::string("T");
                newMZstr[irow] = exp2Ref->mzStr[(*index2table)[irow].second];
                newMZ(irow) = exp2Ref->mz((*index2table)[irow].second);
            }
            //non-reference experiment
            else
            {
                chmat(irow) = std::string("M") + expNonRef->mzStr[(*index2table)[irow].first] + std::string("T");
                newMZstr[irow] = expNonRef->mzStr[(*index2table)[irow].first];
                newMZ(irow) = expNonRef->mz((*index2table)[irow].first);
            }
        }

    }
    Rcpp::rownames(mat) = chmat;
    //add first nonreftable
    std::vector<std::string> newRTstr(index2table->size());
    Rcpp::NumericVector newRT(index2table->size(), 0);


    for(unsigned long irow = 0; irow < index2table->size(); ++irow)
    {
        if((*index2table)[irow].first == -1 && (*index2table)[irow].second == -1)
            continue;
        for(int icol = 0; icol < expNonRef->experiments.ncol(); ++icol)
        {
            if((*index2table)[irow].first == -1)
                mat(irow, icol) = NA_REAL;
            else
                mat(irow, icol) = expNonRef->experiments((*index2table)[irow].first, icol);
            chcolnames(icol) = expNonRef->colnames(icol);
        }
        if((*index2table)[irow].first != -1)
            chmat(irow) += std::string("/") + expNonRef->rtStr[(*index2table)[irow].first];

        if((*index2table)[irow].first != -1 && (*index2table)[irow].second != -1)
        {
            newRT(irow) = (expNonRef->rt((*index2table)[irow].first) + exp2Ref->rt((*index2table)[irow].second))/2.0;
        }
        else if((*index2table)[irow].first != -1)
        {
            newRT(irow) = expNonRef->rt((*index2table)[irow].first);
        }
        else
        {
            newRT(irow) = exp2Ref->rt((*index2table)[irow].second);
        }
        newRTstr[irow] = std::to_string(newRT(irow));
    }
    //add secondly reftable
    int matnewsize = 0;
    for(unsigned long irow = 0; irow < index2table->size(); ++irow)
    {
        if((*index2table)[irow].first == -1 && (*index2table)[irow].second == -1)
            continue;
        ++matnewsize;

        for(int icol = expNonRef->experiments.ncol(); icol < mat.ncol(); ++icol)
        {
            if((*index2table)[irow].second == -1)
                mat(irow, icol) = NA_REAL;
            else
                mat(irow, icol) = exp2Ref->experiments((*index2table)[irow].second, icol-expNonRef->experiments.ncol());
            chcolnames(icol) = exp2Ref->colnames(icol-expNonRef->experiments.ncol());
        }
        if((*index2table)[irow].second != -1)
            chmat(irow) +=  std::string("/") + exp2Ref->rtStr[(*index2table)[irow].second];
    }

     Rcpp::NumericMatrix mat2(matnewsize, exp2Ref->experiments.ncol() + expNonRef->experiments.ncol());
     Rcpp::CharacterVector chmat2(matnewsize);
     std::vector<std::string> newRTstr2(matnewsize);
     Rcpp::NumericVector newRT2(matnewsize);
     std::vector<std::string> newMZstr2(matnewsize);
     Rcpp::NumericVector newMZ2(matnewsize);

     int inewrow = 0;
     for(unsigned long irow = 0; irow < index2table->size(); ++irow)
     {
         if((*index2table)[irow].first != -1 || (*index2table)[irow].second != -1)
         {
             chmat2[inewrow] = (chmat(irow));
             mat2(inewrow, Rcpp::_) = mat(irow, Rcpp::_);
             newRTstr2[inewrow] = newRTstr[irow];
             newRT2[inewrow] = newRT2[irow];
             newMZstr2[inewrow] = newMZstr[irow];
             newMZ2[inewrow] = newMZ[irow];
            ++inewrow;
         }
     }
     Rcpp::colnames(mat2) = chcolnames;
     Rcpp::rownames(mat2) = chmat2;
     experiment toReturn;
     toReturn.experiments = mat2;
     toReturn.mz = newMZ2;
     toReturn.mzStr = newMZstr2;
     toReturn.rownames = chmat2;
     toReturn.rt = newRT2;
     toReturn.rtStr = newRTstr2;
     toReturn.colnames = chcolnames;
     toReturn.experimentID = -1;	//minus one -> multiple alingnment
     toReturn.mz_lower_bound = Rcpp::NumericVector(newRT.size());
     toReturn.mz_upper_bound = Rcpp::NumericVector(newRT.size());

    return toReturn;
}


//functiont that combines two alignments into the final experiment -RT
/**
 * @brief rtcorrectedAlignment for pairwise alignment
 * @param window_size a value of the window size
 * @param s_match a score for matched features
 * @param s_del a penalization for deletion of features
 * @param s_ins a penalization for insertion of features
 * @param expNonRef the non-reference experiment
 * @param exp2Ref the reference experiment
 * @return the final pairwise aligned experiment
 */
experiment metaboCombineR::rtcorrectedAlignment(int window_size, double s_match, double s_del, double s_ins, experiment *expNonRef, experiment *exp2Ref)
{
    //the algorithm as is presented in the paper
    mySolution rowAlign = globalAlignment(expNonRef, exp2Ref, s_match, s_del, s_ins);
    alingmentIndexes resAlign = processAlignment(&rowAlign, &(exp2Ref->mzStr));
    resAlign.expNonref = expNonRef;
    resAlign.expRef = exp2Ref;
    resAlign.startRow = 0;

    std::vector<std::pair<int,int> > index2table = rtCorrection(&resAlign, window_size); //nonref, ref
    experiment finalMatrix = buildMatrixRT(&index2table, exp2Ref, expNonRef);
    return finalMatrix;
}


//compute similarity function
/**
 * @brief compute a similarity of two experiments
 * @param first the first experiment
 * @param second the second expriment
 * @return similarity score
 */
double metaboCombineR::getSimilarity(experiment *first, experiment *second)
{
    std::vector<int> hist1(first->mzStr.size(), 0);
    std::vector<int> hist2(second->mzStr.size(), 0);

    //first experiment
    for(int ifirst = 0; ifirst < first->mzStr.size(); ++ifirst)
    {
        //iterate in the triangular format, without reflexive rel.
        for(int ifirst2 = ifirst + 1; ifirst2 < first->mzStr.size(); ++ifirst2)
        {
            //mzs are matched
            if(isMzInMatch(first->mz_lower_bound[ifirst], first->mz_upper_bound[ifirst], first->mz_lower_bound[ifirst2], first->mz_upper_bound[ifirst2]))
            {
                hist1[ifirst]++;
                hist1[ifirst2]++;
            }
        }
    }

    //second experiment
    for(int isecond = 0; isecond < second->mzStr.size(); ++isecond)
    {
        //iterate in the triangular format, without reflexive rel.
        for(int isecond2 = isecond + 1; isecond2 < second->mzStr.size(); ++isecond2)
        {
            //mzs are matched
            if(isMzInMatch(second->mz_lower_bound[isecond], second->mz_upper_bound[isecond], second->mz_lower_bound[isecond2], second->mz_upper_bound[isecond2]))
            {
                hist2[isecond]++;
                hist2[isecond2]++;
            }
        }
    }

    //first vs second experiment
    double subscore = 0;
    double intersectionCount = 0;
    for(int ifirst = 0; ifirst < first->mzStr.size(); ++ifirst)
    {
        //iterate in the triangular format, without reflexive rel.
        int occur1 = 0;
        int occur2 = 0;
        for(int isecond = 0; isecond < second->mzStr.size(); ++isecond)
        {
            //mzs are matched
            if(isMzInMatch(first->mz_lower_bound[ifirst], first->mz_upper_bound[ifirst], second->mz_lower_bound[isecond], second->mz_upper_bound[isecond]))
            {
                occur1 += hist1[ifirst];
                occur2 += hist2[isecond];
                ++intersectionCount;
            }
        }
        subscore += std::abs(occur1 - occur2);
    }
    double score = first->mzStr.size() + second->mzStr.size() - (2*intersectionCount);
    score +=  subscore;
    return score;
}


//compute the similarity for each experiment in vector and returns positions of the best
/**
 * @brief find the most similar experiment according to getSimilarity function
 * @param multiAlignmentExperiment a vector of experiments
 * @return positons of the most similar experiments
 */
std::pair<unsigned long, unsigned long> metaboCombineR::findTheMostSimilar(std::vector<experiment> *multiAlignmentExperiment)
{
    double bestscore = DBL_MAX;
    int indNonREF = 0;
    int indREF = 0;
    for(unsigned long irow = 0; irow < multiAlignmentExperiment->size(); ++irow)
    {
        for(unsigned long icol = irow+1; icol < multiAlignmentExperiment->size(); ++icol)
        {
            double score = getSimilarity(&((*multiAlignmentExperiment)[irow]), &((*multiAlignmentExperiment)[icol]));            
            if(score < bestscore)
            {
                bestscore = score;
                if((*multiAlignmentExperiment)[irow].mzStr.size() < (*multiAlignmentExperiment)[icol].mzStr.size())
                {
                    indNonREF = irow;
                    indREF = icol;
                }
                else
                {
                    indNonREF = icol;
                    indREF = irow;
                }
            }
        }
    }
    std::pair<int, int> toreturn;
    toreturn.first = indNonREF;
    toreturn.second = indREF;
    return toreturn;
}


/**
 * @brief kmersAlignment for pairwise alignment
 * @param window_size a k parameter of k-mer
 * @param s_match a score for matched features
 * @param s_del a penalization for deletion of features
 * @param s_ins a penalization for insertion of features
 * @param expNonRef a non-reference experiment
 * @param exp2Ref a reference experiment
 * @return the final matrix of pairwise alignment
 */
experiment metaboCombineR::kmersAlignment(int window_size, double s_match, double s_del, double s_ins, experiment *expNonRef, experiment *exp2Ref)
{
    //the algorithm as is presented in the paper
    std::vector<experiment_kmer> allkmers = createKmers(window_size, expNonRef);
    std::vector<alingmentIndexes> allAlignments;
    for(unsigned long isample = 0; isample < allkmers.size(); ++isample)
    {
        mySolution rowAlign = localAlignment(&allkmers[isample], exp2Ref, s_match, s_del, s_ins);
        alingmentIndexes resAlign = findThebestAlignment(&rowAlign, &(exp2Ref->mzStr));
        resAlign.expNonref = expNonRef;
        resAlign.expRef = exp2Ref;
        resAlign.startRow = isample;
        allAlignments.push_back(resAlign);
    }
    std::vector<std::pair<int,int> > index2table = buildAlignment(&allAlignments); //combine various alignment into one

    experiment finalMatrix = buildMatrix(&index2table, exp2Ref, expNonRef);
    return finalMatrix;
}


//main function, this  combines all experiments
/**
 * @brief runKmersAlignment the core function
 * @param allExps experiments to be processed
 * @param matchvalue a value for the feature matching boundaries
 * @param windowsize a k parameter of k-mers algorithm
 * @param matchmode type of m/z value matching (ppm, abs, or trunc)
 * @param s_match a score for matched features
 * @param s_del a penalization for deletion of features
 * @param s_ins a penalization for insertion of features
 * @param mzres mz values reported in pairwise alignment
 * @return the final multi-alignment numeric matrix
 */
Rcpp::NumericMatrix metaboCombineR::runKmersAlignment(Rcpp::List allExps, double matchvalue, int windowsize, std::string matchmode, double s_match, double s_del, double s_ins, std::string mzres)
{
    this->mzres = mzres;
    // Make function callable from C++
    Rcpp::Function getRTsR("getRTs");
    Rcpp::Function getMZsR("getMZs");
    Rcpp::Function sortDataFrameByRt("sortDataFrameByRt");

    //iterate over all experiments and prepare data structures
    for(int iexp = 0; iexp < allExps.size(); ++iexp)
    {
        experiment newExp;        
        Rcpp::NumericMatrix originalDataFrame = sortDataFrameByRt(allExps[iexp]);
        Rcpp::NumericMatrix originalFrame = originalDataFrame;

        newExp.rt = getRTsR(originalDataFrame);
        newExp.mz = getMZsR(originalDataFrame);

        newExp.mz_lower_bound = Rcpp::NumericVector(newExp.mz.size());
        newExp.mz_upper_bound = Rcpp::NumericVector(newExp.mz.size());

        newExp.experimentID = iexp;
        newExp.mzStr.resize(newExp.rt.size());
        newExp.rtStr.resize(newExp.rt.size());
        newExp.profile_indexes.resize(newExp.rt.size());


        for(int irt = 0; irt < newExp.rt.size(); ++irt)
        {
            newExp.profile_indexes[irt] = irt;
            newExp.mzStr[irt] = std::to_string(newExp.mz(irt));
            newExp.rtStr[irt] = std::to_string(newExp.rt(irt));
        }
        if(updateMzBoundary(&newExp, matchvalue, matchmode))
        {
            Rcpp::Rcerr << "match values is not specified. Terminating...";
            return Rcpp::NumericMatrix();
        }

        newExp.experiments = originalFrame;
        Rcpp::CharacterVector dataColnames = Rcpp::colnames(originalDataFrame);
        newExp.colnames = dataColnames;
        for(int icolname = 0; icolname < dataColnames.size(); ++icolname)
        {
            std::ostringstream expname;
            expname << "subExp" << iexp << "_" << dataColnames(icolname);
            newExp.colnames(icolname) = expname.str();
        }
        this->experiments.push_back(newExp);
    }

    std::vector<experiment> multiAlignmentExperiment;
    Rcpp::Rcout << "Successfully loaded " << this->experiments.size()  << " datasets!" << std::endl;

    //push all experiments
    for(int iexp = 0; iexp < allExps.size(); ++iexp)
    {
        multiAlignmentExperiment.push_back(this->experiments[iexp]);
    }
    int window_size = windowsize;

    //the main loop, each iteration makes one pairwise alignment
    Rcpp::Rcout << "Working .";
    do
    {
        std::pair<unsigned long, unsigned long> mostsimilar = findTheMostSimilar(&multiAlignmentExperiment);
        Rcpp::Rcout << ".";        
        std::vector<experiment> newExp;
        for(unsigned long iel = 0; iel < multiAlignmentExperiment.size(); ++iel)
        {
            if(iel != mostsimilar.first && iel != mostsimilar.second)
            {
                newExp.push_back(multiAlignmentExperiment[iel]);
            }
        }
        //make a pairwise alingment using kmersAlignment algorithm
        experiment newalignment = kmersAlignment(window_size, s_match, s_del, s_ins, &(multiAlignmentExperiment[mostsimilar.first]), &(multiAlignmentExperiment[mostsimilar.second]));
        //update mz boundary
        if(updateMzBoundary(&newalignment, matchvalue, matchmode))
        {
            Rcpp::Rcerr << "match values is not specified. Terminating...";
            return Rcpp::NumericMatrix();
        }
        newExp.push_back(newalignment);
        
        //store the order of experiments in which are aligned
        if(multiAlignmentExperiment[mostsimilar.first].experimentID >= 0)
        {
            this->alignmentOrder.push_back(multiAlignmentExperiment[mostsimilar.first].experimentID);
        }
        if(multiAlignmentExperiment[mostsimilar.second].experimentID >= 0)
        {
            this->alignmentOrder.push_back(multiAlignmentExperiment[mostsimilar.second].experimentID);
        }
        multiAlignmentExperiment = newExp;        
    } while(multiAlignmentExperiment.size() > 1);

    Rcpp::Rcout << " Done" << std::endl;
    return multiAlignmentExperiment.back().experiments;
}


//main function, this  combines all experiments
/**
 * @brief runRtcorrectedAlignment the core function
 * @param allExps experiments to be processed
 * @param matchvalue a value for the feature matching boundaries
 * @param rtwindowsize a window size for the rt correction
 * @param matchmode type of m/z value matching (ppm, abs, or trunc)
 * @param s_match a score for matched features
 * @param s_del a penalization for deletion of features
 * @param s_ins a penalization for insertion of features
 * @param mzres mz values reported in pairwise alignment
 * @return the final multi-alignment numeric matrix
 */
Rcpp::NumericMatrix metaboCombineR::runRtcorrectedAlignment(Rcpp::List allExps, double matchvalue, int rtwindowsize, std::string matchmode, double s_match, double s_del, double s_ins, std::string mzres)
{
    this->mzres = mzres;
    // Make function callable from C++
    Rcpp::Function getRTsR("getRTs");
    Rcpp::Function getMZsR("getMZs");
    Rcpp::Function sortDataFrameByRt("sortDataFrameByRt");

    //iterate over all experiments and prepare data structures
    for(int iexp = 0; iexp < allExps.size(); ++iexp)
    {
        experiment newExp;        
        Rcpp::NumericMatrix originalDataFrame = sortDataFrameByRt(allExps[iexp]);
        Rcpp::NumericMatrix originalFrame = originalDataFrame;

        newExp.rt = getRTsR(originalDataFrame);
        newExp.mz = getMZsR(originalDataFrame);

        newExp.mz_lower_bound = Rcpp::NumericVector(newExp.mz.size());
        newExp.mz_upper_bound = Rcpp::NumericVector(newExp.mz.size());

        newExp.experimentID = iexp;
        newExp.mzStr.resize(newExp.rt.size());
        newExp.rtStr.resize(newExp.rt.size());
        newExp.profile_indexes.resize(newExp.rt.size());
        for(int irt = 0; irt < newExp.rt.size(); ++irt)
        {
            newExp.profile_indexes[irt] = irt;
            newExp.mzStr[irt] = std::to_string(newExp.mz(irt));
            newExp.rtStr[irt] = std::to_string(newExp.rt(irt));

        }
        //update mz boundary
        if(updateMzBoundary(&newExp, matchvalue, matchmode))
        {
            Rcpp::Rcerr << "match values is not specified. Terminating...";
            return Rcpp::NumericMatrix();
        }

        newExp.experiments = originalFrame;
        Rcpp::CharacterVector dataColnames = Rcpp::colnames(originalDataFrame);
        newExp.colnames = dataColnames;
        for(int icolname = 0; icolname < dataColnames.size(); ++icolname)
        {
            std::ostringstream expname;
            expname << "subExp" << iexp << "_" << dataColnames(icolname);
            newExp.colnames(icolname) = expname.str();
        }
        this->experiments.push_back(newExp);
    }

    std::vector<experiment> multiAlignmentExperiment;
    Rcpp::Rcout << "Successfully loaded " << this->experiments.size()  << " datasets!" << std::endl;

    //push all experiments
    for(int iexp = 0; iexp < allExps.size(); ++iexp)
    {
        multiAlignmentExperiment.push_back(this->experiments[iexp]);
    }
    int window_size = rtwindowsize;

    //the main loop, each iteration makes one pairwise alignment
    Rcpp::Rcout << "Working .";
    do
    {
        std::pair<unsigned long, unsigned long> mostsimilar = findTheMostSimilar(&multiAlignmentExperiment);
        Rcpp::Rcout << ".";
        std::vector<experiment> newExp;
        for(unsigned long iel = 0; iel < multiAlignmentExperiment.size(); ++iel)
        {
            if(iel != mostsimilar.first && iel != mostsimilar.second)
            {
                newExp.push_back(multiAlignmentExperiment[iel]);
            }
        }
        experiment newalignment;
        //reference sequence will be the longer sequence
        if(multiAlignmentExperiment[mostsimilar.first].mzStr.size() < multiAlignmentExperiment[mostsimilar.second].mzStr.size())
        {
            //make a pairwise alingment using rtcorrectedAlignment algorithm
            newalignment = rtcorrectedAlignment(window_size, s_match, s_del, s_ins, &(multiAlignmentExperiment[mostsimilar.first]), &(multiAlignmentExperiment[mostsimilar.second]));
        }
        else
        {
            //make a pairwise alingment using rtcorrectedAlignment algorithm
            newalignment = rtcorrectedAlignment(window_size, s_match, s_del, s_ins, &(multiAlignmentExperiment[mostsimilar.second]), &(multiAlignmentExperiment[mostsimilar.first]));
        }

        if(updateMzBoundary(&newalignment, matchvalue, matchmode))
        {
            Rcpp::Rcerr << "Match values are not specified. Terminating...";
            return Rcpp::NumericMatrix();
        }
        newExp.push_back(newalignment);

        //store the order in which experiments are aligned
        if(multiAlignmentExperiment[mostsimilar.first].experimentID >= 0)
        {
            this->alignmentOrder.push_back(multiAlignmentExperiment[mostsimilar.first].experimentID);
        }
        if(multiAlignmentExperiment[mostsimilar.second].experimentID >= 0)
        {
            this->alignmentOrder.push_back(multiAlignmentExperiment[mostsimilar.second].experimentID);
        }
        multiAlignmentExperiment = newExp;
    } while(multiAlignmentExperiment.size() > 1);

    Rcpp::Rcout << " Done" << std::endl;
    return multiAlignmentExperiment.back().experiments;
}

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


#define SILENT 0

typedef enum {MATCH, ONLYREF, ONLYNONREF} alignIndexSTATES;

using namespace std;

class metaboCombineR {
public:
    Rcpp::NumericMatrix run(Rcpp::List allExps, int mzprecision, int windowsize);
    Rcpp::NumericMatrix runRT(Rcpp::List allExps, int mzprecision, int rtwindowsize);
    void combineExperiments();
    Rcpp::NumericMatrix getfinalMatrix();

private:
    std::vector<experiment> experiments;
    Rcpp::NumericMatrix finalMatrix;
    std::vector<std::pair<int,int> > aggregateAlignments(std::vector<alingmentIndexes> *allAlignments);
    std::vector<std::pair<int,int> > rtcorrection(alingmentIndexes *allAlignments, int rtwindowsize);
    experiment generateAlignmentExperiment(int window_size, experiment *expNonRef, experiment *exp2Ref);
    experiment generateAlignmentExperimentRT(int window_size, experiment *expNonRef, experiment *exp2Ref);
    double getSimilarity(experiment *first, experiment *second);
    std::pair<int, int> findTheMostSimilar(std::vector<experiment> *multiAlignmentExperiment);
};


RCPP_MODULE(metaboCombineR){
    using namespace Rcpp;
    class_<metaboCombineR>("metaboCombineR")
    // expose the default constructor
    .constructor()
    .method("run", &metaboCombineR::run, "get the message")
    .method("runRT", &metaboCombineR::runRT, "get the message")
    ;
}

//subfunction that aggregates two alignment (indexes)
std::vector<std::pair<int,int> > metaboCombineR::aggregateAlignments(std::vector<alingmentIndexes> *allAlignments)
{

    std::vector<std::pair<int,int> > resultsInds;   //nonref, ref
    std::vector<std::vector<int> > commonIndexesRefCandidates((*allAlignments)[0].expRef->mz.size());
    std::vector<int> commonIndexesRef((*allAlignments)[0].expRef->mz.size(), -1);

    for(int iali = 0; iali < allAlignments->size(); ++iali)
    {
        for(int irow = 0; irow < (*allAlignments)[iali].indexRef.size(); ++irow)
        {
            if((*allAlignments)[iali].indexRef[irow] != -1 && (*allAlignments)[iali].indexNonref[irow] != -1)
            {
                commonIndexesRef[(*allAlignments)[iali].indexRef[irow]] = (*allAlignments)[iali].indexNonref[irow]+(*allAlignments)[iali].startRow;
                commonIndexesRefCandidates[(*allAlignments)[iali].indexRef[irow]].push_back((*allAlignments)[iali].indexNonref[irow]+(*allAlignments)[iali].startRow);     
            }

        }
    }
    for(int iref = 0; iref < commonIndexesRefCandidates.size(); ++iref)
    {
        if(!commonIndexesRefCandidates[iref].empty())
        {
            std::map<std::string,std::vector<int> > majorityVote; //mz value, commonIndexesRef
            for(int icand = 0; icand < commonIndexesRefCandidates[iref].size(); ++icand)
            {
                majorityVote[(*allAlignments)[0].expNonref->rtStr[commonIndexesRefCandidates[iref][icand]]].push_back(commonIndexesRefCandidates[iref][icand]);
            }
            //look for votes
            std::map<std::string,std::vector<int> >::iterator it;

            std::string bestCand;
            int bestCandvotes = 0;
            double mindistance = 100000;
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
    std::vector<int> missingNonref(commonIndexesRef.size(), -1);
    std::map<int, bool> missingNonrefOccurance;
    for(int ires = 0; ires < commonIndexesRef.size(); ++ires)
    {
        if(commonIndexesRef[ires] != -1)
        {         
            missingNonref[commonIndexesRef[ires]] = 1;
            missingNonrefOccurance[commonIndexesRef[ires]] = true;
        }
    }

    std::vector<int> missingNonrefInd;
    for(int ic = 0; ic < (*allAlignments)[0].expNonref->mzStr.size(); ++ic)
    {
        if(missingNonrefOccurance.find(ic) == missingNonrefOccurance.end())
        {
            missingNonrefInd.push_back(ic);
        }
    }
    int mnmislastid = 0;//0;
    for(int ires = 0; ires < commonIndexesRef.size(); ++ires)
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
    for(int ires = mnmislastid; ires < missingNonrefInd.size(); ++ires)
    {
        resultsInds.push_back(std::pair<int,int>(missingNonrefInd[ires], -1));    
    }
    return resultsInds;
}

//subfunction that aggregates two alignment (indexes)
std::vector<std::pair<int,int> > metaboCombineR::rtcorrection(alingmentIndexes *allAlignments, int rtwindowsize)
{

    std::vector<std::pair<int,int> > index;

    Rcpp::Rcout << "rtcorrection 1#" << std::endl;
    std::vector<std::string> onerowalign(allAlignments->indexRef.size());
    std::vector<std::pair<int,int>> pairindex(allAlignments->indexRef.size());
    Rcpp::Rcout << "inderef: " << allAlignments->indexRef.size() << " nonref: " << allAlignments->indexNonref.size() << std::endl;

    std::vector<std::string> a1(allAlignments->indexRef.size());
    std::vector<std::string> a2(allAlignments->indexNonref.size());

    int realindex1 = 0;
    int realindex2 = 0;
    std::vector<int> reals1index(allAlignments->indexRef.size());
    std::vector<int> reals2index(allAlignments->indexRef.size());

    for(int ial = 0; ial < allAlignments->indexRef.size(); ++ial)
    {
        if(allAlignments->indexRef[ial] != -1)
        {
            onerowalign[ial] = allAlignments->expRef->mzStr[allAlignments->indexRef[ial]];
            a1[ial] = allAlignments->expRef->mzStr[allAlignments->indexRef[ial]];
            reals1index[ial] = realindex1;
            ++realindex1;
            //Rcpp::Rcout << allAlignments->expRef->mzStr[allAlignments->indexRef[ial]] << " ";
        }
        else
        {
            a1[ial] = "-1";
            onerowalign[ial] = allAlignments->expNonref->mzStr[allAlignments->indexNonref[ial]];
            //reals1index[ial] = -1;
            ;//Rcpp::Rcout << "-1 ";
        }        
    }

    for(int ial = 0; ial < allAlignments->indexNonref.size(); ++ial)
    {
        if(allAlignments->indexNonref[ial] != -1)
        {
            a2[ial] = allAlignments->expNonref->mzStr[allAlignments->indexNonref[ial]];
            reals2index[ial] = realindex2;
            ++realindex2;
        }
        else
        {
            //reals2index[ial] = -1;
            a2[ial] = "-1";
        }
    }




    Rcpp::Rcout << "size onerowaling: " << onerowalign.size() << std::endl;

    for(int iii = 0; iii < onerowalign.size(); ++iii)
    {
        Rcpp::Rcout << onerowalign[iii] << std::endl;
    }
    Rcpp::Rcout << std::endl;


   // Rcpp::Rcout << "AAAAAAAAAAAAAAAAA" << std::endl;
   /* for(int ial = 0; ial < allAlignments->indexNonref.size(); ++ial)
    {
        if(allAlignments->indexNonref[ial] != -1)
        {
            //onerowalign[ial] = allAlignments->expRef->mzStr[allAlignments->indexRef[ial]];
            ;//Rcpp::Rcout << allAlignments->expNonref->mzStr[allAlignments->indexNonref[ial]] << " ";
        }
        else
        {
            ;//Rcpp::Rcout << "-1 ";
        }
       // Rcpp::Rcout << allAlignments->indexNonref[ial] << " ";
       // Rcpp::Rcout << allAlignments->expRef->mzStr[allAlignments->indexRef[ial]] << " ";
    }
    */
    Rcpp::Rcout << std::endl;
   /*
    for(int ial = 0; ial < a1.size(); ++ial)
    {
        Rcpp::Rcout << a1[ial] << " ";
    }
    Rcpp::Rcout << std::endl;
    for(int ial = 0; ial < a2.size(); ++ial)
    {
        Rcpp::Rcout << a2[ial] << " ";
    }

*/
    Rcpp::Rcout << std::endl;

    Rcpp::Function correctClassiqAling("correctClassiqAling");
    Rcpp::List listRTcorrected = correctClassiqAling(onerowalign, rtwindowsize);
    Rcpp::Rcout << "rtcorrection 2#" << std::endl;
    Rcpp::CharacterVector xxx0 = listRTcorrected[0];
    Rcpp::CharacterVector xxx1 = listRTcorrected[1];
    for(int ii = 0; ii < xxx0.size(); ++ii)
    {
        Rcpp::Rcout << "(" << xxx1[ii] << "," << xxx0[ii] << ") ";
    }

    Rcpp::Rcout << "rtcorrection 4#" << std::endl;

    Rcpp::Rcout << "correcting verssion 2" << std::endl;
    Rcpp::Function correctClassiqAling2("correctClassiqAling2");
    Rcpp::List listRTcorrected2 = correctClassiqAling2(a1, a2, rtwindowsize);
    Rcpp::CharacterVector mz1 = listRTcorrected2[0];
    Rcpp::CharacterVector mz2 = listRTcorrected2[1];

    //std::vector<std::string> reals1index;
    //std::vector<std::string> reals2index;
    for(int ii = 0; ii < mz1.size(); ++ii)
    {
        Rcpp::Rcout << "(" << mz1[ii] << "," << mz2[ii] << " # " << reals1index[atoi(mz1[ii])-1] << "," << reals2index[atoi(mz2[ii])-1] << ") ";
    }

    //common elements
    for(int icom = 0; icom < allAlignments->indexRef.size(); ++icom)
    {
        if(allAlignments->indexRef[icom] != -1 && allAlignments->indexNonref[icom] != -1)
        {
            std::pair<int,int> mypair;
            //mypair.first = allAlignments->indexRef[icom];
            //mypair.second = allAlignments->indexNonref[icom];
            mypair.second = allAlignments->indexRef[icom];
            mypair.first = allAlignments->indexNonref[icom];
            index.push_back(mypair);
            //Rcpp::Rcout << "first: " << mypair.first << " texxt: " << allAlignments->expRef->mzStr[mypair.first] <<  " second: " << mypair.second << " text: " << allAlignments->expNonref->mzStr[mypair.second] << std::endl;
        }
        else if(allAlignments->indexRef[icom] == -1)
        {
            std::pair<int,int> mypair;
            mypair.second = -1;//allAlignments->indexRef[icom];
            mypair.first = allAlignments->indexNonref[icom];
            index.push_back(mypair);
        }
        else
        {
            std::pair<int,int> mypair;
            mypair.second = allAlignments->indexRef[icom];
            mypair.first = -1;// allAlignments->indexNonref[icom];
            index.push_back(mypair);
            //Rcpp::Rcout << "n first: " <<allAlignments->indexRef[icom] << " nsecond: " << allAlignments->indexNonref[icom] << std::endl;
        }
    }

    int iindexLast = 0;
    for(int iswap = 0; iswap < mz1.size(); ++iswap)
    {
        //mz1 ref
        //Rcpp::Rcout << "test: " << index[iswap].second << " test2: " << index[iswap].first << std::endl; //ref
        //Rcpp::Rcout << "ref swap: " << mz1[iswap] << " nonrefswap: " << mz2[iswap]; //ref
        //Rcpp::Rcout << "z ref swap: " << index[reals1index[atoi(mz1[iswap])-1]].first << " bude: " << index[reals2index[atoi(mz2[iswap])-1]].first << std::endl;
        Rcpp::Rcout << "index1 real: " << index[atoi(mz1[iswap])-1].second << " index2 real: " << index[atoi(mz1[iswap])-1].second << " index ref: " << index[reals1index[atoi(mz1[iswap])-1]].first << " a " << index[reals1index[atoi(mz1[iswap])-1]].second << " | " << index[reals2index[atoi(mz2[iswap])-1]].first << " a " << index[reals2index[atoi(mz2[iswap])-1]].second  << std::endl;

         index[atoi(mz1[iswap])-1].first = index[atoi(mz2[iswap])-1].first;
         index[atoi(mz2[iswap])-1].first = -1;
         Rcpp::Rcout << "bla: " << index[atoi(mz2[iswap])-1].first << " bla2: " << index[atoi(mz2[iswap])-1].second << std::endl;
        //Rcpp::Rcout << "reference: " << index[atoi(mz1[iswap])-1].first << " a " << index[atoi(mz1[iswap])-1].second << " non-reference: " << index[atoi(mz2[iswap])-1].first << " a " << index[atoi(mz2[iswap])-1].second << std::endl;

        //reference


        //reals1index[atoi(mz1[iswap])-1];
        //nonreference
        //reals2index[atoi(mz2[iswap])-1];

        //index[reals1index[atoi(mz1[iswap])-1]].first = index[reals2index[atoi(mz2[iswap])-1]].first;
        //index[reals2index[atoi(mz2[iswap])-1]].first = -1;
        //index[atoi(mz2[iswap])].first = -1;
        //mz2[iswap]; //nonref
    }

    std::vector<std::pair<int,int> > finalindex;
    for(int ii = 0; ii < index.size(); ++ii)
    {
        Rcpp::Rcout << "A: " << index[ii].first << " B: " << index[ii].second << std::endl;
        //if(index[ii].first != -1 && index[ii].second != -1)
        //{
            finalindex.push_back(index[ii]);
        //}
        //if(index[ii].first != -1 && index[ii].first ==
        ;//index[ii].second;
    }
    return finalindex;
    //            std::vector<std::pair<int,int> > resultsInds;   //nonref, ref
    //            return resultsInds;


    //alignment index to sequence index


    /*
    std::vector<std::pair<int,int> > resultsInds;   //nonref, ref
    std::vector<std::vector<int> > commonIndexesRefCandidates((*allAlignments)[0].expRef->mz.size());
    std::vector<int> commonIndexesRef((*allAlignments)[0].expRef->mz.size(), -1);

    for(int iali = 0; iali < allAlignments->size(); ++iali)
    {
        for(int irow = 0; irow < (*allAlignments)[iali].indexRef.size(); ++irow)
        {
            if((*allAlignments)[iali].indexRef[irow] != -1 && (*allAlignments)[iali].indexNonref[irow] != -1)
            {
                commonIndexesRef[(*allAlignments)[iali].indexRef[irow]] = (*allAlignments)[iali].indexNonref[irow]+(*allAlignments)[iali].startRow;
                commonIndexesRefCandidates[(*allAlignments)[iali].indexRef[irow]].push_back((*allAlignments)[iali].indexNonref[irow]+(*allAlignments)[iali].startRow);
            }

        }
    }
    for(int iref = 0; iref < commonIndexesRefCandidates.size(); ++iref)
    {
        if(!commonIndexesRefCandidates[iref].empty())
        {
            std::map<std::string,std::vector<int> > majorityVote; //mz value, commonIndexesRef
            for(int icand = 0; icand < commonIndexesRefCandidates[iref].size(); ++icand)
            {
                majorityVote[(*allAlignments)[0].expNonref->rtStr[commonIndexesRefCandidates[iref][icand]]].push_back(commonIndexesRefCandidates[iref][icand]);
            }
            //look for votes
            std::map<std::string,std::vector<int> >::iterator it;

            std::string bestCand;
            int bestCandvotes = 0;
            double mindistance = 100000;
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
    std::vector<int> missingNonref(commonIndexesRef.size(), -1);
    std::map<int, bool> missingNonrefOccurance;
    for(int ires = 0; ires < commonIndexesRef.size(); ++ires)
    {
        if(commonIndexesRef[ires] != -1)
        {
            missingNonref[commonIndexesRef[ires]] = 1;
            missingNonrefOccurance[commonIndexesRef[ires]] = true;
        }
    }

    std::vector<int> missingNonrefInd;
    for(int ic = 0; ic < (*allAlignments)[0].expNonref->mzStr.size(); ++ic)
    {
        if(missingNonrefOccurance.find(ic) == missingNonrefOccurance.end())
        {
            missingNonrefInd.push_back(ic);
        }
    }
    int mnmislastid = 0;//0;
    for(int ires = 0; ires < commonIndexesRef.size(); ++ires)
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
    for(int ires = mnmislastid; ires < missingNonrefInd.size(); ++ires)
    {
        resultsInds.push_back(std::pair<int,int>(missingNonrefInd[ires], -1));
    }
    return resultsInds;
    */

}

//functiont that combines two alignments into the final experiment
experiment metaboCombineR::generateAlignmentExperiment(int window_size, experiment *expNonRef, experiment *exp2Ref)
{
    std::vector<std::vector<std::string > > mysamples_windows;
    for(int istart = window_size-1; istart < expNonRef->mz.size(); ++istart)
    {
        std::vector<std::string > actSample;
        actSample.resize(window_size);
        for(int iwin = 0; iwin < window_size; ++iwin)
        {
            actSample[(window_size-1)-iwin] = expNonRef->mzStr[istart-iwin];
        }
        mysamples_windows.push_back(actSample);
    }


    std::vector<alingmentIndexes> allAlignments;
    for(int isample = 0; isample < mysamples_windows.size(); ++isample)
    {
        std::list<std::string*> rowsrc_profile_mz;
        for(unsigned long irow = 0; irow < mysamples_windows[isample].size(); ++irow)
        {
            rowsrc_profile_mz.push_back(&mysamples_windows[isample][irow]);
        }
        mySolution rowAlign = computeSimilarityMatrix4CompleteLocal(rowsrc_profile_mz, &(exp2Ref->mzStr));
        alingmentIndexes resAlign = processLocalAlignments(&rowAlign, &(exp2Ref->mzStr));
        resAlign.expNonref = expNonRef;
        resAlign.expRef = exp2Ref;
        resAlign.startRow = isample;
        allAlignments.push_back(resAlign);
    }
    std::vector<std::pair<int,int> > index2table = aggregateAlignments(&allAlignments); //nonref, ref
    Rcpp::NumericMatrix mat(index2table.size(), exp2Ref->experiments.ncol() + expNonRef->experiments.ncol());
    Rcpp::CharacterVector chmat(index2table.size());
    Rcpp::CharacterVector chcolnames(exp2Ref->experiments.ncol() + expNonRef->experiments.ncol());
    //add rownames
    std::vector<std::string> newMZstr(index2table.size());
    Rcpp::NumericVector newMZ(index2table.size(), 0);
    for(int irow = 0; irow < index2table.size(); ++irow)
    {
        if(index2table[irow].first < index2table[irow].second)
        {
            chmat(irow) = std::string("M") + exp2Ref->mzStr[index2table[irow].second] + std::string("T");
            newMZstr[irow] = exp2Ref->mzStr[index2table[irow].second];
            newMZ(irow) = exp2Ref->mz(index2table[irow].second);
        }
        else
        {
            chmat(irow) = std::string("M") + expNonRef->mzStr[index2table[irow].first] + std::string("T");
            newMZstr[irow] = expNonRef->mzStr[index2table[irow].first];
            newMZ(irow) = expNonRef->mz(index2table[irow].first);
        }
    }

    Rcpp::rownames(mat) = chmat;
    //add first nonreftable
    std::vector<std::string> newRTstr(index2table.size());
    Rcpp::NumericVector newRT(index2table.size(), 0);
    for(int irow = 0; irow < index2table.size(); ++irow)
    {
        for(int icol = 0; icol < expNonRef->experiments.ncol(); ++icol)
        {
            if(index2table[irow].first == -1)
                mat(irow, icol) = NA_REAL;
            else
                mat(irow, icol) = expNonRef->experiments(index2table[irow].first, icol);
            chcolnames(icol) = expNonRef->colnames(icol);
        }
        if(index2table[irow].first != -1)
            chmat(irow) += std::string("/") + expNonRef->rtStr[index2table[irow].first];

        if(index2table[irow].first != -1 && index2table[irow].second != -1)
        {
            newRT(irow) = (expNonRef->rt(index2table[irow].first) + exp2Ref->rt(index2table[irow].second))/2.0;
        }
        else if(index2table[irow].first != -1)
        {
            newRT(irow) = expNonRef->rt(index2table[irow].first);
        }
        else
        {
            newRT(irow) = exp2Ref->rt(index2table[irow].second);
        }
        std::ostringstream rt2str;
        rt2str << newRT(irow);
        newRTstr[irow] = rt2str.str();
    }
    //add secondly reftable
    for(int irow = 0; irow < index2table.size(); ++irow)
    {
        for(int icol = expNonRef->experiments.ncol(); icol < mat.ncol(); ++icol)
        {
            if(index2table[irow].second == -1)
                mat(irow, icol) = NA_REAL;
            else
                mat(irow, icol) = exp2Ref->experiments(index2table[irow].second, icol-expNonRef->experiments.ncol());
            chcolnames(icol) = exp2Ref->colnames(icol-expNonRef->experiments.ncol());
        }
        if(index2table[irow].second != -1)
            chmat(irow) +=  std::string("/") + exp2Ref->rtStr[index2table[irow].second];
    }

    Rcpp::colnames(mat) = chcolnames;
    experiment toReturn;
    toReturn.experiments = mat;
    toReturn.mz = newMZ;
    toReturn.mzStr = newMZstr;
    toReturn.rownames = chmat;
    toReturn.rt = newRT;
    toReturn.rtStr = newRTstr;
    toReturn.colnames = chcolnames;
    return toReturn;
}

//functiont that combines two alignments into the final experiment -RT
experiment metaboCombineR::generateAlignmentExperimentRT(int window_size, experiment *expNonRef, experiment *exp2Ref)
{
    //std::vector<std::vector<std::string > > mysamples_windows;
    //for(int istart = window_size-1; istart < expNonRef->mz.size(); ++istart)
    //{
    //    std::vector<std::string > actSample;
    //    actSample.resize(window_size);
    //    for(int iwin = 0; iwin < window_size; ++iwin)
    //    {
    //       actSample[(window_size-1)-iwin] = expNonRef->mzStr[istart-iwin];
    //    }
    //    mysamples_windows.push_back(actSample);
    //}


    //std::vector<alingmentIndexes> allAlignments;
    //for(int isample = 0; isample < mysamples_windows.size(); ++isample)
    //{
        //std::list<std::string*> rowsrc_profile_mz;
        //for(unsigned long irow = 0; irow < mysamples_windows[isample].size(); ++irow)
        //{
        //    rowsrc_profile_mz.push_back(&mysamples_windows[isample][irow]);
        //}
        std::list<std::string *> nonrefmzStr;
 //       Rcpp::Rcout << "FIRST:" << std::endl;
        for(int il = 0; il < expNonRef->mzStr.size(); ++il)
        {
            nonrefmzStr.push_back(&(expNonRef->mzStr[il]));
 //           Rcpp::Rcout << *(nonrefmzStr.back()) << " ";
        }
        Rcpp::Rcout << std::endl;
//         Rcpp::Rcout << "SECOND:" << std::endl;
        for(int il = 0; il < exp2Ref->mzStr.size(); ++il)
        {
            //nonrefmzStr.push_back(&(expNonRef->mzStr[il]));
//            Rcpp::Rcout << exp2Ref->mzStr[il] << " ";
        }
// Rcpp::Rcout << std::endl;

        Rcpp::Rcout << "1#" << std::endl;
        mySolution rowAlign = computeSimilarityMatrix4Complete((nonrefmzStr), &(exp2Ref->mzStr));//computeSimilarityMatrix4CompleteSemiglobal((nonrefmzStr), &(exp2Ref->mzStr));
        //mySolution rowAlign = computeSimilarityMatrix4CompleteSemiglobal((nonrefmzStr), &(exp2Ref->mzStr));
        std::string file= "myfile.log";
        printSemiAlignmentToFile(&rowAlign, &(exp2Ref->mzStr), file);
        Rcpp::Rcout << "2#" << std::endl;
        Rcpp::Rcout << "scorerow: " << rowAlign.score << std::endl;

        //NAHRADIT SEMIGLOBAL rekonstrukci!!!!!!!!
        alingmentIndexes resAlign = processSemiAlignment(&rowAlign, &(exp2Ref->mzStr));



        //      alingmentIndexes resAlign;
        resAlign.expNonref = expNonRef;
        resAlign.expRef = exp2Ref;
        resAlign.startRow = 0;
        //allAlignments.push_back(resAlign);
    //}
    //std::vector<std::pair<int,int> > index2table = aggregateAlignments(&allAlignments); //nonref, ref

 Rcpp::Rcout << "3#" << std::endl;
    std::vector<std::pair<int,int> > index2table = rtcorrection(&resAlign, window_size); //nonref, ref
     Rcpp::Rcout << "4#" << std::endl;

    Rcpp::NumericMatrix mat(index2table.size(), exp2Ref->experiments.ncol() + expNonRef->experiments.ncol());
    Rcpp::CharacterVector chmat(index2table.size());
    Rcpp::CharacterVector chcolnames(exp2Ref->experiments.ncol() + expNonRef->experiments.ncol());
    //add rownames
    std::vector<std::string> newMZstr(index2table.size());
    Rcpp::NumericVector newMZ(index2table.size(), 0);


    for(int irow = 0; irow < index2table.size(); ++irow)
    {
        if(index2table[irow].first == -1 && index2table[irow].second == -1)
            continue;
        if(index2table[irow].first < index2table[irow].second)
        {
            chmat(irow) = std::string("M") + exp2Ref->mzStr[index2table[irow].second] + std::string("T");
            newMZstr[irow] = exp2Ref->mzStr[index2table[irow].second];
            newMZ(irow) = exp2Ref->mz(index2table[irow].second);
        }
        else
        {
            chmat(irow) = std::string("M") + expNonRef->mzStr[index2table[irow].first] + std::string("T");
            newMZstr[irow] = expNonRef->mzStr[index2table[irow].first];
            newMZ(irow) = expNonRef->mz(index2table[irow].first);
        }
    }

    Rcpp::rownames(mat) = chmat;
    //add first nonreftable
    std::vector<std::string> newRTstr(index2table.size());
    Rcpp::NumericVector newRT(index2table.size(), 0);


    for(int irow = 0; irow < index2table.size(); ++irow)
    {
        if(index2table[irow].first == -1 && index2table[irow].second == -1)
            continue;
        for(int icol = 0; icol < expNonRef->experiments.ncol(); ++icol)
        {
            if(index2table[irow].first == -1)
                mat(irow, icol) = NA_REAL;
            else
                mat(irow, icol) = expNonRef->experiments(index2table[irow].first, icol);
            chcolnames(icol) = expNonRef->colnames(icol);
        }
        if(index2table[irow].first != -1)
            chmat(irow) += std::string("/") + expNonRef->rtStr[index2table[irow].first];

        if(index2table[irow].first != -1 && index2table[irow].second != -1)
        {
            newRT(irow) = (expNonRef->rt(index2table[irow].first) + exp2Ref->rt(index2table[irow].second))/2.0;
        }
        else if(index2table[irow].first != -1)
        {
            newRT(irow) = expNonRef->rt(index2table[irow].first);
        }
        else
        {
            newRT(irow) = exp2Ref->rt(index2table[irow].second);
        }
        std::ostringstream rt2str;
        rt2str << newRT(irow);
        newRTstr[irow] = rt2str.str();
    }
    //add secondly reftable
    for(int irow = 0; irow < index2table.size(); ++irow)
    {
        if(index2table[irow].first == -1 && index2table[irow].second == -1)
            continue;
        for(int icol = expNonRef->experiments.ncol(); icol < mat.ncol(); ++icol)
        {
            if(index2table[irow].second == -1)
                mat(irow, icol) = NA_REAL;
            else
                mat(irow, icol) = exp2Ref->experiments(index2table[irow].second, icol-expNonRef->experiments.ncol());
            chcolnames(icol) = exp2Ref->colnames(icol-expNonRef->experiments.ncol());
        }
        if(index2table[irow].second != -1)
            chmat(irow) +=  std::string("/") + exp2Ref->rtStr[index2table[irow].second];
    }


    Rcpp::colnames(mat) = chcolnames;
    experiment toReturn;
    toReturn.experiments = mat;
    toReturn.mz = newMZ;
    toReturn.mzStr = newMZstr;
    toReturn.rownames = chmat;
    toReturn.rt = newRT;
    toReturn.rtStr = newRTstr;
    toReturn.colnames = chcolnames;


    //experiment toReturn;
    return toReturn;
}

//compute similarity function
double metaboCombineR::getSimilarity(experiment *first, experiment *second)
{
    std::map<std::string, int> hist1;
    std::map<std::string, int> hist2;
    std::vector<std::string> myintersection;
    int intersectionCount = 0;
    //histogram1
    for(int ife = 0; ife < first->mzStr.size(); ++ife)
    {
        hist1[first->mzStr[ife]]++;
    }
    //histogram2
    for(int ise = 0; ise < second->mzStr.size(); ++ise)
    {
        hist2[second->mzStr[ise]]++;
        if(hist1.find(second->mzStr[ise]) != hist1.end())
        {
            intersectionCount++;
            myintersection.push_back(second->mzStr[ise]);
        }
    }


    double score = first->mzStr.size() + second->mzStr.size() - (2*intersectionCount);
    for(int icommon = 0; icommon < myintersection.size(); ++icommon)
    {
        score += std::abs(hist1[myintersection[icommon]] - hist2[myintersection[icommon]]);
    }
    return score;
}

//compute the similarity for each experiment in vector and returns positions of the best
std::pair<int, int> metaboCombineR::findTheMostSimilar(std::vector<experiment> *multiAlignmentExperiment)
{
    double bestscore = DBL_MAX;
    int indNonREF;
    int indREF;
    for(int irow = 0; irow < multiAlignmentExperiment->size(); ++irow)
    {
        for(int icol = irow+1; icol < multiAlignmentExperiment->size(); ++icol)
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

//main function, this  combines all experiments
Rcpp::NumericMatrix metaboCombineR::run(Rcpp::List allExps, int mzprecision, int windowsize)
{
    //iterate over all experiments and prepare data structures
    for(int iexp = 0; iexp < allExps.size(); ++iexp)
    {
        experiment newExp;
        // Make function callable from C++
        Rcpp::Function getRTsR("getRTs");
        Rcpp::Function getMZsR("getMZs");
        Rcpp::Function sortDataFrameByRt("sortDataFrameByRt");
        Rcpp::NumericMatrix originalDataFrame = sortDataFrameByRt(allExps[iexp], mzprecision);
        Rcpp::NumericMatrix originalFrame = originalDataFrame;

        newExp.rt = getRTsR(originalDataFrame);
        newExp.mz = getMZsR(originalDataFrame, mzprecision);

        newExp.experimentID.push_back(iexp);
        newExp.mzStr.resize(newExp.rt.size());
        newExp.rtStr.resize(newExp.rt.size());
        newExp.profile_indexes.resize(newExp.rt.size());
        for(int irt = 0; irt < newExp.rt.size(); ++irt)
        {
            newExp.profile_indexes[irt] = irt;
            std::ostringstream mzstrstream;
            std::ostringstream rtstrstream;
            mzstrstream << newExp.mz(irt);
            rtstrstream << newExp.rt(irt);
            newExp.mzStr[irt] = mzstrstream.str();
            newExp.rtStr[irt] = rtstrstream.str();

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
        std::pair<int, int> mostsimilar = findTheMostSimilar(&multiAlignmentExperiment);
        Rcpp::Rcout << ".";        
        std::vector<experiment> newExp;
        for(int iel = 0; iel < multiAlignmentExperiment.size(); ++iel)
        {
            if(iel != mostsimilar.first && iel != mostsimilar.second)
            {
                newExp.push_back(multiAlignmentExperiment[iel]);
            }
        }
        newExp.push_back(generateAlignmentExperiment(window_size, &(multiAlignmentExperiment[mostsimilar.first]), &(multiAlignmentExperiment[mostsimilar.second])));
        multiAlignmentExperiment = newExp;        
    } while(multiAlignmentExperiment.size() > 1);

    Rcpp::Rcout << " Done" << std::endl;
    return multiAlignmentExperiment.back().experiments;
}

//main function, this  combines all experiments
Rcpp::NumericMatrix metaboCombineR::runRT(Rcpp::List allExps, int mzprecision, int rtwindowsize)
{
    //iterate over all experiments and prepare data structures
    for(int iexp = 0; iexp < allExps.size(); ++iexp)
    {
        experiment newExp;
        // Make function callable from C++
        Rcpp::Function getRTsR("getRTs");
        Rcpp::Function getMZsR("getMZs");
        Rcpp::Function sortDataFrameByRt("sortDataFrameByRt");
        Rcpp::NumericMatrix originalDataFrame = sortDataFrameByRt(allExps[iexp], mzprecision);
        Rcpp::NumericMatrix originalFrame = originalDataFrame;

        newExp.rt = getRTsR(originalDataFrame);
        newExp.mz = getMZsR(originalDataFrame, mzprecision);

        newExp.experimentID.push_back(iexp);
        newExp.mzStr.resize(newExp.rt.size());
        newExp.rtStr.resize(newExp.rt.size());
        newExp.profile_indexes.resize(newExp.rt.size());
        for(int irt = 0; irt < newExp.rt.size(); ++irt)
        {
            newExp.profile_indexes[irt] = irt;
            std::ostringstream mzstrstream;
            std::ostringstream rtstrstream;
            mzstrstream << newExp.mz(irt);
            rtstrstream << newExp.rt(irt);
            newExp.mzStr[irt] = mzstrstream.str();
            newExp.rtStr[irt] = rtstrstream.str();

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
        std::pair<int, int> mostsimilar = findTheMostSimilar(&multiAlignmentExperiment);
        Rcpp::Rcout << ".";
        std::vector<experiment> newExp;
        for(int iel = 0; iel < multiAlignmentExperiment.size(); ++iel)
        {
            if(iel != mostsimilar.first && iel != mostsimilar.second)
            {
                newExp.push_back(multiAlignmentExperiment[iel]);
            }
        }
        newExp.push_back(generateAlignmentExperimentRT(window_size, &(multiAlignmentExperiment[mostsimilar.first]), &(multiAlignmentExperiment[mostsimilar.second])));
        multiAlignmentExperiment = newExp;
    } while(multiAlignmentExperiment.size() > 1);

    Rcpp::Rcout << " Done" << std::endl;
    return multiAlignmentExperiment.back().experiments;
}

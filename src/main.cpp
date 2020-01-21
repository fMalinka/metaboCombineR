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
    metaboCombineR() : msg("hello") {}
    Rcpp::NumericMatrix run(Rcpp::List allExps, int mzprecision, int windowsize);
    void combineExperiments();
    Rcpp::NumericMatrix getfinalMatrix();

private:
    std::string msg;
    std::vector<experiment> experiments;
    Rcpp::NumericMatrix finalMatrix;
};


RCPP_MODULE(metaboCombineR){
    using namespace Rcpp;

    class_<metaboCombineR>("metaboCombineR")
    // expose the default constructor
    .constructor()

    .method("run", &metaboCombineR::run, "get the message")
    ;
}

std::vector<string> loadCSVDataMZ(string path)
{
    std::ifstream myfile(path.c_str());
    std::string line;
    std::getline(myfile, line); //head
    std::vector<string> myvector;
    //myvector.push_back("XXXX");
    while(std::getline(myfile, line))
    {
        if(SILENT)
            Rcpp::Rcout << line << std::endl;
        myvector.push_back(line);
        //myvector.push_back(line);
    }
    return myvector;
}

bool toAddmz(int index, std::vector<int> *array, int lastid)
{
    for(int i = lastid; i < array->size(); ++i)
    {
 //       if(index < (*array)[i])
 //           return false;
        if((*array)[i] == index)
            return true;
    }
    return false;
}

int checkIndexState(std::pair<int,int> *line)
{
    if(line->first == -1)
        return ONLYREF;
    else if(line->second == -1)
        return ONLYNONREF;
    else
        return MATCH;
}

std::vector<std::pair<int,int> > aggregateAlignments(std::vector<alingmentIndexes> *allAlignments)
{

    std::vector<std::pair<int,int> > resultsInds;   //nonref, ref
    std::vector<std::vector<int> > commonIndexesRefCandidates((*allAlignments)[0].expRef->mz.size());
    std::vector<int> commonIndexesRef((*allAlignments)[0].expRef->mz.size(), -1);
//Rcpp::Rcout << "_E0" << std::endl;
    for(int iali = 0; iali < allAlignments->size(); ++iali)
    {
        for(int irow = 0; irow < (*allAlignments)[iali].indexRef.size(); ++irow)
        {
            if((*allAlignments)[iali].indexRef[irow] != -1 && (*allAlignments)[iali].indexNonref[irow] != -1)
            {
                commonIndexesRef[(*allAlignments)[iali].indexRef[irow]] = (*allAlignments)[iali].indexNonref[irow]+(*allAlignments)[iali].startRow;
                commonIndexesRefCandidates[(*allAlignments)[iali].indexRef[irow]].push_back((*allAlignments)[iali].indexNonref[irow]+(*allAlignments)[iali].startRow);
                //Rcpp::Rcout << "as: " << (*allAlignments)[iali].indexNonref[irow]+(*allAlignments)[iali].startRow << std::endl;
            }

        }
    }
//Rcpp::Rcout << "_E1" << std::endl;
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
//Rcpp::Rcout << "_E2" << std::endl;
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
//Rcpp::Rcout << "_E3" << std::endl;
    std::vector<int> missingNonrefInd;
    for(int ic = 0; ic < (*allAlignments)[0].expNonref->mzStr.size(); ++ic)
    {
        if(missingNonrefOccurance.find(ic) == missingNonrefOccurance.end())
        {
            missingNonrefInd.push_back(ic);
            //Rcpp::Rcout << "ic: " << ic << std::endl;
        }
    }
//Rcpp::Rcout << "_E4" << std::endl;
    /*
    for(int imis = 0; imis < missingNonref.size(); ++imis)
    {
        if(missingNonref[imis] == -1)
        {
            missingNonrefInd.push_back(imis);

        }
    }
    */

    int nCommonHits = 0;//0;//0;
    int mnmislastid = 0;//0;
    //Rcpp::Rcout << "size ref: " << commonIndexesRef.size() << " non: " << missingNonrefInd.size() << std::endl;
    for(int ires = 0; ires < commonIndexesRef.size(); ++ires)
    {
        //Rcpp::Rcout << " " << ires;
        if(commonIndexesRef[ires] != -1)
        {            
            //only on nonref
            //Rcpp::Rcout << " minnsin: " << missingNonrefInd[mnmislastid] << " common: " << commonIndexesRef[ires] << std::endl;
            if(!missingNonrefInd.empty())
            {
                //Rcpp::Rcout << "mnmislastid: " << mnmislastid << " missingNonrefInd: " << missingNonrefInd.size() << std::endl;
                if(mnmislastid < (missingNonrefInd.size()))
                {
                    //Rcpp::Rcout << "passs" << std::endl;
                    while(missingNonrefInd[mnmislastid] < commonIndexesRef[ires])
                    {

                        resultsInds.push_back(std::pair<int,int>(missingNonrefInd[mnmislastid], -1));
                        //Rcpp::Rcout << "mmnilastid: " << mnmislastid << " nonrefind: " << missingNonrefInd[mnmislastid] << std::endl;
                        mnmislastid++;
                        if(mnmislastid >= missingNonrefInd.size())
                            break;
                    }
                }
     /*           else
                {
                    for(int ires = mnmislastid; ires < missingNonrefInd.size(); ++ires)
                    {
                        resultsInds.push_back(std::pair<int,int>(missingNonrefInd[ires], -1));
                        mnmislastid++;
                        //Rcpp::Rcout << missingNonrefInd[ires] << std::endl;
                    }
                }
                */
            }

            //common
            resultsInds.push_back(std::pair<int,int>(commonIndexesRef[ires], ires));


/*

            resultsInds.push_back(std::pair<int,int>(nCommonHits, -1));
                Rcpp::Rcout << nCommonHits << std::endl;
            //++nCommonHits;
//resultsInds.push_back(std::pair<int,int>(nCommonHits, -1));
//resultsInds.push_back(std::pair<int,int>(commonIndexesRef[nCommonHits], -1));
            //if(toAddmz(nCommonHits, &missingNonrefInd, mnmislastid))
            if(toAddmz(nCommonHits, &missingNonrefInd, 0))
            {
                resultsInds.push_back(std::pair<int,int>(nCommonHits, -1));
                ++mnmislastid;
            }
            //resultsInds.push_back(std::pair<int,int>(commonIndexesRef[ires], ires));
            ++nCommonHits;
*/



        }
        else
        {
            resultsInds.push_back(std::pair<int,int>(-1, ires));
        }        
    }
//    Rcpp::Rcout << std::endl << "_E5" << std::endl;
    //end

    //Rcpp::Rcout << "END\n" << "mnmislastid: " << mnmislastid << " missingNonrefInd: " << missingNonrefInd.size() << std::endl;
    for(int ires = mnmislastid; ires < missingNonrefInd.size(); ++ires)
    {
        //Rcpp::Rcout << "ADDEDEDED\n";
        resultsInds.push_back(std::pair<int,int>(missingNonrefInd[ires], -1));
        //Rcpp::Rcout << missingNonrefInd[ires] << std::endl;
    }


//Rcpp::Rcout << "_E6" << std::endl;
    return resultsInds;
}

Rcpp::NumericMatrix generateAlignmentMatrix(int window_size, experiment *expNonRef, experiment *exp2Ref)
{
    //test
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
        //printAlignment(&rowAlign, &mysamples[1]);
//        std::string name;
//        std::stringstream ss;
//        ss << isample;
//        name = "align_" + ss.str();//std::to_string(isample);
//        Rcpp::Rcout << name << std::endl;
//        printLocalAlignmentToFile(&rowAlign, &exp2Ref->mzStr, name);
        alingmentIndexes resAlign = processLocalAlignments(&rowAlign, &(exp2Ref->mzStr));
        resAlign.expNonref = expNonRef;
        resAlign.expRef = exp2Ref;
        resAlign.startRow = isample;
        allAlignments.push_back(resAlign);
    }

    std::vector<std::pair<int,int> > index2table = aggregateAlignments(&allAlignments); //nonref, ref

    Rcpp::NumericMatrix mat(index2table.size(), exp2Ref->experiments.ncol() + expNonRef->experiments.ncol());
    Rcpp::CharacterVector chmat(index2table.size());
    //add rownames
    for(int irow = 0; irow < index2table.size(); ++irow)
    {
        if(index2table[irow].first < index2table[irow].second)
        {
            chmat(irow) = exp2Ref->mzStr[index2table[irow].second] + std::string("T");
        }
        else
        {
            chmat(irow) = expNonRef->mzStr[index2table[irow].first] + std::string("T");
        }
    }

    Rcpp::rownames(mat) = chmat;
    //add first nonreftable
    for(int irow = 0; irow < index2table.size(); ++irow)
    {
        for(int icol = 0; icol < exp2Ref->experiments.ncol(); ++icol)
        {
            if(index2table[irow].first == -1)
                mat(irow, icol) = NA_REAL;
            else
                mat(irow, icol) = expNonRef->experiments(index2table[irow].first, icol);
        }
        if(index2table[irow].first != -1)
            chmat(irow) += std::string("/") + expNonRef->rtStr[index2table[irow].first];
    }
    //add secondly reftable
    for(int irow = 0; irow < index2table.size(); ++irow)
    {
        for(int icol = exp2Ref->experiments.ncol(); icol < mat.ncol(); ++icol)
        {
            if(index2table[irow].second == -1)
                mat(irow, icol) = NA_REAL;
            else
                mat(irow, icol) = exp2Ref->experiments(index2table[irow].second, icol-exp2Ref->experiments.ncol());
        }
        if(index2table[irow].second != -1)
            chmat(irow) +=  std::string("/") + exp2Ref->rtStr[index2table[irow].second];
    }

    return mat;
}

experiment generateAlignmentExperiment(int window_size, experiment *expNonRef, experiment *exp2Ref)
{
    //test
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
//        Rcpp::Rcout << "B computeSimilarityMatrix4CompleteLocal" << std::endl;
        mySolution rowAlign = computeSimilarityMatrix4CompleteLocal(rowsrc_profile_mz, &(exp2Ref->mzStr));
//        Rcpp::Rcout << "E computeSimilarityMatrix4CompleteLocal" << std::endl;
        //printAlignment(&rowAlign, &mysamples[1]);
        //std::string name;
        //std::stringstream ss;
        //ss << isample;
        //name = "align_" + ss.str();//std::to_string(isample);
        //Rcpp::Rcout << name << std::endl;
        //printLocalAlignmentToFile(&rowAlign, &(exp2Ref->mzStr), name);
//        Rcpp::Rcout << "B processLocalAlignments" << std::endl;
        alingmentIndexes resAlign = processLocalAlignments(&rowAlign, &(exp2Ref->mzStr));
//        Rcpp::Rcout << "E processLocalAlignments" << std::endl;
        resAlign.expNonref = expNonRef;
        resAlign.expRef = exp2Ref;
        resAlign.startRow = isample;
        allAlignments.push_back(resAlign);
 //       Rcpp::Rcout << "isample: " << isample << "/" << mysamples_windows.size() << std::endl;
    }
//Rcpp::Rcout << "_E_" << std::endl;
    std::vector<std::pair<int,int> > index2table = aggregateAlignments(&allAlignments); //nonref, ref
//Rcpp::Rcout << "_F_" << std::endl;
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
            chmat(irow) = exp2Ref->mzStr[index2table[irow].second] + std::string("T");
            newMZstr[irow] = exp2Ref->mzStr[index2table[irow].second];
            newMZ(irow) = exp2Ref->mz(index2table[irow].second);
        }
        else
        {
            chmat(irow) = expNonRef->mzStr[index2table[irow].first] + std::string("T");
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

std::vector<std::pair<int,int> > generateAlignmentIndexes(int window_size, experiment *expNonRef, experiment *exp2Ref)
{
    //test
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
    return index2table;
}

double getSimilarity(experiment *first, experiment *second)
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

std::pair<int, int> findTheMostSimilar(std::vector<experiment> *multiAlignmentExperiment)
{
    double bestscore = 1000000;
    int indNonREF;
    int indREF;
    for(int irow = 0; irow < multiAlignmentExperiment->size(); ++irow)
    {
        for(int icol = irow+1; icol < multiAlignmentExperiment->size(); ++icol)
        {
            double score = getSimilarity(&((*multiAlignmentExperiment)[irow]), &((*multiAlignmentExperiment)[icol]));
            //Rcpp::Rcout << "irow: " << irow << " icol: " << icol <<" score: " << score << std::endl;
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


Rcpp::NumericMatrix metaboCombineR::run(Rcpp::List allExps, int mzprecision, int windowsize)
{
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
        //init row index
        newExp.profile_indexes.resize(newExp.rt.size());
        for(int irt = 0; irt < newExp.rt.size(); ++irt)
        {
            newExp.profile_indexes[irt] = irt;
            std::ostringstream mzstrstream;
            std::ostringstream rtstrstream;
            mzstrstream << newExp.mz(irt);
            rtstrstream << newExp.rt(irt);
            newExp.mzStr[irt] = mzstrstream.str();
            //Rcpp::Rcout << mzstrstream.str() << " a " << newExp.mz(irt) << std::endl;
            newExp.rtStr[irt] = rtstrstream.str();

        }        
        newExp.experiments = originalFrame;//(Rcpp::_, originalFrame.ncol());
        Rcpp::CharacterVector dataColnames = Rcpp::colnames(originalDataFrame);
        newExp.colnames = dataColnames;
        for(int icolname = 0; icolname < dataColnames.size(); ++icolname)
        {
            std::ostringstream expname;
            expname << "subExp" << iexp << "_" << dataColnames(icolname);
            newExp.colnames(icolname) = expname.str();
        }

        //newExp.experimentName = expname.str();
        this->experiments.push_back(newExp);
    }
    //Rcpp::Rcout << "------\n";

    std::vector<experiment> multiAlignmentExperiment;
    Rcpp::Rcout << "Successfully loaded " << this->experiments.size()  << " datasets!" << std::endl;


//    Rcpp::Rcout << "Finding centroid with the highest number of peaks..." << std::endl;
//    int centroid = 0;
//    int maxcentroid = 0;
    for(int iexp = 0; iexp < allExps.size(); ++iexp)
    {
//        if(maxcentroid < this->experiments[iexp].mz.size())
//        {
//            maxcentroid = this->experiments[iexp].mz.size();
//            centroid = iexp;
//        }
        multiAlignmentExperiment.push_back(this->experiments[iexp]);
    }

//    Rcpp::Rcout << "The centroid is data table number " << centroid << " with " << maxcentroid << " peaks" << std::endl;

    //experiment exp2Ref = this->experiments[centroid];
    //std::vector<Rcpp::NumericMatrix> multiAlignments(this->experiments.size());
    //std::vector<std::vector<std::pair<int,int> > > multiAlignmentIndexes(this->experiments.size());
    int window_size = windowsize;

    //std::vector<std::vector<double> > simmatrix(multiAlignmentExperiment.size());

    //multiAlignmentExperiment.push_back(this->experiments[0]);
    //multiAlignmentExperiment.push_back(this->experiments[1]);
    //multiAlignmentExperiment.push_back(this->experiments[2]);

    Rcpp::Rcout << "Working .";
    do
    {
        std::pair<int, int> mostsimilar = findTheMostSimilar(&multiAlignmentExperiment);
        Rcpp::Rcout << ".";
        //Rcpp::Rcout << "non ref: " << mostsimilar.first  << " (" << multiAlignmentExperiment[mostsimilar.first].mzStr.size() << ", row: " << multiAlignmentExperiment[mostsimilar.first].experiments.nrow()  << ", col: " << multiAlignmentExperiment[mostsimilar.first].experiments.ncol() << ") ref: " << mostsimilar.second << " (" << multiAlignmentExperiment[mostsimilar.second].mzStr.size() << ", row: " << multiAlignmentExperiment[mostsimilar.second].experiments.nrow()  << ", col: " << multiAlignmentExperiment[mostsimilar.second].experiments.ncol()  << ")" << std::endl;
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
        //Rcpp::Rcout << "Making an alignment from subexperiment " <<
        //Rcpp::Rcout << "new size: " << multiAlignmentExperiment.size() << std::endl;
    } while(multiAlignmentExperiment.size() > 1);
    Rcpp::Rcout << " Done" << std::endl;
    return multiAlignmentExperiment.back().experiments;

    /*
    double bestscore = 1000000;
    int indNonREF;
    int indREF;
    for(int irow = 0; irow < multiAlignmentExperiment.size(); ++irow)
    {
        for(int icol = irow+1; icol < multiAlignmentExperiment.size(); ++icol)
        {
            double score = getSimilarity(&(multiAlignmentExperiment[irow]), &(multiAlignmentExperiment[icol]));
            Rcpp::Rcout << "irow: " << irow << " icol: " << icol <<" score: " << score << std::endl;
            if(score < bestscore)
            {
                bestscore = score;
                if(multiAlignmentExperiment[irow].mzStr.size() < multiAlignmentExperiment[icol].mzStr.size())
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
    */
 //   while(!multiAlignmentExperiment.empty())
 //   {
 //       ;
 //   }


    /*
    for(int iali = 0; iali < this->experiments.size(); ++iali)
    {
        if(iali != centroid)
        {
            multiAlignmentExperiment.push_back(this->experiments[iali]);
        }
    }

    experiment aaaa = generateAlignmentExperiment(window_size, &(this->experiments[0]), &exp2Ref);
    experiment bbbb = generateAlignmentExperiment(window_size, &(this->experiments[1]), &aaaa);
    return bbbb.experiments;
*/
//    while(!multiAlignmentExperiment.empty())
//    {
//        Rcpp::Rcout << "size: " << multiAlignmentExperiment.size() << std::endl;
//        experiment toadd = multiAlignmentExperiment.front();
//        multiAlignmentExperiment.pop_front();
//        exp2Ref = generateAlignmentExperiment(window_size, &toadd, &exp2Ref);
//    }


    //generateAlignmentMatrix(window_size, &expNonRef, &exp2Ref);

    //#pragma omp parallel for
    //for(int iali = 0; iali < this->experiments.size(); ++iali)
    //{
    //    if(iali != centroid)
    //    {
    //        experiment expNonRef = this->experiments[0];

//            experiment newexp = generateAlignmentExperiment(window_size, &expNonRef, &exp2Ref);
//            multiAlignmentExperiment[0] =  newexp;
//            return exp2Ref.experiments;
    //    }
    //}


    /*
    #pragma omp parallel for
    for(int iali = 0; iali < this->experiments.size(); ++iali)
    {
        if(iali != centroid)
        {
            experiment expNonRef = this->experiments[iali];
            std::vector<std::pair<int,int> > results = generateAlignmentIndexes(window_size, &expNonRef, &exp2Ref);
            multiAlignmentIndexes[iali] =  results;
        }
    }

    for(int ia = 0; ia < multiAlignmentIndexes[0].size(); ++ia)
    {
        Rcpp::Rcout << "first: " << multiAlignmentIndexes[0][ia].first << " second: " << multiAlignmentIndexes[0][ia].second << std::endl;
    }
    Rcpp::Rcout << "---------------------------------------------" << std::endl;
    for(int ia = 0; ia < multiAlignmentIndexes[1].size(); ++ia)
    {
        Rcpp::Rcout << "first: " << multiAlignmentIndexes[1][ia].first << " second: " << multiAlignmentIndexes[1][ia].second << std::endl;
    }

    std::vector<std::vector<int> > alignIndexREF(multiAlignmentIndexes.size()-1);   //align/row
    std::vector<std::vector<std::vector<int> > > alignIndexNONREF(multiAlignmentIndexes.size()-1);   //align/ref//row
    for(int icombine = 0; icombine < multiAlignmentIndexes.size()-1; ++icombine)
    {
        if(!multiAlignmentIndexes[icombine].empty())
        {
            //looking for refids
            std::vector<int> nonrefacts;
            for(int irow = 0; irow < multiAlignmentIndexes[icombine].size(); ++irow)
            {
                int state = checkIndexState(&(multiAlignmentIndexes[icombine][irow]));
                if(state == MATCH || state == ONLYREF)
                {
                    alignIndexREF[icombine].push_back(irow);
                    //if(!nonrefacts.empty())
                    //{
                        alignIndexNONREF[icombine].push_back(nonrefacts);
                        Rcpp::Rcout  << "ref..." << std::endl;

                        for(int itest = 0; itest < alignIndexNONREF[icombine].back().size(); ++itest)
                        {
                            experiment expNonRef = this->experiments[icombine];
                            Rcpp::Rcout << "NONREF: " << expNonRef.mzStr[multiAlignmentIndexes[icombine][alignIndexNONREF[icombine].back()[itest]].first] << " ";
                        }
                        Rcpp::Rcout << std::endl;
                        Rcpp::Rcout << "REF! " << exp2Ref.mzStr[multiAlignmentIndexes[icombine][irow].second] << std::endl;
                        nonrefacts.clear();
                    //}
                    //alignIndexNONREF[icombine].push_back(nonrefacts);
                        //experiment expNonRef = this->experiments[icombine];
                        //Rcpp::Rcout << "REF! " << multiAlignmentIndexes[icombine][irow].second << " - "  << expNonRef.mzStr[multiAlignmentIndexes[icombine][irow].second] << std::endl;
                }
                else
                {
                    Rcpp::Rcout  << "nonref..." << std::endl;
                    //experiment expNonRef = this->experiments[icombine];
                    //Rcpp::Rcout << "NONREF! " << expNonRef.mzStr[multiAlignmentIndexes[icombine][irow].first] << std::endl;
                    nonrefacts.push_back(irow);
                }
            }
        }
    }


    for(int irow = 0; irow < alignIndexREF[0].size(); ++irow)
    {
        for(int icombine = 0; icombine < multiAlignmentIndexes.size()-1; ++icombine)
        {
            if(!multiAlignmentIndexes[icombine].empty())
            {
                Rcpp::Rcout << " |# " << alignIndexNONREF[icombine][irow].size();
                experiment expNonRef = this->experiments[icombine];
                for(int ii = 0; ii < alignIndexNONREF[icombine][irow].size(); ++ii)
                {
                    Rcpp::Rcout << "!";
                    //Rcpp::Rcout << expNonRef.mzStr[multiAlignmentIndexes[icombine][alignIndexNONREF[icombine][irow][ii]].second]  << " ";
                    //Rcpp::Rcout << alignIndexNONREF[icombine][irow][ii] << " ";
                    Rcpp::Rcout << expNonRef.mzStr[multiAlignmentIndexes[icombine][alignIndexNONREF[icombine][irow][ii]].first] << " ";
                    //Rcpp::Rcout << expNonRef.mzStr[alignIndexNONREF[icombine][irow][ii]] << " ";
                }

                Rcpp::Rcout << " #} ";
                Rcpp::Rcout << exp2Ref.mzStr[multiAlignmentIndexes[icombine][alignIndexREF[icombine][irow]].second]  << " ";
            }
        }
        Rcpp::Rcout << " (" << exp2Ref.mzStr[irow]  << ") ";
        Rcpp::Rcout << std::endl;
    }

*/


    /*

    //combine all pairwise alignments into one table
    //1. look at each alignment
    std::vector<int> lastmatchIndex(multiAlignmentIndexes.size(), -1);  //multiAlignmentIndexes - nonref,ref
    std::vector<int> lastIndex(multiAlignmentIndexes.size(), 0);

    std::vector<std::vector<int> > refIds(multiAlignmentIndexes.size());
    std::vector<std::vector<std::vector<int> > > norefCandidatasIds_interval(multiAlignmentIndexes.size()); //experiment/number/indexes
    //looking for refids
    for(int ialign = 0; ialign < multiAlignmentIndexes.size(); ++ialign)
    {
        //checkIndexState(&((multiAlignmentIndexes[ialign])[lastIndex[ialign]]));
        if(ialign != centroid)
        {
            //Rcpp::Rcout << "size: " << multiAlignmentIndexes[ialign].size() << std::endl;
            int statmatch = 0; int statonlyref = 0; int statonlynonref = 0; int statmatchref = 0;
            std::vector<int> nonrefInterval;
            for(int iind = 0; iind < multiAlignmentIndexes[ialign].size(); ++iind)
            {
                int state = checkIndexState(&(multiAlignmentIndexes[ialign][iind]));
                if(state == MATCH)
                    statmatch++;
                else if(state == ONLYREF)
                    statonlyref++;
                else
                    statonlynonref++;
                if(state == MATCH || state == ONLYREF)
                {
                    //exp2Ref.mzStr[multiAlignmentIndexes[ialign][iind].first]
                    Rcpp::Rcout << exp2Ref.mzStr[multiAlignmentIndexes[ialign][iind].second]  << " ";

                    statmatchref++;
                    refIds[ialign].push_back(iind);
                    //Rcpp::Rcout << " " << iind;

                    //close interval
                    if(!nonrefInterval.empty())
                    {
                        if(iind-1 != nonrefInterval.back())
                            nonrefInterval.push_back(iind-1);

                        norefCandidatasIds_interval[ialign].push_back(nonrefInterval);
                        nonrefInterval.clear();
                    }
                }
                else
                {
                    nonrefInterval.push_back(iind);
                    //Rcpp::Rcout << this->experiments[ialign].mzStr[multiAlignmentIndexes[ialign][iind].first]  << "(" << multiAlignmentIndexes[ialign][iind].first << ")" << multiAlignmentIndexes[ialign][iind].first <<  " ";
                }
            }
            //Rcpp::Rcout << std::endl <<"ialign: " << ialign << " match:" << statmatch << " onlyref:" << statonlyref << " onlynonref:" << statonlynonref << " statmatchRef: " << statmatchref << std::endl;
        }
        //lastIndex[ialign]
        //
        for(int iii=0; iii < norefCandidatasIds_interval[ialign].size(); ++iii)
        {
            for(int iiii = 0; iiii <norefCandidatasIds_interval[ialign][iii].size(); ++iiii)
            {
                //Rcpp::Rcout << norefCandidatasIds_interval[ialign][iii][iiii] << " ";
            }
            //Rcpp::Rcout << "***" << std::endl;
            int sec = multiAlignmentIndexes[ialign][iii].second;
            int frst = multiAlignmentIndexes[ialign][iii].first;
            if(sec == -1)
                ;//Rcpp::Rcout << " f: " << this->experiments[ialign].mzStr[frst] << " s: -1" ;
            else if(frst == -1)
                ;//Rcpp::Rcout << " f: -1  s: " << exp2Ref.mzStr[sec] ;
            else
                ;//Rcpp::Rcout << " f: " << this->experiments[ialign].mzStr[multiAlignmentIndexes[ialign][iii].first] << " s: " << exp2Ref.mzStr[multiAlignmentIndexes[ialign][iii].second] ;
            //Rcpp::Rcout << std::endl;

        }
        Rcpp::Rcout << "--------------" << std::endl;
    }
    */

//    Rcpp::NumericMatrix results;
//    Rcpp::Rcout << "END" << std::endl;
//    return results;
}


double computeScore(std::vector<std::string> sample1, std::vector<std::string> sample2, unsigned long FILTER_LIMIT)
{
    //std::vector<std::string> sample2 = loadCSVDataMZ(tgtpath);  //tgt
    //std::vector<std::string> sample1 = loadCSVDataMZ(srcpath);  //src

    //prepare data
    if(SILENT)
        Rcpp::Rcout << "generating initial set..";
    std::map<std::string, bool> commonSeqs;
    std::vector<std::string> src_profile = sample1;   //obsahuje puvodni src retezec
    std::vector<std::string> tgt_profile = sample2;   //obsahuje puvodni src retezec
    std::multimap<std::string, unsigned long> peakPosTgz; //obsahuje pozici peaku u target - pocitano pro vzdalenost
    double normalizator = std::max(src_profile.size(), tgt_profile.size());
    //multipleAlingmentR tgt;
    mySolutionShort init = generateInitionSet(&sample1, &sample2, &commonSeqs, &peakPosTgz);
    if(SILENT)
        Rcpp::Rcout << "OK" << std::endl;
    //queue alignments
    std::priority_queue<mySolutionShort> newTmp;
    newTmp.push(init);

    clock_t begin = clock();
    //iterate over rows/src sequences
    for(unsigned long irow = 0; irow < src_profile.size(); ++irow)
    {
        //newly generated candidates
        std::priority_queue<mySolutionShort> newTmp2;
       //325
        if(irow == 32)
        {
            Rcpp::Rcout << "ttttt";
        }

        //for each myALIGN in OPEN
        unsigned long opentmpSIZE = newTmp.size();
        for(unsigned long ialign = 0; ialign < opentmpSIZE; ++ialign)
        {
            mySolutionShort actSol = newTmp.top();
            newTmp.pop();
            //look at the row irow
            //make new combinations
            std::vector<mySolutionShort> mycands = makeCombinations(irow, &actSol, &commonSeqs, &src_profile, &peakPosTgz);
            //evaluate candidates
            for(unsigned long isol = 0; isol < mycands.size(); ++isol)
            {
                mySolutionShort newscore = computeSimilarityMatrix4(mycands[isol], irow, &tgt_profile, &normalizator);
                newTmp2.push(newscore);
            }
            //filter candidates
            while(newTmp2.size() > FILTER_LIMIT)
            {
                mySolutionShort todelete = newTmp2.top();
                freeMatrix(&todelete);
                newTmp2.pop();
            }
            freeMatrix(&actSol);
        }
        //continue
        newTmp = newTmp2;
        while(newTmp2.empty())
        {
            mySolutionShort todelete = newTmp2.top();
            freeMatrix(&todelete);
            newTmp2.pop();
        }
       // mySolution best = newTmp2.top();
        while(newTmp.size() > FILTER_LIMIT)
        {
            mySolutionShort todelete = newTmp.top();
            freeMatrix(&todelete);
            newTmp.pop();
        }
        if(SILENT)
            Rcpp::Rcout << "row " << irow+1 << " / " << src_profile.size() << std::endl;
    }
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    //find the best
    while(newTmp.size() != 1)
    {
        mySolutionShort todelete = newTmp.top();
        freeMatrix(&todelete);
        newTmp.pop();
    }

    mySolutionShort best = newTmp.top();
    double finalscore = best.score;
    newTmp.pop();
    freeMatrix(&best);
    return finalscore;
}

std::vector<std::string> makeAlignment(std::vector<std::string> sample1, std::vector<std::string> sample2, unsigned long FILTER_LIMIT)
{
    //std::vector<std::string> sample2 = loadCSVDataMZ(tgtpath);  //tgt
    //std::vector<std::string> sample1 = loadCSVDataMZ(srcpath);  //src

    //prepare data
    if(SILENT)
        Rcpp::Rcout << "generating initial set..";
    std::map<std::string, bool> commonSeqs;
    std::vector<std::string> src_profile = sample1;   //obsahuje puvodni src retezec
    std::vector<std::string> tgt_profile = sample2;   //obsahuje puvodni src retezec
    std::multimap<std::string, unsigned long> peakPosTgz; //obsahuje pozici peaku u target - pocitano pro vzdalenost
    double normalizator = std::max(src_profile.size(), tgt_profile.size());
    //multipleAlingmentR tgt;
    mySolutionShort init = generateInitionSet(&sample1, &sample2, &commonSeqs, &peakPosTgz);
    if(SILENT)
        Rcpp::Rcout << "OK" << std::endl;
    //queue alignments
    std::priority_queue<mySolutionShort> newTmp;
    newTmp.push(init);

    clock_t begin = clock();
    //iterate over rows/src sequences
    for(unsigned long irow = 0; irow < src_profile.size(); ++irow)
    {
        //newly generated candidates
        std::priority_queue<mySolutionShort> newTmp2;

        //for each myALIGN in OPEN
        unsigned long opentmpSIZE = newTmp.size();
        for(unsigned long ialign = 0; ialign < opentmpSIZE; ++ialign)
        {
            mySolutionShort actSol = newTmp.top();
            newTmp.pop();
            //look at the row irow
            //make new combinations
            std::vector<mySolutionShort> mycands = makeCombinations(irow, &actSol, &commonSeqs, &src_profile, &peakPosTgz);
            //evaluate candidates
            for(unsigned long isol = 0; isol < mycands.size(); ++isol)
            {
                mySolutionShort newscore = computeSimilarityMatrix4(mycands[isol], irow, &tgt_profile, &normalizator);
                newTmp2.push(newscore);
            }
            //filter candidates
            while(newTmp2.size() > FILTER_LIMIT)
            {
                mySolutionShort todelete = newTmp2.top();
                freeMatrix(&todelete);
                newTmp2.pop();
            }
            freeMatrix(&actSol);
        }
        //continue
        newTmp = newTmp2;
        while(newTmp2.empty())
        {
            mySolutionShort todelete = newTmp2.top();
            freeMatrix(&todelete);
            newTmp2.pop();
        }
       // mySolution best = newTmp2.top();
        while(newTmp.size() > FILTER_LIMIT)
        {
            mySolutionShort todelete = newTmp.top();
            freeMatrix(&todelete);
            newTmp.pop();
        }
        if(SILENT)
            Rcpp::Rcout << "row " << irow+1 << " / " << src_profile.size() << std::endl;
    }
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    //find the best
    while(newTmp.size() != 1)
    {
        mySolutionShort todelete = newTmp.top();
        freeMatrix(&todelete);
        newTmp.pop();
    }
    mySolutionShort bestnoncomplete = newTmp.top();
    mySolution best = computeSimilarityMatrix4Complete(bestnoncomplete.src.profile_mz, &tgt_profile);
    std::vector<std::string> resalign = mergeAlignments(&best, &sample2);
    freeMatrix(&bestnoncomplete);
    return resalign;
}

void make2ColAlignment(std::vector<std::string> sample1, std::vector<std::string> sample2, unsigned long FILTER_LIMIT)
{
    //std::vector<std::string> sample2 = loadCSVDataMZ(tgtpath);  //tgt
    //std::vector<std::string> sample1 = loadCSVDataMZ(srcpath);  //src

    //prepare data
    if(SILENT)
        Rcpp::Rcout << "generating initial set..";
    std::map<std::string, bool> commonSeqs;
    std::vector<std::string> src_profile = sample1;   //obsahuje puvodni src retezec
    std::vector<std::string> tgt_profile = sample2;   //obsahuje puvodni src retezec
    std::multimap<std::string, unsigned long> peakPosTgz; //obsahuje pozici peaku u target - pocitano pro vzdalenost
    double normalizator = std::max(src_profile.size(), tgt_profile.size());
    //multipleAlingmentR tgt;
    mySolutionShort init = generateInitionSet(&sample1, &sample2, &commonSeqs, &peakPosTgz);
    if(SILENT)
        Rcpp::Rcout << "OK" << std::endl;
    //queue alignments
    std::priority_queue<mySolutionShort> newTmp;
    newTmp.push(init);

    clock_t begin = clock();
    //iterate over rows/src sequences
    for(unsigned long irow = 0; irow < src_profile.size(); ++irow)
    {
        //newly generated candidates
        std::priority_queue<mySolutionShort> newTmp2;

        int candEval = 0; int candNonEval = 0;
        //for each myALIGN in OPEN
        unsigned long opentmpSIZE = newTmp.size();
        for(unsigned long ialign = 0; ialign < opentmpSIZE; ++ialign)
        {
            bool ban = false;

            mySolutionShort actSol = newTmp.top();
            newTmp.pop();
            //look at the row irow
            //make new combinations
            std::vector<mySolutionShort> mycands = makeCombinations(irow, &actSol, &commonSeqs, &src_profile, &peakPosTgz);
            //evaluate candidates
            mySolutionShort newscorezero = computeSimilarityMatrix4(mycands[0], irow, &tgt_profile, &normalizator);
            double newscorezero_maxLocalScore = newscorezero.maxLocalScore;
            //freeMatrix(&newscorezero);
            double maxscore = 0;//newscore.maxLocalScore;
            for(unsigned long isol = 0; isol < mycands.size(); ++isol)
            {

                std::vector<unsigned long> orderPenalization(mycands[isol].src.orderPenalty.size());
                unsigned long iitOrder = 0;
                for(std::list<unsigned long>::iterator itOrderPenalization = mycands[isol].src.orderPenalty.begin(); itOrderPenalization != mycands[isol].src.orderPenalty.end(); ++itOrderPenalization)
                {
                    orderPenalization[iitOrder] = (*itOrderPenalization);
                    ++iitOrder;
                }


                double param1 = orderPenalization[mycands[isol].startROW-1];//std::abs((long int)x - (long int)y);
                double sd = 300;
                double coef1 = erfc(param1 / (sd*sqrt(2)) );

                double param2 = orderPenalization[mycands[isol].endROW-1];//std::abs((long int)x - (long int)y);
                double coef2 = erfc(param2 / (sd*sqrt(2)) );
                double predictedScore = newscorezero_maxLocalScore + 10*coef1 + 10*coef2;


                if(ban)
                {
                    ++candNonEval;
                    freeMatrix(&mycands[isol]);
                }
                else
                {
                    mySolutionShort newscore = computeSimilarityMatrix4(mycands[isol], irow, &tgt_profile, &normalizator);
                    if(newscore.maxLocalScore > maxscore)
                    {
                        maxscore = newscore.maxLocalScore;

                        newTmp2.push(newscore);
                        ++candEval;
                        Rcpp::Rcout << " " << isol;
                    }
                    else
                    {
                        //for(unsigned long isolRest = isol; isolRest < mycands.size(); ++isolRest)
                        //{
                            ++candNonEval;
                            freeMatrix(&mycands[isol]);
                            maxscore = 1000000;
                            ban = true;
                        //}
                       //break;
                    }
                }

            }
            //filter candidates
            while(newTmp2.size() > FILTER_LIMIT)
            {
                mySolutionShort todelete = newTmp2.top();
                freeMatrix(&todelete);
                newTmp2.pop();
            }
            freeMatrix(&actSol);
        }
        //continue
        newTmp = newTmp2;
        while(newTmp.size() > FILTER_LIMIT)
        {
            mySolutionShort todelete = newTmp.top();
            freeMatrix(&todelete);
            newTmp.pop();
        }
        if(SILENT)
            Rcpp::Rcout << "row " << irow+1 << " / " << src_profile.size() << " evaluated: "<< candEval << " nonEvaluated: " << candNonEval << std::endl;
    }
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    //find the best
    while(newTmp.size() != 1)
    {
        mySolutionShort todelete = newTmp.top();
        freeMatrix(&todelete);
        newTmp.pop();
    }
    mySolutionShort bestnoncomplete = newTmp.top();
    mySolution best = computeSimilarityMatrix4Complete(bestnoncomplete.src.profile_mz, &tgt_profile);
    printAlignment(&best, &tgt_profile);
    //std::vector<std::string> resalign = mergeAlignments(&best, &sample2);
    freeMatrix(&bestnoncomplete);
}


int main(int argc, char* argv[])
{
    //load data
    //std::vector<std::string> sample2 = loadCSVDataMZ("/home/frantisek/QT_projects/mydamerau4/main/sample2head20_sub_sub.csv");
    //std::vector<std::string> sample1 = loadCSVDataMZ("/home/frantisek/QT_projects/mydamerau4/main/sample1head20_sub_sub.csv");


    //std::vector<std::string> sample2 = loadCSVDataMZ("/home/frantisek/QT_projects/mydamerau4/main/sample1head20_sub.csv");
    //std::vector<std::string> sample1 = loadCSVDataMZ("/home/frantisek/QT_projects/mydamerau4/main/sample2head20_sub.csv");


    unsigned long FILTER_LIMIT = 1;
    bool twocolumnResult = false;
    bool normalAling = false;
    bool normalAling2 = false;
    bool newalg = false;
     std::vector<std::string > samplePaths;

    for(int iparam = 1; iparam < argc; ++iparam)
    {

        if(strcmp(argv[iparam], "-i") == 0 && iparam+1 < argc)
        {
            Rcpp::Rcout << "-i: " << argv[iparam+1] << std::endl;
            FILTER_LIMIT = strtoul(argv[iparam+1], NULL, 0);
            ++iparam;
        }
        else if(strcmp(argv[iparam], "-v") == 0)
        {
            twocolumnResult = true;
            Rcpp::Rcout << "Two column results..." << std::endl;
        }
        else if(strcmp(argv[iparam], "-normal") == 0)
        {
            normalAling = true;
            //Rcpp::Rcout << "Two column results..." << std::endl;
        }
        else if(strcmp(argv[iparam], "-new") == 0)
        {
            newalg = true;
            //Rcpp::Rcout << "Two column results..." << std::endl;
        }
        else if(strcmp(argv[iparam], "-normal2") == 0)
        {
            normalAling2 = true;
            //Rcpp::Rcout << "Two column results..." << std::endl;
        }
        else
        {
            samplePaths.push_back(argv[iparam]);
            Rcpp::Rcout << argv[iparam] << std::endl;
        }
    }

    std::vector<std::vector<std::string > > mysamples;
    for(int ipath = 0; ipath < samplePaths.size(); ++ipath)
    {
        mysamples.push_back(loadCSVDataMZ(std::string(samplePaths[ipath])));
    }

    //vypocet comulative probability for gaussian
    //vztah z https://www.boost.org/doc/libs/1_42_0/libs/math/doc/sf_and_dist/html/math_toolkit/dist/dist_ref/dists/normal_dist.html
    //erfc z http://www.cplusplus.com/reference/cmath/erfc/
    //applet for check https://homepage.divms.uiowa.edu/~mbognar/applets/normal.html
    //double param = 15;
    //double sd = 100/3;
    //double res = erfc(param / (sd*sqrt(2)) );
    //Rcpp::Rcout << "res: " << res << std::endl;

    //return 1;


    if(newalg)
    {
        int window_size = 5;
        std::vector<std::vector<std::string > > mysamples_windows;
        for(int istart = window_size-1; istart < mysamples[0].size(); ++istart)
        {
            std::vector<std::string > actSample;
            actSample.resize(window_size);
            for(int iwin = 0; iwin < window_size; ++iwin)
            {
                actSample[(window_size-1)-iwin] = mysamples[0][istart-iwin];
            }
            mysamples_windows.push_back(actSample);
        }


        for(int isample = 0; isample < mysamples_windows.size(); ++isample)
        {
            std::list<std::string*> rowsrc_profile_mz;
            for(unsigned long irow = 0; irow < mysamples_windows[isample].size(); ++irow)
            {
                rowsrc_profile_mz.push_back(&mysamples_windows[isample][irow]);
            }
            mySolution rowAlign = computeSimilarityMatrix4CompleteLocal(rowsrc_profile_mz, &mysamples[1]);
            //printAlignment(&rowAlign, &mysamples[1]);
            std::string name;            
            std::stringstream ss;
            ss << isample;
            name = "align_" + ss.str();//std::to_string(isample);
            Rcpp::Rcout << name << std::endl;
            printLocalAlignmentToFile(&rowAlign, &mysamples[1], name);
        }


        return 1;
    }

    if(!normalAling && !normalAling2)
    {
        if(!twocolumnResult)
        {
            int round = 1;
            while(mysamples.size() != 1)
            {
                double bestscore = 0;
                int bestscoreI;
                int bestScoreJ;
                Rcpp::Rcout << "round " << round << std::endl;
                #pragma omp parallel for
                for(int i = 0; i < mysamples.size(); ++i)
                {
                    for(int j = i+1; j < mysamples.size(); ++j)
                    {
                        double score;
                        score = computeScore(mysamples[i], mysamples[j], FILTER_LIMIT);
                        Rcpp::Rcout << "score: " << samplePaths[i] << " j: " << samplePaths[j] << " : " << score << std::endl;
                        if(score > bestscore)
                        {
                            bestscore = score;
                            bestscoreI = i;
                            bestScoreJ = j;
                        }
                    }
                }

                Rcpp::Rcout << "best samples: " << samplePaths[bestscoreI] << " and " << samplePaths[bestScoreJ] << " ... will be concatenated" << std::endl;

                std::vector<std::vector<std::string > > new_mysamples;
                std::vector<std::string > new_samplePaths;
                for(int isam = 0; isam < mysamples.size(); ++isam)
                {
                    if(!(isam == bestscoreI || isam == bestScoreJ))
                    {
                        new_mysamples.push_back(mysamples[isam]);
                        new_samplePaths.push_back(samplePaths[isam]);
                    }
                }
                std::vector<std::string> newalign = makeAlignment(mysamples[bestscoreI], mysamples[bestScoreJ], FILTER_LIMIT);
                Rcpp::Rcout << "new alignment:" << std::endl;
                for(int ialign = 0; ialign < newalign.size(); ++ialign)
                {
                    Rcpp::Rcout << newalign[ialign] << std::endl;
                }

                std::string newpath = samplePaths[bestscoreI] + std::string("_") + samplePaths[bestScoreJ];
                mysamples = new_mysamples;
                samplePaths = new_samplePaths;
                mysamples.push_back(newalign);
                samplePaths.push_back(newpath);
                ++round;


                //Rcpp::Rcout << "########### END ROUND ###########" << std::endl << std::endl;
            }
        }
        else
        {
            Rcpp::Rcout << "Two columns result" << std::endl;
            make2ColAlignment(mysamples[0], mysamples[1], FILTER_LIMIT);
        }
    }
    else if(normalAling)
    {
        //Rcpp::Rcout << "=====RAW alignment: =====" << std::endl;
        std::list<std::string*> rowsrc_profile_mz;
        for(unsigned long irow = 0; irow < mysamples[0].size(); ++irow)
        {
            rowsrc_profile_mz.push_back(&mysamples[0][irow]);
        }
        mySolution rowAlign = computeSimilarityMatrix4Complete(rowsrc_profile_mz, &mysamples[1]);
        //Rcpp::Rcout << "Printing..." << std::endl;
        //printAlignment(&rowAlign, &mysamples[1]);
        printAlignmentOneColumn(&rowAlign, &mysamples[1]);
        //Rcpp::Rcout << "finished in " << elapsed_secs << "seconds" << std::endl;
        //Rcpp::Rcout << "SCORE wrt order: " << best.score << " only alignment: " << rowAlign.score << std::endl;
        //freeMatrix(&bestnoncomplete);
    }
    else
    {
        //Rcpp::Rcout << "=====RAW alignment: =====" << std::endl;
        std::list<std::string*> rowsrc_profile_mz;
        for(unsigned long irow = 0; irow < mysamples[0].size(); ++irow)
        {
            rowsrc_profile_mz.push_back(&mysamples[0][irow]);
        }
        mySolution rowAlign = computeSimilarityMatrix4Complete(rowsrc_profile_mz, &mysamples[1]);
        printAlignment(&rowAlign, &mysamples[1]);
        //printAlignmentOneColumn(&rowAlign, &mysamples[1]);
    }




    return 1;
}

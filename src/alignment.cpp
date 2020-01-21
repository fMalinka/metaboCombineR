#include "alignment.h"

#define N_NEIGHBORS 100
#define MAX_PEAK_DISTANCE 150

int s_match = 100;
int s_dismatch = -500000;//std::numeric_limits<int>::min();
int s_del = 1;
int s_ins = 1;

void printMatrix(double *matrix, int nrow, int ncol)
{
    for(int x = 0; x < nrow*ncol; ++x)
    {
        if(x % ncol ==0)
            Rcpp::Rcout << std::endl;

            Rcpp::Rcout << matrix[x] <<",";

    }
}

void freeMatrix(mySolutionShort *tofree)
{
    delete[] tofree->mymatrix;
}

mySolutionShort generateInitionSet(std::vector<std::string> *src, std::vector<std::string> *tgt, std::map<std::string, bool> *commonSeqs, std::multimap<std::string, unsigned long> *peakPosTgz)
{
    std::vector<multipleAlingmentR> newAlign(2);
    multipleAlingmentR sampleSrc;
    //sampleSrc.alignedExp = std::vector<std::vector<double> >(src.size());  //number of rows
    //create table with 0s

    sampleSrc.n_alignments = 1;
    //sampleSrc.profile_gap = std::list<double>(src->size(), 0);
    //sampleSrc.profile_rt = std::list<double>(src->size(), 0);
    //std::copy(src.begin(), src.end(), std::back_inserter(sampleSrc.profile_mz));

    //sampleSrc.profile_mz.resize(src->size());
    std::map<std::string, bool> commonSeq1;
    for(unsigned long isrc = 0; isrc < src->size(); ++isrc)
    {
        commonSeq1[(*src)[isrc]] = true;
        //sampleSrc.profile_mz.push_back(&(*src)[isrc]);
    }
    //sampleSrc.it_profile_mz = sampleSrc.profile_mz.begin();


    //std::copy(tgt->begin(), tgt->end(), std::back_inserter(tgtAlignment->profile_mz);
    newAlign[0] = sampleSrc;
    //newAlign[1] = sampleTgt;


    //generate inition alignment
    //myALIGN initAlign = generateInitialAlign(&newAlign);
    mySolutionShort init;
    init.src = newAlign[0];
    //init.tgt = newAlign[1];
    init.mymatrix = new double[(tgt->size()+1)*(N_NEIGHBORS+1)];//generateInitialAlign(&sampleSrc, tgt->size());
    //first row inicializace
    for(int x = 0; x < (tgt->size()+1); ++x)
    {
        init.mymatrix[x] = -x; //global
        //init.mymatrix[x] = 0;   //semiglobal
    }
    //first col inicializace
    double fcol = 0;
    for(int x = 0; x < (tgt->size()+1)*(N_NEIGHBORS+1); x += (tgt->size()+1))
    {
        init.mymatrix[x] = -fcol; //global
        //init.mymatrix[x] = 0;   //semiglobal
        ++fcol;
    }
    init.sizeMymatrix = (tgt->size()+1)*(N_NEIGHBORS+1);
    //init orderPenalization
    //for(int ipenalty = 0; ipenalty <= N_NEIGHBORS; ++ipenalty)
    //{
    //    init.src.orderPenalty.push_back(0); //becouse +1
    //}

//    printMatrix(init.mymatrix, (N_NEIGHBORS+1), (tgt->size()+1));


    //init.src_profile = *src;
    init.score = 0;
    //initSet.mymatrix = initAlign;

    for(unsigned long itgt = 0; itgt < tgt->size(); ++itgt)
    {
        if(commonSeq1.count((*tgt)[itgt]))
        {
            (*commonSeqs)[(*tgt)[itgt]] = true;
        }
        peakPosTgz->insert(std::pair<std::string, unsigned long>((*tgt)[itgt], itgt));
    }
    init.startROW = 0;
    init.endROW = N_NEIGHBORS + 1;

    return  init;

}

std::vector<mySolutionShort> makeCombinations(unsigned long actRow, mySolutionShort *mysol, std::map<std::string, bool> *commonSeqs, std::vector<std::string> *src_profile, std::multimap<std::string, unsigned long> *peakPosTgz)
{
    std::vector<mySolutionShort> combinations;
    std::list<std::string> mycomb;
    //make combinations

/*
    bool canSwap = true;
    //cout distances
    auto rangePeakDist = peakPosTgz->equal_range((*src_profile)[actRow]);  // It returns a pair representing the range of elements with key equal to 'c'
    unsigned long minPeakDistance = 500000;
    for (auto itPDist = rangePeakDist.first; itPDist != rangePeakDist.second; itPDist++)
    {
        unsigned long dist = std::abs((signed int)actRow - (signed int)itPDist->second);
        if(dist < minPeakDistance)
        {
            minPeakDistance = dist;
        }
    }

    if(minPeakDistance >= MAX_PEAK_DISTANCE)
        canSwap = false;
*/
    if(commonSeqs->count((*src_profile)[actRow]))// && canSwap)
    {
        //std::vector<std::string*> srcprofileVect(mysol->src.profile_mz.begin(), mysol->src.profile_mz.end());
        //zaprve jednoduse pridam nakonec, potom iterace
        double firstColOffset = 0;
        mySolutionShort myactsol = *mysol;
        unsigned long offset = 0;

        myactsol.src.profile_mz.push_back(&(*src_profile)[actRow]);

        myactsol.startROW = actRow+1;//myactsol.src.profile_mz.size();
        myactsol.endROW = actRow+1;//myactsol.src.profile_mz.size();
        //if(myactsol.src.profile_mz.size() > N_NEIGHBORS)
        if(actRow+1 > N_NEIGHBORS)
        {
            //myactsol.src.profile_mz.pop_front();
            //myactsol.mymatrix.push_back(matrix_col);
            //myactsol.mymatrix.erase(myactsol.mymatrix.begin());
            myactsol.startROW = N_NEIGHBORS;
            myactsol.endROW = N_NEIGHBORS;
            offset = mysol->sizeMymatrix/(N_NEIGHBORS+1);
            firstColOffset = mysol->mymatrix[0]-1;
        }
        myactsol.mymatrix = new double[mysol->sizeMymatrix];
        memmove(myactsol.mymatrix, mysol->mymatrix + offset, sizeof(double)*(mysol->sizeMymatrix - offset));

        if(actRow+1 > N_NEIGHBORS)
        {
            unsigned long tgtSize = mysol->sizeMymatrix/(N_NEIGHBORS+1);
            //first col inicializace
            double fcol = firstColOffset;
            for(int x = 0; x < tgtSize*(N_NEIGHBORS+1); x += tgtSize)
            {
                myactsol.mymatrix[x] = fcol;
                --fcol;
            }
        }


        //Rcpp::Rcout << "actual: " << **myactsol.src.it_profile_mz << std::endl;
        //Rcpp::Rcout << "actual2: " << **myactsol.src.profile_mz.begin() << std::endl;

        //++myactsol.src.it_profile_mz;
        myactsol.src.orderPenalty.push_back(0);
        //order penalization
        if(actRow+1 > N_NEIGHBORS)
        {
            myactsol.src.orderPenalty.erase(myactsol.src.orderPenalty.begin());
        }
        combinations.push_back(myactsol);



        for(unsigned long icomb = 1; icomb < N_NEIGHBORS && icomb < actRow+1; ++icomb)
        {
            mySolutionShort myactsol = *mysol;
            std::list<std::string*>::iterator it =  myactsol.src.profile_mz.end();

            if(actRow+1 > N_NEIGHBORS)
                myactsol.src.orderPenalty.erase(myactsol.src.orderPenalty.begin());

            std::list<unsigned long>::iterator itOrderPenalty = myactsol.src.orderPenalty.end();
            signed long indSwap = 0;
            for(unsigned long iit = 0; iit < icomb; ++iit)
            {
                --it;
                ++indSwap;
                --itOrderPenalty;
            }

            signed long swapIndex = actRow+1 - indSwap -1;//srcprofileVect.size() - indSwap -1;
            //-1 protoze chci davat novy prvek prek begin()
            //Rcpp::Rcout << "test: " << **it << std::endl;
            if(swapIndex == -1 || commonSeqs->count(**it))
            {
                myactsol.src.profile_mz.insert(it, &(*src_profile)[actRow]);
                unsigned long end = N_NEIGHBORS;
                unsigned long offset = 0;
                if(actRow+1 > N_NEIGHBORS)
                {
                    //myactsol.src.profile_mz.pop_front();
                    ;//myactsol.mymatrix.push_back(matrix_col);
                    //myactsol.mymatrix.erase(myactsol.mymatrix.begin());
                    offset = (mysol->sizeMymatrix/(N_NEIGHBORS+1));//indSwap * (mysol->sizeMymatrix/(N_NEIGHBORS+1));
                }
                else
                    end = myactsol.src.profile_mz.size();

                myactsol.startROW = end - icomb;
                myactsol.endROW = end;
                myactsol.mymatrix = new double[mysol->sizeMymatrix];
                memmove(myactsol.mymatrix, mysol->mymatrix + offset, sizeof(double)*(mysol->sizeMymatrix - offset));

                unsigned long tgtSize = mysol->sizeMymatrix/(N_NEIGHBORS+1);
                //first col inicializace
                if(actRow+1 > N_NEIGHBORS)
                {
                    double fcol = firstColOffset;
                    for(int x = 0; x < tgtSize*(N_NEIGHBORS+1); x += tgtSize)
                    {
                        myactsol.mymatrix[x] = fcol;
                        --fcol;
                    }
                }

                myactsol.src.orderPenalty.insert(itOrderPenalty, indSwap);
                //aktualizuj orderPenalty
                bool reachedPenaltylimit = false;
                for(std::list<unsigned long>::iterator itOrderPenaltyAct = itOrderPenalty; itOrderPenaltyAct != myactsol.src.orderPenalty.end(); ++itOrderPenaltyAct)
                {
                    *itOrderPenaltyAct += 1;
                    if(*itOrderPenaltyAct > MAX_PEAK_DISTANCE)
                    {
                        reachedPenaltylimit = true;
                        break;
                    }
                }


                if(!reachedPenaltylimit)
                    combinations.push_back(myactsol);
                else
                {
                    freeMatrix(&myactsol);
                }
            }
        }
    }
    else    //seq neni ve tgt sekvenci - nema smysl delat vice variant
    {
        mySolutionShort myactsol = *mysol;
        myactsol.src.profile_mz.push_back(&((*src_profile)[actRow]));
        unsigned long end = actRow+1;
        unsigned long offset = 0;
        double firstColOffset = 0;

        if(actRow+1 > N_NEIGHBORS)
        {
            offset = mysol->sizeMymatrix/(N_NEIGHBORS+1);
            firstColOffset = mysol->mymatrix[0]-1;
            end = N_NEIGHBORS;
        }
        else
            end = myactsol.src.profile_mz.size();

        myactsol.mymatrix = new double[mysol->sizeMymatrix];
        memmove(myactsol.mymatrix, mysol->mymatrix + offset, sizeof(double)*(mysol->sizeMymatrix - offset));

         unsigned long tgtSize = mysol->sizeMymatrix/(N_NEIGHBORS+1);
         //first col inicializace
         if(actRow+1 > N_NEIGHBORS)
         {
             double fcol = firstColOffset;
             for(int x = 0; x < tgtSize*(N_NEIGHBORS+1); x += tgtSize)
             {
                 myactsol.mymatrix[x] = fcol;
                 --fcol;
             }
         }

         myactsol.src.orderPenalty.push_back(0);
         //order penalization
         if(actRow+1 > N_NEIGHBORS)
         {
             myactsol.src.orderPenalty.erase(myactsol.src.orderPenalty.begin());
         }

        myactsol.endROW = end;
        myactsol.startROW = end;
        //++myactsol.src.it_profile_mz;
        combinations.push_back(myactsol);
    }

    return combinations;
}

mySolution computeSimilarityMatrix4Complete(std::list<std::string *> src_profile_mz, std::vector<std::string> *tgt_profile)
{
    //KONTROLA CO JE X A CO Y???!!!!! - dle vysledku prvni souradnice je horizontal (sloupce)


    //int dimX = first->profile_mz.size()+1;
    mySolution newSolution;// =  *src_profile_mz;
    newSolution.src.profile_mz = src_profile_mz;
    std::vector<std::vector<m_elemnt > > matrix(src_profile_mz.size()+1);
    for(unsigned long x = 0; x < src_profile_mz.size()+1; ++x)
    {
        for(unsigned long y = 0; y < tgt_profile->size()+1; ++y)
        {
            m_elemnt tmp;
            //std::vector<m_elemnt *> tmp_anc;
            //tmp.ancestor = tmp_anc;
            tmp.score = 0;
            tmp.index2s1 = 0;   //-1
            tmp.index2s2 = 0;   //-1
            tmp.direction = std::vector<int>();//from_none;
            matrix[x].push_back(tmp);//dp_matrix[y][x] = 0;
        }
    }    
    double val = 0;
    newSolution.mymatrix = matrix;
    //matrix[0][0].direction.push_back(from_none);
    for(unsigned long y = 0; y < matrix[0].size(); ++y)
    {

        //pridano
        if(y >= 1)
        {
            newSolution.mymatrix[0][y].score = val;
            newSolution.mymatrix[0][y].direction.push_back(from_top);
            newSolution.mymatrix[0][y].ancestor.push_back(&(newSolution.mymatrix[0][y-1]));
            newSolution.mymatrix[0][y].index2s2 = y-1;
        }

        --val;
    }
    val = 0;
    for(unsigned long x = 0; x < matrix.size(); ++x)
    {

        //pridano
        if(x >= 1)
        {
            newSolution.mymatrix[x][0].score = val;
            newSolution.mymatrix[x][0].direction.push_back(from_left);
            newSolution.mymatrix[x][0].ancestor.push_back(&(newSolution.mymatrix[x-1][0]));
            newSolution.mymatrix[x][0].index2s1 = x-1;
        }

        --val;
    }


    unsigned long dimY = tgt_profile->size();//second->profile_mz.size()+1;
    unsigned long dimX = newSolution.src.profile_mz.size();//second->profile_mz.size()+1;

    //std::vector<std::string> first;
    std::vector<std::string> first;
    for(std::list<std::string*>::iterator it = newSolution.src.profile_mz.begin(); it != newSolution.src.profile_mz.end(); ++it)
    {
        first.push_back(**it);
    }

    //std::vector<std::string> second(newSolution.tgt.profile_mz.begin(), newSolution.tgt.profile_mz.end());

    int shoda = 0;
    //actRow += 2;

    unsigned long x = 1;

    for(x; x <= dimX; ++x)
    {
        for(unsigned long y = 1; y <= dimY; ++y)
        {

            double plus_match = 0;// ((firstVector[x-1]==secondVector[y-1])?s_match:s_dismatch);
             if(first[x-1] == (*tgt_profile)[y-1])
             {
                 shoda++;
             }
             else
             {
                plus_match += s_dismatch;
             }

            double match = (newSolution.mymatrix)[x-1][y-1].score + plus_match;
            double del = (newSolution.mymatrix)[x][y-1].score - s_del;
            double ins = (newSolution.mymatrix)[x-1][y].score - s_ins;

            double max = std::max(match, std::max(del, ins));

            //int indcharx = x-1;
            //int indchary = y-1;

            (newSolution.mymatrix)[x][y].index2s1 = x-1;
            (newSolution.mymatrix)[x][y].index2s2 = y-1;

                (newSolution.mymatrix)[x][y].score = max;
                //match
                if(max == match)
                {
                    (newSolution.mymatrix)[x][y].direction.push_back(from_oposit);
                    (newSolution.mymatrix)[x][y].ancestor.push_back(&((newSolution.mymatrix)[x-1][y-1]));
                }

                //del
                if(max == del)
                {
                    (newSolution.mymatrix)[x][y].direction.push_back(from_top);
                    (newSolution.mymatrix)[x][y].ancestor.push_back(&((newSolution.mymatrix)[x][y-1]));
                }
                //ins
                else if(max == ins)
                {
                    (newSolution.mymatrix)[x][y].direction.push_back(from_left);
                    (newSolution.mymatrix)[x][y].ancestor.push_back(&((newSolution.mymatrix)[x-1][y]));
                }

        }
    }
    newSolution.score = (newSolution.mymatrix)[dimX][dimY].score;
    return newSolution;

}

mySolution computeSimilarityMatrix4CompleteSemiglobal(std::list<std::string *> src_profile_mz, std::vector<std::string> *tgt_profile)
{
    //KONTROLA CO JE X A CO Y???!!!!! - dle vysledku prvni souradnice je horizontal (sloupce)


    //int dimX = first->profile_mz.size()+1;
    mySolution newSolution;// =  *src_profile_mz;
    newSolution.src.profile_mz = src_profile_mz;
    std::vector<std::vector<m_elemnt > > matrix(src_profile_mz.size()+1);
    for(unsigned long x = 0; x < src_profile_mz.size()+1; ++x)
    {
        for(unsigned long y = 0; y < tgt_profile->size()+1; ++y)
        {
            m_elemnt tmp;
            //std::vector<m_elemnt *> tmp_anc;
            //tmp.ancestor = tmp_anc;
            tmp.score = 0;
            tmp.index2s1 = 0;   //-1
            tmp.index2s2 = 0;   //-1
            tmp.direction = std::vector<int>();//from_none;
            matrix[x].push_back(tmp);//dp_matrix[y][x] = 0;
        }
    }
    double val = 0;
    newSolution.mymatrix = matrix;
    //matrix[0][0].direction.push_back(from_none);
    for(unsigned long y = 0; y < matrix[0].size(); ++y)
    {

        //pridano
        if(y >= 1)
        {
            newSolution.mymatrix[0][y].score = 0;//val;
            newSolution.mymatrix[0][y].direction.push_back(from_top);
            newSolution.mymatrix[0][y].ancestor.push_back(&(newSolution.mymatrix[0][y-1]));
            newSolution.mymatrix[0][y].index2s2 = y-1;
        }

        --val;
    }
    val = 0;
    for(unsigned long x = 0; x < matrix.size(); ++x)
    {

        //pridano
        if(x >= 1)
        {
            newSolution.mymatrix[x][0].score = 0;//val;
            newSolution.mymatrix[x][0].direction.push_back(from_left);
            newSolution.mymatrix[x][0].ancestor.push_back(&(newSolution.mymatrix[x-1][0]));
            newSolution.mymatrix[x][0].index2s1 = x-1;
        }

        --val;
    }


    unsigned long dimY = tgt_profile->size();//second->profile_mz.size()+1;
    unsigned long dimX = newSolution.src.profile_mz.size();//second->profile_mz.size()+1;

    //std::vector<std::string> first;
    std::vector<std::string> first;
    for(std::list<std::string*>::iterator it = newSolution.src.profile_mz.begin(); it != newSolution.src.profile_mz.end(); ++it)
    {
        first.push_back(**it);
    }

    //std::vector<std::string> second(newSolution.tgt.profile_mz.begin(), newSolution.tgt.profile_mz.end());

    int shoda = 0;
    //actRow += 2;

    unsigned long x = 1;

    for(x; x <= dimX; ++x)
    {
        for(unsigned long y = 1; y <= dimY; ++y)
        {

            double plus_match = 0;// ((firstVector[x-1]==secondVector[y-1])?s_match:s_dismatch);
             if(first[x-1] == (*tgt_profile)[y-1])
             {
                 shoda++;
             }
             else
             {
                plus_match += s_dismatch;
             }

            double match = (newSolution.mymatrix)[x-1][y-1].score + plus_match;
            double del = (newSolution.mymatrix)[x][y-1].score - s_del;
            double ins = (newSolution.mymatrix)[x-1][y].score - s_ins;
            if(x == dimX || y == dimY)
            {
                del = (newSolution.mymatrix)[x][y-1].score;
                ins = (newSolution.mymatrix)[x-1][y].score;
            }

            double max = std::max(match, std::max(del, ins));

            //int indcharx = x-1;
            //int indchary = y-1;

            (newSolution.mymatrix)[x][y].index2s1 = x-1;
            (newSolution.mymatrix)[x][y].index2s2 = y-1;

                (newSolution.mymatrix)[x][y].score = max;
                //match
                if(max == match)
                {
                    (newSolution.mymatrix)[x][y].direction.push_back(from_oposit);
                    (newSolution.mymatrix)[x][y].ancestor.push_back(&((newSolution.mymatrix)[x-1][y-1]));
                }

                //del
                if(max == del)
                {
                    (newSolution.mymatrix)[x][y].direction.push_back(from_top);
                    (newSolution.mymatrix)[x][y].ancestor.push_back(&((newSolution.mymatrix)[x][y-1]));
                }
                //ins
                else if(max == ins)
                {
                    (newSolution.mymatrix)[x][y].direction.push_back(from_left);
                    (newSolution.mymatrix)[x][y].ancestor.push_back(&((newSolution.mymatrix)[x-1][y]));
                }

                if((newSolution.mymatrix)[x][y].score < 0)
                    (newSolution.mymatrix)[x][y].score = 0;

        }
    }
    newSolution.score = (newSolution.mymatrix)[dimX][dimY].score;
    return newSolution;

}

mySolution computeSimilarityMatrix4CompleteLocal(std::list<std::string *> src_profile_mz, std::vector<std::string> *tgt_profile)
{
    //KONTROLA CO JE X A CO Y???!!!!! - dle vysledku prvni souradnice je horizontal (sloupce)


    //int dimX = first->profile_mz.size()+1;
    mySolution newSolution;// =  *src_profile_mz;
    newSolution.src.profile_mz = src_profile_mz;
    std::vector<std::vector<m_elemnt > > matrix(src_profile_mz.size()+1);
    for(unsigned long x = 0; x < src_profile_mz.size()+1; ++x)
    {
        for(unsigned long y = 0; y < tgt_profile->size()+1; ++y)
        {
            m_elemnt tmp;
            //std::vector<m_elemnt *> tmp_anc;
            //tmp.ancestor = tmp_anc;
            tmp.score = 0;
            tmp.index2s1 = 0;   //-1
            tmp.index2s2 = 0;   //-1
            tmp.direction = std::vector<int>();//from_none;
            matrix[x].push_back(tmp);//dp_matrix[y][x] = 0;
        }
    }
    double val = 0;
    newSolution.mymatrix = matrix;
    //matrix[0][0].direction.push_back(from_none);
    for(unsigned long y = 0; y < matrix[0].size(); ++y)
    {

        //pridano
        if(y >= 1)
        {
            newSolution.mymatrix[0][y].score = 0;//val;
            newSolution.mymatrix[0][y].direction.push_back(from_top);
            newSolution.mymatrix[0][y].ancestor.push_back(&(newSolution.mymatrix[0][y-1]));
            newSolution.mymatrix[0][y].index2s2 = y-1;
        }

        --val;
    }
    val = 0;
    for(unsigned long x = 0; x < matrix.size(); ++x)
    {

        //pridano
        if(x >= 1)
        {
            newSolution.mymatrix[x][0].score = 0;//val;
            newSolution.mymatrix[x][0].direction.push_back(from_left);
            newSolution.mymatrix[x][0].ancestor.push_back(&(newSolution.mymatrix[x-1][0]));
            newSolution.mymatrix[x][0].index2s1 = x-1;
        }

        --val;
    }


    unsigned long dimY = tgt_profile->size();//second->profile_mz.size()+1;
    unsigned long dimX = newSolution.src.profile_mz.size();//second->profile_mz.size()+1;

    //Rcpp::Rcout << "dim: tgt: " << dimY << " src:" << dimX << std::endl;

    //std::vector<std::string> first;
    std::vector<std::string> first;
    for(std::list<std::string*>::iterator it = newSolution.src.profile_mz.begin(); it != newSolution.src.profile_mz.end(); ++it)
    {
        first.push_back(**it);
    }

    //std::vector<std::string> second(newSolution.tgt.profile_mz.begin(), newSolution.tgt.profile_mz.end());

    int shoda = 0;
    //actRow += 2;

    unsigned long x = 1;

    for(x; x <= dimX; ++x)
    {
        for(unsigned long y = 1; y <= dimY; ++y)
        {

            double plus_match = 0;// ((firstVector[x-1]==secondVector[y-1])?s_match:s_dismatch);
             if(first[x-1] == (*tgt_profile)[y-1])
             {
                 shoda++;
                 plus_match += s_match;
             }
             else
             {
                plus_match += s_dismatch;
             }

            double match = (newSolution.mymatrix)[x-1][y-1].score + plus_match;
            double del = (newSolution.mymatrix)[x][y-1].score - s_del;
            double ins = (newSolution.mymatrix)[x-1][y].score - s_ins;

            double max = std::max(match, std::max(del, ins));

            //int indcharx = x-1;
            //int indchary = y-1;

            (newSolution.mymatrix)[x][y].index2s1 = x-1;
            (newSolution.mymatrix)[x][y].index2s2 = y-1;

                (newSolution.mymatrix)[x][y].score = max;
                //match
                if(max == match)
                {
                    (newSolution.mymatrix)[x][y].direction.push_back(from_oposit);
                    (newSolution.mymatrix)[x][y].ancestor.push_back(&((newSolution.mymatrix)[x-1][y-1]));
                }

                //del
/*                if(max == del)
                {
                    (newSolution.mymatrix)[x][y].direction.push_back(from_top);
                    (newSolution.mymatrix)[x][y].ancestor.push_back(&((newSolution.mymatrix)[x][y-1]));
                }
                //ins
                else if(max == ins)
                {
                    (newSolution.mymatrix)[x][y].direction.push_back(from_left);
                    (newSolution.mymatrix)[x][y].ancestor.push_back(&((newSolution.mymatrix)[x-1][y]));
                }
  */            if(max == ins)
                {
                    (newSolution.mymatrix)[x][y].direction.push_back(from_left);
                    (newSolution.mymatrix)[x][y].ancestor.push_back(&((newSolution.mymatrix)[x-1][y]));
                }
                else if(max == del)
                {
                    (newSolution.mymatrix)[x][y].direction.push_back(from_top);
                    (newSolution.mymatrix)[x][y].ancestor.push_back(&((newSolution.mymatrix)[x][y-1]));
                }


                //if((newSolution.mymatrix)[x][y].score > 0)
                //    Rcpp::Rcout << "score: " << (newSolution.mymatrix)[x][y].score << std::endl;
                if((newSolution.mymatrix)[x][y].score < 0)
                    (newSolution.mymatrix)[x][y].score = 0;

        }
    }
    //newSolution.score = (newSolution.mymatrix)[dimX][dimY].score;
    return newSolution;

}

mySolutionShort computeSimilarityMatrix4(mySolutionShort mysol, unsigned long actRow, std::vector<std::string> *tgt_profile, double *normalizator)
{
    unsigned long dimY = tgt_profile->size();//tgt->profile_mz.size();//second->profile_mz.size()+1;

    double maxLocalScore = -(double)(mysol.src.profile_mz.size() + tgt_profile->size());
    double maxRowSum = 0;
    //Rcpp::Rcout << "BEFORE:" << std::endl;
    //printMatrix(mysol.mymatrix, N_NEIGHBORS+1, tgt_profile->size()+1);
    //Rcpp::Rcout << "#####" << std::endl;

    std::list<std::string*> profile = mysol.src.profile_mz;
    //ponechame jenom maximalne 10 poslednich
    while(profile.size() > N_NEIGHBORS)
        profile.pop_front();

    std::vector<std::string> first(profile.size());
    unsigned long iit = 0;
    for(std::list<std::string*>::iterator it = profile.begin(); it != profile.end(); ++it)
    {
        first[iit] = (**it);
        ++iit;
    }

    std::vector<unsigned long> orderPenalization(mysol.src.orderPenalty.size());
    unsigned long iitOrder = 0;
    for(std::list<unsigned long>::iterator itOrderPenalization = mysol.src.orderPenalty.begin(); itOrderPenalization != mysol.src.orderPenalty.end(); ++itOrderPenalization)
    {
        orderPenalization[iitOrder] = (*itOrderPenalization);
        ++iitOrder;
    }
    actRow += 2;
    unsigned long ysize = dimY+1;

    /*
    Rcpp::Rcout << "PRINT matrix:" << std::endl;
    printMatrix(mysol.mymatrix, (N_NEIGHBORS+1), (tgt_profile->size()+1));
    Rcpp::Rcout << std::endl << std::endl;
*/
    //unsigned long x = 1;
    unsigned long xzero = (mysol.startROW-1)%N_NEIGHBORS;
    for(unsigned long x = mysol.startROW; x <= mysol.endROW; ++x)
    {
        for(unsigned long y = 1; y <= dimY; ++y)
        {
            //Rcpp::Rcout << "x: " << x << " y: " << y << std::endl;
            bool matchFlag = false;
            double orderPenalization_norm = orderPenalization[xzero]/(((double)*normalizator));
            double plus_match = 0;// ((firstVector[x-1]==secondVector[y-1])?s_match:s_dismatch);
             if(first[x-1] == (*tgt_profile)[y-1])
             {
                 matchFlag = true;
                 //plus_match = s_match;// - std::abs((long int)x - (long int)y)/(*normalizator);
                 //plus_match = s_match - orderPenalization_norm;
                 // plus_match = s_match - (std::abs((long int)x - (long int)y)/(MAX_PEAK_DISTANCE))*s_match;

                 double param = orderPenalization[xzero];//std::abs((long int)x - (long int)y);
                 double sd = 300;
                 double coef = erfc(param / (sd*sqrt(2)) );

                 plus_match = s_match * coef;
             }
             else
             {
                plus_match += s_dismatch;
             }


//orderPenalization_norm /= 100;
             double match = (mysol.mymatrix)[((x-1)*ysize)+(y-1)] + plus_match;
             double del = (mysol.mymatrix)[(x*ysize) + (y-1)] - s_del - orderPenalization_norm;
             double ins = (mysol.mymatrix)[((x-1)*ysize) + y] - s_ins - orderPenalization_norm;

            unsigned long index = (x * ysize) + y;
            double max = std::max(match, std::max(del, ins));
            //if(matchFlag)
            //    max = match;
            (mysol.mymatrix)[index] = max;

            if(x ==mysol.endROW)
            {
                if(max > maxLocalScore)
                {
                    maxLocalScore = max;
                }
                maxRowSum += max * erfc(std::abs((long int)x - (long int)y) / (300*sqrt(2)) );
            }
        }
        ++xzero;
    }

    mysol.score = (mysol.mymatrix)[(mysol.endROW) * ysize + dimY];
    mysol.maxLocalScore = maxLocalScore;//maxRowSum;//mysol.score;//maxLocalScore;

    /*
    Rcpp::Rcout << "PRINT matrix2222:" << std::endl;
        printMatrix(mysol.mymatrix, (N_NEIGHBORS+1), (tgt_profile->size()+1));
    Rcpp::Rcout << std::endl << std::endl;
    */

    //Rcpp::Rcout << "AFTER:" << std::endl;
    //printMatrix(mysol.mymatrix, N_NEIGHBORS+1, tgt_profile->size()+1);
    //Rcpp::Rcout << "####" << std::endl;
    return mysol;

}


void printAlignment(mySolution *mysol, std::vector<std::string> *tgt)
{
    std::vector<std::string> seq1;
    std::vector<std::string> seq2;
    std::vector<std::string> src;
    for(std::list<std::string*>::iterator it = mysol->src.profile_mz.begin(); it != mysol->src.profile_mz.end(); ++it)
    {
        src.push_back(**it);
    }

    m_elemnt element = mysol->mymatrix[mysol->src.profile_mz.size()][tgt->size()];

    int iancestor = 0;
    while(!element.ancestor.empty())
    {
        if(element.direction[iancestor] == from_oposit)
        {
            seq1.push_back(src[element.index2s1]);
            seq2.push_back((*tgt)[element.index2s2]);
        }
        else if(element.direction[iancestor] == from_top)
        {
            seq1.push_back("-");
            seq2.push_back((*tgt)[element.index2s2]);
        }
        else
        {
            seq1.push_back(src[element.index2s1]);
            seq2.push_back("-");
        }

        if((element.ancestor.empty()))
            break;
        else
            element = *(element.ancestor[iancestor]);
    }
    std::reverse(seq1.begin(), seq1.end());
    std::reverse(seq2.begin(), seq2.end());
//horizontal
//    for(int ii = 0; ii < seq1.size(); ++ii)
//    {
//        Rcpp::Rcout << seq1[ii] << " ";
//    }
//    Rcpp::Rcout << std::endl;
//    for(int ii = 0; ii < seq2.size(); ++ii)
//    {
//        Rcpp::Rcout << seq2[ii] << " ";
//    }
//    Rcpp::Rcout << std::endl;

    //vertikal
//    Rcpp::Rcout << "VERTIKAL: " << std::endl;
    for(int ii = 0; ii < seq1.size(); ++ii)
    {
        Rcpp::Rcout << seq1[ii] << "," << seq2[ii];
        Rcpp::Rcout << std::endl;
    }

}

void printAlignmentToFile(mySolution *mysol, std::vector<std::string> *tgt, std::string filename)
{
    std::vector<std::string> seq1;
    std::vector<std::string> seq2;
    std::vector<std::string> src;

    std::ofstream myfile;
    myfile.open(filename.c_str());

    for(std::list<std::string*>::iterator it = mysol->src.profile_mz.begin(); it != mysol->src.profile_mz.end(); ++it)
    {
        src.push_back(**it);
    }

    m_elemnt element = mysol->mymatrix[mysol->src.profile_mz.size()][tgt->size()];

    int iancestor = 0;
    while(!element.ancestor.empty())
    {
        if(element.direction[iancestor] == from_oposit)
        {
            seq1.push_back(src[element.index2s1]);
            seq2.push_back((*tgt)[element.index2s2]);
        }
        else if(element.direction[iancestor] == from_top)
        {
            seq1.push_back("-");
            seq2.push_back((*tgt)[element.index2s2]);
        }
        else
        {
            seq1.push_back(src[element.index2s1]);
            seq2.push_back("-");
        }


        if((element.ancestor.empty()))
            break;
        else
        {
            element = *(element.ancestor[iancestor]);
        }
    }
    std::reverse(seq1.begin(), seq1.end());
    std::reverse(seq2.begin(), seq2.end());
//horizontal
//    for(int ii = 0; ii < seq1.size(); ++ii)
//    {
//        Rcpp::Rcout << seq1[ii] << " ";
//    }
//    Rcpp::Rcout << std::endl;
//    for(int ii = 0; ii < seq2.size(); ++ii)
//    {
//        Rcpp::Rcout << seq2[ii] << " ";
//    }
//    Rcpp::Rcout << std::endl;

    //vertikal
//    Rcpp::Rcout << "VERTIKAL: " << std::endl;
    for(int ii = 0; ii < seq1.size(); ++ii)
    {
        myfile << seq1[ii] << "," << seq2[ii];
        myfile << std::endl;
    }
    myfile.close();

}

void printSemiAlignmentToFile(mySolution *mysol, std::vector<std::string> *tgt, std::string filename)
{
    std::vector<std::string> seq1;
    std::vector<std::string> seq2;
    std::vector<std::string> src;

    std::ofstream myfile;
    myfile.open(filename.c_str());

    for(std::list<std::string*>::iterator it = mysol->src.profile_mz.begin(); it != mysol->src.profile_mz.end(); ++it)
    {
        src.push_back(**it);
    }

    m_elemnt element = mysol->mymatrix[mysol->src.profile_mz.size()][tgt->size()];

    int iancestor = 0;
    while(!element.ancestor.empty())
    {
        if(element.direction[iancestor] == from_left)
        {
            seq1.push_back(src[element.index2s1]);
            seq2.push_back("-");
        }
        else if(element.direction[iancestor] == from_oposit)
        {
            seq1.push_back(src[element.index2s1]);
            seq2.push_back((*tgt)[element.index2s2]);
        }
        else if(element.direction[iancestor] == from_top)
        {
            seq1.push_back("-");
            seq2.push_back((*tgt)[element.index2s2]);
        }
        else
        {
            //seq1.push_back(src[element.index2s1]);
            //seq2.push_back("-");
        }


        if((element.ancestor.empty()))
            break;
        else
        {
            element = *(element.ancestor[iancestor]);
        }
    }
    std::reverse(seq1.begin(), seq1.end());
    std::reverse(seq2.begin(), seq2.end());
//horizontal
//    for(int ii = 0; ii < seq1.size(); ++ii)
//    {
//        Rcpp::Rcout << seq1[ii] << " ";
//    }
//    Rcpp::Rcout << std::endl;
//    for(int ii = 0; ii < seq2.size(); ++ii)
//    {
//        Rcpp::Rcout << seq2[ii] << " ";
//    }
//    Rcpp::Rcout << std::endl;

    //vertikal
//    Rcpp::Rcout << "VERTIKAL: " << std::endl;
    for(int ii = 0; ii < seq1.size(); ++ii)
    {
        myfile << seq1[ii] << "," << seq2[ii];
        myfile << std::endl;
    }
    myfile.close();

}

void printLocalAlignmentToFile(mySolution *mysol, std::vector<std::string> *tgt, std::string filename)
{
    std::vector<std::string> seq1;
    std::vector<std::string> seq2;
    std::vector<std::string> src;

    std::ofstream myfile;
    myfile.open(filename.c_str());

    for(std::list<std::string*>::iterator it = mysol->src.profile_mz.begin(); it != mysol->src.profile_mz.end(); ++it)
    {
        src.push_back(**it);
    }

    //find the best
    double localbest = -1000000;
    unsigned long besttgt = 0;
    m_elemnt element;
    for(int ibest = 0; ibest <= tgt->size(); ++ibest)
    {
        if(mysol->mymatrix[mysol->src.profile_mz.size()][ibest].score > localbest)
        {
            localbest = mysol->mymatrix[mysol->src.profile_mz.size()][ibest].score;
            element = mysol->mymatrix[mysol->src.profile_mz.size()][ibest];
            besttgt = ibest;
        }
    }

    //m_elemnt element = mysol->mymatrix[mysol->src.profile_mz.size()][tgt->size()];
    //print
    //for(unsigned int itgt = tgt->size()-1; itgt >= besttgt; --itgt)
    for(unsigned int itgt = tgt->size()-1; itgt >= besttgt; --itgt)
    {
        seq1.push_back("--");
        seq2.push_back((*tgt)[itgt]);
    }

    int iancestor = 0;
    unsigned int printedSrc = 0;
    while(!element.ancestor.empty())
    {
        if(element.direction[iancestor] == from_oposit)
        {
            seq1.push_back(src[element.index2s1]);
            seq2.push_back((*tgt)[element.index2s2]);
            ++printedSrc;
        }
        else if(element.direction[iancestor] == from_top)
        {
            seq1.push_back("-");
            seq2.push_back((*tgt)[element.index2s2]);
        }
        else
        {
            seq1.push_back(src[element.index2s1]);
            seq2.push_back("-");
            ++printedSrc;
        }

        if((element.ancestor.empty()))
            break;
        else
        {
            element = *(element.ancestor[iancestor]);
        }
    }
    std::reverse(seq1.begin(), seq1.end());
    std::reverse(seq2.begin(), seq2.end());
    for(int ii = 0; ii < seq1.size(); ++ii)
    {
        myfile << seq1[ii] << "," << seq2[ii];
        myfile << std::endl;
    }
    myfile.close();

}

alingmentIndexes processLocalAlignments(mySolution *mysol, std::vector<std::string> *tgt)
{
    std::vector<int> seq1;
    std::vector<int> seq2;
    std::vector<std::string> src;

    //std::ofstream myfile;
    //myfile.open(filename.c_str());

    for(std::list<std::string*>::iterator it = mysol->src.profile_mz.begin(); it != mysol->src.profile_mz.end(); ++it)
    {
        src.push_back(**it);
    }

    //find the best
    double localbest = -1000000;
    int besttgt = 0;
    m_elemnt element;
//    Rcpp::Rcout << "A_";
    //for(int ibest = 0; ibest < tgt->size(); ++ibest)
    for(int ibest = 0; ibest <= tgt->size(); ++ibest)
    {
        if(mysol->mymatrix[mysol->src.profile_mz.size()][ibest].score > localbest)
        {
            localbest = mysol->mymatrix[mysol->src.profile_mz.size()][ibest].score;
            element = mysol->mymatrix[mysol->src.profile_mz.size()][ibest];
            besttgt = ibest;
        }
    }
//    Rcpp::Rcout << "_B, best2: " << besttgt << std::endl;


    //m_elemnt element = mysol->mymatrix[mysol->src.profile_mz.size()][tgt->size()];
    //print
    //for(unsigned int itgt = tgt->size()-1; itgt >= besttgt; --itgt)
    //for(unsigned int itgt = tgt->size()-1; itgt > besttgt; --itgt)
    //rovna se zpusobuje zpomaleni!!!
    //Rcpp::Rcout << "my besttgt: " << besttgt << " tgt: " << tgt->size()-1 << std::endl;
    for(int itgt = tgt->size()-1; itgt >= besttgt; --itgt)
    {
        //seq1.push_back("-");
        //seq2.push_back((*tgt)[itgt]);
    //    Rcpp::Rcout << "itg: " << itgt << std::endl;
        seq1.push_back(-1);
        seq2.push_back(itgt);
    }
    //Rcpp::Rcout << "-------" << std::endl;

//    Rcpp::Rcout << "_C" << std::endl;

    //match have not found
    if(besttgt != 0)
    {
        int iancestor = 0;
        while(!element.ancestor.empty())
        {
            if(element.direction[iancestor] == from_oposit)
            {
                //seq1.push_back(src[element.index2s1]);
                //seq2.push_back((*tgt)[element.index2s2]);
                seq1.push_back(element.index2s1);
                seq2.push_back(element.index2s2);
            }
            else if(element.direction[iancestor] == from_top)
            {
                //seq1.push_back("-");
                //seq2.push_back((*tgt)[element.index2s2]);
                seq1.push_back(-1);
                seq2.push_back(element.index2s2);
            }
            else
            {
                //seq1.push_back(src[element.index2s1]);
                //seq2.push_back("-");
                seq1.push_back(element.index2s1);
                seq2.push_back(-1);
            }

            if((element.ancestor.empty()))
                break;
            else
            {
                element = *(element.ancestor[iancestor]);
            }
        }
    }

    //Rcpp::Rcout << "_D_" << std::endl;
    std::reverse(seq1.begin(), seq1.end());
    std::reverse(seq2.begin(), seq2.end());
    alingmentIndexes myAlign;
    myAlign.indexNonref = seq1;
    myAlign.indexRef = seq2;

 //   for(int ii = 0; ii < seq1.size(); ++ii)
 //   {
        //myfile << seq1[ii] << "," << seq2[ii];
        //myfile << std::endl;
 //      Rcpp::Rcout << "test: " << seq1[ii] << "," << seq2[ii] << std::endl;
   // }

//    Rcpp::Rcout << "_D_" << std::endl;
    return myAlign;
 /*   for(int ii = 0; ii < seq1.size(); ++ii)
    {
        //myfile << seq1[ii] << "," << seq2[ii];
        //myfile << std::endl;
        Rcpp::Rcout << "test: " << seq1[ii] << "," << seq2[ii] << std::endl;
    }
    */
    //myfile.close();

}

void printAlignmentOneColumn(mySolution *mysol, std::vector<std::string> *tgt)
{
    std::vector<std::string> seq1;
    std::vector<std::string> seq2;
    std::vector<std::string> src;
    for(std::list<std::string*>::iterator it = mysol->src.profile_mz.begin(); it != mysol->src.profile_mz.end(); ++it)
    {
        src.push_back(**it);
    }

    m_elemnt element = mysol->mymatrix[mysol->src.profile_mz.size()][tgt->size()];

    int iancestor = 0;
    while(!element.ancestor.empty())
    {
        if(element.direction[iancestor] == from_oposit)
        {
            seq1.push_back(src[element.index2s1]);
            seq2.push_back((*tgt)[element.index2s2]);
        }
        else if(element.direction[iancestor] == from_top)
        {
            seq1.push_back("-");
            seq2.push_back((*tgt)[element.index2s2]);
        }
        else
        {
            seq1.push_back(src[element.index2s1]);
            seq2.push_back("-");
        }

        if((element.ancestor.empty()))
            break;
        else
            element = *(element.ancestor[iancestor]);
    }
    std::reverse(seq1.begin(), seq1.end());
    std::reverse(seq2.begin(), seq2.end());
//horizontal
//    for(int ii = 0; ii < seq1.size(); ++ii)
//    {
//        Rcpp::Rcout << seq1[ii] << " ";
//    }
//    Rcpp::Rcout << std::endl;
//    for(int ii = 0; ii < seq2.size(); ++ii)
//   {
//        Rcpp::Rcout << seq2[ii] << " ";
//    }
//    Rcpp::Rcout << std::endl;

    //vertikal
//    Rcpp::Rcout << "VERTIKAL: " << std::endl;
    for(int ii = 0; ii < seq1.size(); ++ii)
    {
        if(seq1[ii] == "-")
            Rcpp::Rcout << seq2[ii] << std::endl;
        else
            Rcpp::Rcout << seq1[ii] << std::endl;
        //Rcpp::Rcout << seq1[ii] << "," << seq2[ii];
        //Rcpp::Rcout << std::endl;
    }

}

std::vector<std::string> mergeAlignments(mySolution *mysol, std::vector<std::string> *tgt)
{
    std::vector<std::string> seq1;
    std::vector<std::string> seq2;
    std::vector<std::string> src;
    for(std::list<std::string*>::iterator it = mysol->src.profile_mz.begin(); it != mysol->src.profile_mz.end(); ++it)
    {
        src.push_back(**it);
    }

    m_elemnt element = mysol->mymatrix[mysol->src.profile_mz.size()][tgt->size()];

    int iancestor = 0;
    while(!element.ancestor.empty())
    {
        if(element.direction[iancestor] == from_oposit)
        {
            seq1.push_back(src[element.index2s1]);
            seq2.push_back((*tgt)[element.index2s2]);
        }
        else if(element.direction[iancestor] == from_top)
        {
            seq1.push_back("-");
            seq2.push_back((*tgt)[element.index2s2]);
        }
        else
        {
            seq1.push_back(src[element.index2s1]);
            seq2.push_back("-");
        }

        if((element.ancestor.empty()))
            break;
        else
            element = *(element.ancestor[iancestor]);
    }
    std::reverse(seq1.begin(), seq1.end());
    std::reverse(seq2.begin(), seq2.end());
    //vertikal

    std::vector<std::string> final(seq1.size());
    for(int ii = 0; ii < seq1.size(); ++ii)
    {
        if(seq1[ii] == "-")
            final[ii] = seq2[ii];
        else
            final[ii] = seq1[ii];
    }
    return final;
}

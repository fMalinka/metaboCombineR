/**
    metaboCombineR, alignment.cpp
    Purpose: for alignment representation

    @author Frantisek Malinka
    @version 0.99.0 22/01/2020
*/
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

//computes similarity matrix and returns it
mySolution computeSimilarityMatrix4CompleteLocal(std::list<std::string *> src_profile_mz, std::vector<std::string> *tgt_profile)
{
    mySolution newSolution;// =  *src_profile_mz;
    newSolution.src.profile_mz = src_profile_mz;
    std::vector<std::vector<m_elemnt > > matrix(src_profile_mz.size()+1);
    for(unsigned long x = 0; x < src_profile_mz.size()+1; ++x)
    {
        for(unsigned long y = 0; y < tgt_profile->size()+1; ++y)
        {
            m_elemnt tmp;
            tmp.score = 0;
            tmp.index2s1 = 0;   //-1
            tmp.index2s2 = 0;   //-1
            tmp.direction = std::vector<int>();//from_none;
            matrix[x].push_back(tmp);
        }
    }
    double val = 0;
    newSolution.mymatrix = matrix;
    for(unsigned long y = 0; y < matrix[0].size(); ++y)
    {
        //added
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
        //added
        if(x >= 1)
        {
            newSolution.mymatrix[x][0].score = 0;//val;
            newSolution.mymatrix[x][0].direction.push_back(from_left);
            newSolution.mymatrix[x][0].ancestor.push_back(&(newSolution.mymatrix[x-1][0]));
            newSolution.mymatrix[x][0].index2s1 = x-1;
        }
        --val;
    }


    unsigned long dimY = tgt_profile->size();
    unsigned long dimX = newSolution.src.profile_mz.size();

    std::vector<std::string> first;
    for(std::list<std::string*>::iterator it = newSolution.src.profile_mz.begin(); it != newSolution.src.profile_mz.end(); ++it)
    {
        first.push_back(**it);
    }


    int shoda = 0;   
    unsigned long x = 1;

    for(x; x <= dimX; ++x)
    {
        for(unsigned long y = 1; y <= dimY; ++y)
        {
            double plus_match = 0;
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
            (newSolution.mymatrix)[x][y].index2s1 = x-1;
            (newSolution.mymatrix)[x][y].index2s2 = y-1;

                (newSolution.mymatrix)[x][y].score = max;
                //match
                if(max == match)
                {
                    (newSolution.mymatrix)[x][y].direction.push_back(from_oposit);
                    (newSolution.mymatrix)[x][y].ancestor.push_back(&((newSolution.mymatrix)[x-1][y-1]));
                }
                //insertion
                if(max == ins)
                {
                    (newSolution.mymatrix)[x][y].direction.push_back(from_left);
                    (newSolution.mymatrix)[x][y].ancestor.push_back(&((newSolution.mymatrix)[x-1][y]));
                }
                //deletion
                else if(max == del)
                {
                    (newSolution.mymatrix)[x][y].direction.push_back(from_top);
                    (newSolution.mymatrix)[x][y].ancestor.push_back(&((newSolution.mymatrix)[x][y-1]));
                }

                if((newSolution.mymatrix)[x][y].score < 0)
                    (newSolution.mymatrix)[x][y].score = 0;
        }
    }
    return newSolution;

}

//find the best path, follow the score
alingmentIndexes processLocalAlignments(mySolution *mysol, std::vector<std::string> *tgt)
{
    std::vector<int> seq1;
    std::vector<int> seq2;
    std::vector<std::string> src;

    for(std::list<std::string*>::iterator it = mysol->src.profile_mz.begin(); it != mysol->src.profile_mz.end(); ++it)
    {
        src.push_back(**it);
    }

    //find the best
    double localbest = DBL_MIN;
    int besttgt = 0;
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
    for(int itgt = tgt->size()-1; itgt >= besttgt; --itgt)
    {
        seq1.push_back(-1);
        seq2.push_back(itgt);
    }    

    //match have not found
    if(besttgt != 0)
    {
        int iancestor = 0;
        while(!element.ancestor.empty())
        {
            if(element.direction[iancestor] == from_oposit)
            {
                seq1.push_back(element.index2s1);
                seq2.push_back(element.index2s2);
            }
            else if(element.direction[iancestor] == from_top)
            {
                seq1.push_back(-1);
                seq2.push_back(element.index2s2);
            }
            else
            {
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

    std::reverse(seq1.begin(), seq1.end());
    std::reverse(seq2.begin(), seq2.end());
    alingmentIndexes myAlign;
    myAlign.indexNonref = seq1;
    myAlign.indexRef = seq2;
    return myAlign;
}

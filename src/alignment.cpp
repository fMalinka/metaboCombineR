/**
    metaboCombineR, alignment.cpp
    Purpose: for alignment representation

    @author Frantisek Malinka
    @version 0.99.0 22/01/2020
*/
#include "alignment.h"


/**
 * @brief Return true if features overlap, otherwise false
 * @param mz1_1 bottom mz value of first feature
 * @param mz1_2 upper mz value of first feature
 * @param mz2_1 bottom mz value of second feature
 * @param mz2_2 upper mz value of second feature
 * @return boolean
 */
bool isMzInMatch(double mz1_1, double mz1_2, double mz2_1, double mz2_2)
{
    return((mz1_1 <= mz2_2 && mz2_1 <= mz1_2)? true : false);
}

/**
 * @brief Make a global alignment
 * @param src source experiment
 * @param tgt target experiment
 * @param s_match a score for matched features
 * @param s_del a penalization for deletion of features
 * @param s_ins a penalization for insertion of features
 * @return global alignment
 */
mySolution globalAlignment(experiment *src, experiment *tgt, double s_match, double s_del, double s_ins)
{
    mySolution newSolution;
    for(int istring = 0; istring < src->mzStr.size(); ++istring)
    {
        newSolution.src.profile_mz.push_back(&(src->mzStr[istring]));
    }
    std::vector<std::vector<m_elemnt > > matrix(src->mzStr.size()+1);
    for(unsigned long x = 0; x < src->mzStr.size()+1; ++x)
    {
        for(unsigned long y = 0; y < tgt->mzStr.size()+1; ++y)
        {
            m_elemnt tmp;
            tmp.score = 0;
            tmp.index2s1 = 0;
            tmp.index2s2 = 0;
            tmp.direction = std::vector<int>();
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
        //added
        if(x >= 1)
        {
            newSolution.mymatrix[x][0].score = val;
            newSolution.mymatrix[x][0].direction.push_back(from_left);
            newSolution.mymatrix[x][0].ancestor.push_back(&(newSolution.mymatrix[x-1][0]));
            newSolution.mymatrix[x][0].index2s1 = x-1;
        }
        --val;
    }


    unsigned long dimY = tgt->mzStr.size();
    unsigned long dimX = newSolution.src.profile_mz.size();

    //int shoda = 0;
    unsigned long x = 1;

    for(; x <= dimX; ++x)
    {
        for(unsigned long y = 1; y <= dimY; ++y)
        {            
            double match = -DBL_MAX;
            if(isMzInMatch(src->mz_lower_bound[x-1],  src->mz_upper_bound[x-1], tgt->mz_lower_bound[y-1], tgt->mz_upper_bound[y-1]))
            {
             //in match
             match = (newSolution.mymatrix)[x-1][y-1].score;
            }

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


/**
 * @brief Make a local alignment
 * @param kmer kmer experiment
 * @param tgt_profile target experiment
 * @param s_match a score for matched features
 * @param s_del a penalization for deletion of features
 * @param s_ins a penalization for insertion of features
 * @return global alignment
 */
mySolution localAlignment(experiment_kmer *kmer, experiment *tgt_profile, double s_match, double s_del, double s_ins)
{
    mySolution newSolution;
    for(int istring = 0; istring < kmer->mzStr.size(); ++istring)
    {
        newSolution.src.profile_mz.push_back(&(kmer->mzStr[istring]));
    }
    std::vector<std::vector<m_elemnt > > matrix(kmer->mzStr.size()+1);
    for(unsigned long x = 0; x < kmer->mzStr.size()+1; ++x)
    {
        for(unsigned long y = 0; y < tgt_profile->mzStr.size()+1; ++y)
        {
            m_elemnt tmp;
            tmp.score = 0;
            tmp.index2s1 = 0;
            tmp.index2s2 = 0;
            tmp.direction = std::vector<int>();
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

    unsigned long dimY = tgt_profile->mzStr.size();
    unsigned long dimX = newSolution.src.profile_mz.size();


    unsigned long x = 1;
    for(; x <= dimX; ++x)
    {
        for(unsigned long y = 1; y <= dimY; ++y)
        {
            double match = -DBL_MAX;    //dismatch
            if(isMzInMatch(kmer->mz_lower_bound[x-1],  kmer->mz_upper_bound[x-1], tgt_profile->mz_lower_bound[y-1], tgt_profile->mz_upper_bound[y-1]))
            {
                //in match
                match = (newSolution.mymatrix)[x-1][y-1].score + s_match;
            }

            //double match = (newSolution.mymatrix)[x-1][y-1].score + plus_match;
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
/**
 * @brief Find a path with the best score in terms of DP, for kmersAlignment algorithm
 * @param mysol solution structure
 * @param tgt target mz sequence
 * @return alignment indexed
 */
alingmentIndexes findThebestAlignment(mySolution *mysol, std::vector<std::string> *tgt)
{
    std::vector<int> seq1;
    std::vector<int> seq2;
    //find the best
    double localbest = DBL_MIN;
    int besttgt = 0;
    m_elemnt element;
    for(unsigned long ibest = 0; ibest <= tgt->size(); ++ibest)
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

/**
 * @brief Find a path with the best score in terms of DP, for rtCorrectedAlignment algorithm
 * @param mysol solution structure
 * @param tgt target mz sequence
 * @return alignment indexed
 */
alingmentIndexes processAlignment(mySolution *mysol, std::vector<std::string> *tgt)
{
    std::vector<int> seq1;
    std::vector<int> seq2;    
    m_elemnt element = mysol->mymatrix[mysol->src.profile_mz.size()][tgt->size()];

    int iancestor = 0;
    while(!element.ancestor.empty())
    {
        if(element.direction[iancestor] == from_left)
        {
            seq1.push_back(element.index2s1);
            seq2.push_back(-1);
        }
        else if(element.direction[iancestor] == from_oposit)
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
            break;
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

    alingmentIndexes myAlign;
    myAlign.indexNonref = seq1;
    myAlign.indexRef = seq2;
    return myAlign;
}


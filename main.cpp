#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <chrono>

using namespace std;

int n = 0, InstanceNum = -1, QualTreshold = -1, CutRange = -1, MaxSeries = 1;
bool NonSense = false;

vector<string> SEQ;///all sequences from a file
vector<vector<int>> QUAL;///all qualities from a file
vector<vector<string>> SeqGraph;///all substrings with deletions
vector<map<int, int>> CliquePos;///seqnum: pos(counted from 1)
vector<vector<map<int, int>>> CliqueSeries;///every clique in different vector

class Node
{
public:
    int seq;
    int vert;
    vector<tuple<int, int>> ConPos;///connected verts
    Node(int s, int v, vector<tuple<int, int>> HeadL)///constructor
    {
        seq = s;
        vert = v;
        ConPos = HeadL;
    }
};
vector<Node> GraphData;

void ReadFile(int num)
{
    ifstream SeqFile, QualFile;
    string num_string = to_string(num), line;
    int LineCounter = -1;

    SeqFile.open("FASTy/" + num_string + ".fasta", ios::in | ios::binary);
    QualFile.open("QUALe/" + num_string + ".qual", ios::in | ios::binary);

    if(!SeqFile.good())///parsing .fasta
    {
        cout<<"Blad otwarcia"<<endl;
        exit(0); ///when something goes wrong
    }
    while(SeqFile.good())
    {
        getline(SeqFile, line);
        if(line[0] == '>' || line[0] == '\n');///do nothing
        else
        {
            n++;
            if(n < 7) line.pop_back();///omit newline char
            SEQ.push_back(line);
        }
    }
    if(!QualFile.good())///parsing .qual
    {
        cout<<"Blad otwarcia"<<endl;
        exit(0); ///when something goes wrong
    }
    while(QualFile.good())
    {
        QualFile>>line;///cause w want to omit whitespaces
        int k = atoi(line.c_str());
        if(line[0] == '>')
        {
            LineCounter++;
            QUAL.push_back(vector<int>());
        }
        else if(isdigit(line[0])) QUAL[LineCounter].push_back(k);
    }
    SeqFile.close();
    QualFile.close();
};

void InsertDeletion(int del)
{
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < QUAL[i].size(); j++)
        {
            if(QUAL[i][j] < del) SEQ[i][j] = '-';
        }
    }
};

void MakeSubstr(int cut)
{
    int x = 0;
    for(int i = 0; i < n; i++)
    {
        SeqGraph.push_back(vector<string>());
        string line = SEQ[i];
        x = 0;
        while((x + cut) != (line.size() + 1))
        {
            string sub = line.substr(x, cut);
            SeqGraph[i].push_back(sub);
            x++;
        }
    }
};

int HowMany(int cut)///define how many deletions are permitted
{
    switch(cut)
    {
    case 4:
    {
        int permit = 1;
        return permit;
    }
    case 5:
    {
        int permit = 2;
        return permit;
    }
    case 6:
    {
        int permit = 2;
        return permit;
    }
    case 7:
    {
        int permit = 3;
        return permit;
    }
    }
};

void MakeGraph()
{
    vector<tuple<int, int>> tmp;

    for(int i = 0; i < n; i++)
    {
        for(auto it = SeqGraph[i].begin(); it != SeqGraph[i].end(); ++it)
        {
            string vert1 = *it;
            int cnt1 = count(vert1.begin(), vert1.end(), '-');
            for(int k = 0; k < cnt1; k++) vert1.erase(vert1.begin() + vert1.find('-'));///make deletions
            for(int j = 0; j < n; j++)
            {
                for(auto it2 = SeqGraph[j].begin(); it2 != SeqGraph[j].end(); ++it2)
                {
                    string vert2 = *it2;
                    int cnt2 = count(vert2.begin(), vert2.end(), '-');
                    for(int k = 0; k < cnt2; k++) vert2.erase(vert2.begin() + vert2.find('-'));///make deletions
                    if(cnt1 <= HowMany(CutRange) && cnt2 <= HowMany(CutRange) && i != j)///we concern only those with proper howmany
                    {
                        int checker = 0;
                        if(vert1.size() == vert2.size())///same size, possible identical
                        {
                            for(int k = 0; k < vert1.size(); k++)
                            {
                                if(vert1[k] == vert2[k]) checker++;
                            }
                            if(checker == vert1.size() || checker == vert1.size() - 1)
                            {
                                auto pos = make_tuple(j, distance(SeqGraph[j].begin(), it2));
                                tmp.push_back(pos);
                            }
                        }
                        else if(vert1.find(vert2) != string::npos)///contains
                        {
                            auto pos = make_tuple(j, distance(SeqGraph[j].begin(), it2));
                            tmp.push_back(pos);
                        }
                        else if(vert2.find(vert1) != string::npos)///contains
                        {
                            auto pos = make_tuple(j, distance(SeqGraph[j].begin(), it2));
                            tmp.push_back(pos);
                        }
                    }
                }
            }
            Node *v = new Node(i, distance(SeqGraph[i].begin(), it), tmp);
            GraphData.push_back(*v);
            tmp.clear();
        }
    }
/*
    cout<<endl<<GraphData.size()<<endl;
    for(auto it = GraphData.begin(); it != GraphData.end(); ++it)
    {
        auto pos = &(*it);
        cout<<pos -> seq<<"\t";
    }*/
};

void FindClique(int admiss)
{
    int CliqueSize = 1;///obvious
    map<int, int> temp;
    string TmpSeq(CutRange, 'N');///to check if a clique makes a seq
///MAIN LOOP///
    for(int i = 0; i < (n - 1); i++)
    {
        for(auto it = SeqGraph[i].begin(); it != SeqGraph[i].end(); ++it)///parse upper sequence
        {
            CliqueSize = 1;
            temp.clear();
            temp[i + 1] = distance(SeqGraph[i].begin(), it);
            string vert1 = *it;
            int cnt1 = count(vert1.begin(), vert1.end(), '-');///count deletions in both vertices
            for(int j = 1; j < n; j++)
            {
                for(auto it2 = SeqGraph[j].begin(); it2 != SeqGraph[j].end(); ++it2)///parse lower sequence
                {
                    string vert2 = *it2;
                    int cnt2 = count(vert2.begin(), vert2.end(), '-');///same here
                    if(cnt1 <= admiss && cnt2 <= admiss)
                    {
                        int checker = 0;
                        TmpSeq.clear();
                        for(int i = 0; i < vert1.size(); i++)
                        {
                            if(vert1[i] == vert2[i] || vert1[i] == '-' || vert2[i] == '-') checker++;///check similarity, if it has value of substring, its ok
                            if(vert1[i] != '-') TmpSeq.push_back(vert1[i]);
                            else if(vert2[i] != '-') TmpSeq.push_back(vert2[i]);///thats why we have small things in small quals
                        }
                        if(checker == CutRange && i < j && TmpSeq.size() == CutRange)///despite all it may happen, so we dont want clique in range of the same seq!
                        {
                            CliqueSize++;///it is obvious
                            temp[j + 1] = distance(SeqGraph[j].begin(), it2);
                            break;///exit sequence, start another
                        }
                    }
                }
            }
            if(CliqueSize >= 4) CliquePos.push_back(temp);
        }
    }
///END///
};

int FindSeries()
{
    vector<int> temp, temp2;
    map<int, int> tmp, tmp2;
    int checker = 0, SeriesSize = 1, SeriesNum = -1;
    bool HasSame = true;
///BEGIN///
    for(int i = 0; i < (CliquePos.size() - 1); i++)
    {
        tmp.clear();
        for(auto it = CliquePos[i].begin(); it != CliquePos[i].end(); ++it) tmp[it->first] = it->second;///just push back a referential
        if(SeriesSize <= 1 && SeriesNum >= 0)
        {
            SeriesNum--;
            CliqueSeries.pop_back();
        }
        if(SeriesSize > MaxSeries) MaxSeries = SeriesSize;
        else;///do nothing
        SeriesSize = 1;
        SeriesNum++;
        CliqueSeries.push_back(vector<map<int, int>>());
        CliqueSeries[SeriesNum].push_back(tmp);///reference
        for(int j = 1; j < CliquePos.size(); j++)
        {
            tmp2.clear();
            HasSame = true;
            for(auto it2 = CliquePos[j].begin(); it2 != CliquePos[j].end(); ++it2) tmp2[it2->first] = it2->second;///comparable
            map<int, int>::iterator it3, it4;///check if seq are same
            for(it3 = tmp.begin(), it4 = tmp2.begin(); it3 != tmp.end(); ++it3, ++it4)
            {
                if(it3->first == it4->first);///pass
                else HasSame = false;
            }
            if(tmp.size() == tmp2.size() && HasSame && i < j)
            {
                checker = 0;
                temp.clear();
                temp2.clear();
                for(auto it = tmp.begin(); it != tmp.end(); ++it) temp.push_back(it->second);
                for(auto it = tmp2.begin(); it != tmp2.end(); ++it) temp2.push_back(it->second);
                for(int k = 0; k < temp.size(); k++)
                {
                    if(temp[k] == temp2[k] - 1) checker++;///check if it is a next clique
                }
                if(checker >= temp.size() - HowMany(temp.size()))
                {
                    while(tmp2 != CliqueSeries[SeriesNum].back())
                    {
                        SeriesSize++;
                        CliqueSeries[SeriesNum].push_back(tmp2);
                        tmp = tmp2;
                    }
                }
            }
        }
    }
    if(MaxSeries == 1) NonSense = true;///it may happen that there is no sense of searching
};

int ChooseOne()
{
    int chosen;
    vector<int> ToChoose;///which is best, if more than 1 randomize a choice

    for(int i = 0; i < CliqueSeries.size(); i++)
    {
        if(CliqueSeries[i].size() == MaxSeries) ToChoose.push_back(i);
    }

    unsigned seed = chrono::system_clock::now().time_since_epoch().count();///randomizing
    shuffle(ToChoose.begin(), ToChoose.end(), default_random_engine(seed));

    chosen = ToChoose[0];///thats the one we choose
    return chosen;
};

void ShowAlign(int chosen)
{
    string Answer, tmp2(CutRange, 'X'), tmp3;
    vector<string> tmp;

    cout<<endl<<"Wybrano serie o dlugosci "<<MaxSeries;
    cout<<endl<<"Ogolna dlugosc dopasowania w takim przypadku wynosi "<<(CutRange + MaxSeries) - 1;

    for(int i = 0; i < CliqueSeries[chosen].size(); i++)///print a clique series that whas choosen
    {
        cout<<endl;
        for(auto it = CliqueSeries[chosen][i].begin(); it != CliqueSeries[chosen][i].end(); ++it) cout<<it->first<<": "<<it->second<<"\t";
    }

    for(int i = 0; i < CliqueSeries[chosen].size(); i++)///iterate on cliques
    {
        cout<<endl<<"Klika "<<i + 1<<":";
        for(auto it = CliqueSeries[chosen][i].begin(); it != CliqueSeries[chosen][i].end(); ++it)
        {
            tmp.push_back(SeqGraph[it->first - 1][it->second]);
            cout<<endl<<SeqGraph[it->first - 1][it->second];
        }
        for(int k = 0; k < tmp[0].size(); k++)
        {
            for(int j = 0; j < tmp.size(); j++)
            {
                if(tmp[j][k] != '-') tmp3.push_back(tmp[j][k]);
            }
            for(int l = 0; l < tmp3.size(); l++)
            {
                float cnt = count(tmp3.begin(), tmp3.end(), tmp3[l]);
                float cap = tmp3.size();
                if(cnt >= ceil(cap / 2))
                {
                    tmp2[k] = tmp3[l];///choose one of majority
                    break;
                }
                else tmp2[k] = tmp3[(rand() % tmp3.size()) + 0];///if all are once, randomize
            }
            tmp3.clear();
        }
        cout<<endl<<"Konsensus: "<<tmp2;
        if(Answer.size() == 0) Answer = tmp2;
        else Answer += tmp2.back();///add just one, its obvious
        tmp3.clear();
        tmp.clear();
    }
    cout<<endl<<endl<<"Wyjsciowa sekwencja: "<<Answer;
};

int main()
{

    while(InstanceNum > 5 || InstanceNum < 1)
    {
        cout<<"Prosze podac numer instancji (1-5): ";
        cin>>InstanceNum;
    }
    while(CutRange > 7 || CutRange < 4)
    {
        cout<<"Prosze podac dlugosc podciagu (4-7): ";
        cin>>CutRange;
    }
    while(QualTreshold > 30 || QualTreshold < 0)
    {
        cout<<"Prosze podac ograniczenie jakosci (0-30): ";
        cin>>QualTreshold;
    }

    ReadFile(InstanceNum);
    InsertDeletion(QualTreshold);
    MakeSubstr(CutRange);
    MakeGraph();
    FindClique(HowMany(CutRange));///function in function,cause we want to know how many deletions are permitted

    if(NonSense == false)
    {
        FindSeries();
        ShowAlign(ChooseOne());
        return 0;
    }
    else
    {
        cout<<endl<<"Nie znaleziono zadowalajacego rozwiazania. Koncze prace...";
        return 1;
    }
};

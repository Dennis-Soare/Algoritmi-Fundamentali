#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <stack>
#include <bits/stdc++.h>

using namespace std;

const int Nmax=100005;

ifstream in("sortaret.in");
ofstream out("sortaret.out");

class Graf {
    int nrNoduri, nrMuchii;
    int Start;
    vector<vector<int>> vecini;
    int culoare[Nmax];          // pentru problema DFS
    stack<int> stiva;
    vector<int> lowlink;        // pentru problema CTC
    vector<bool> pe_stiva, vizitat;
    int nrCTC = 0;
    vector<int> grad_intern;

public:

    void problema_BFS();
    void DFS(int nod, int color);
    void problema_DFS();
    void problema_CTC();
    void tarjan(int nod);
    void Sortare_topologica();
    void afisare_vecini();
    void citire_graf();
    bool Havel_Hakimi();
};

bool Graf::Havel_Hakimi()
{
    in>>nrNoduri;
    int grad[nrNoduri+1], suma=0;
    for(int i=0; i<nrNoduri; i++)
    {
        in >> grad[i];
        suma += grad[i];
        if(grad[i] > (nrNoduri - 1))
        {
            out<<"NU";
            return false;
        }
    }
    if(suma % 2 == 1)
    {
        out <<"NU";
        return false;
    }
    sort(grad, grad+nrNoduri);
    while(grad[0] != 0)
    {
        for(int i=1; i<=grad[0]; i++)
        {
            grad[i]--;
            if(grad[i] < 0)
            {
                out<<"NU";
                return false;
            }
        }
        grad[0] = 0;
        sort(grad, grad+nrNoduri);
    }
    out<<"DA";
    return true;
}

void Graf::citire_graf()
{
    in >> nrNoduri >> nrMuchii;
    vecini.resize(nrNoduri+1);
    grad_intern.resize(nrNoduri+1);

    for(int i=0; i<=nrNoduri; i++)
    {
        grad_intern[i] = 0;
    }

    for(int i=0; i<nrMuchii; i++)
    {
        int x,y;
        in>>x>>y;
        vecini[x].push_back(y);
        grad_intern[y]++;
    }
}

void Graf::afisare_vecini()
{
    for(int i=1; i<=nrNoduri; i++)
    {
        cout<< i <<": ";
        int l = vecini[i].size();
        for(int j=0; j<l; j++)
        {
            cout<< vecini[i][j] <<' ';
        }
        cout<<"\n";
    }
}

void Graf::Sortare_topologica()
{
    citire_graf();

    int contor = 0;
    while(contor < nrNoduri)
    {
        for(int nod = 1; nod <= nrNoduri; nod++)
        {
            if(grad_intern[nod] == 0)
            {
                contor++;
                out << nod <<' ';
                int l = vecini[nod].size();
                for(int i=0; i<l; i++)
                {
                    grad_intern[vecini[nod][i]]--;
                }
                grad_intern[nod]--;
            }
        }
    }
}

/*void Graf::tarjan(int nod)
{
    int l = vecini[nod].size();
    for(int i=0; i<l; i++)
    {
        int vec = vecini[nod][i];
        if(!vizitat[vec])
        {
            stiva.push(vec);
            pe_stiva[vec] = true;
            vizitat[vec] = true;
            tarjan(vec);
            if(pe_stiva[nod])
                lowlink[nod] = min(lowlink[nod], lowlink[vec]);
        }
    }
    if(lowlink[nod] = nod)
    {
        nrCTC++;
        while(lowlink[stiva.top()] == nod)
        {
            pe_stiva[stiva.top()] = false;
            stiva.pop();
        }
    }
}

void Graf::problema_CTC()
{
    in >> nrNoduri >> nrMuchii;
    vecini.resize(nrNoduri+1);
    lowlink.resize(nrNoduri+1);
    pe_stiva.resize(nrNoduri+1);
    vizitat.resize(nrNoduri+1);
    nrCTC = 0;

    for(int i=0; i<nrMuchii; i++)
    {
        int x,y;
        in>>x>>y;
        vecini[x].push_back(y);
    }

    for(int i=1; i<=nrNoduri+1; i++)
    {
        vizitat[i] = false;
        pe_stiva[i] = false;
        lowlink[i] = i;
    }
    for(int nod=1; nod<=nrNoduri; nod++)
    {
        if(!vizitat[nod])
        {
            stiva.push(nod);
            pe_stiva[nod] = true;
            vizitat[nod] = true;
            tarjan(nod);
        }
    }
    cout << nrCTC;
}*/

void Graf::DFS(int nod, int color)
{
    culoare[nod] = color;
    int l = vecini[nod].size();
    for(int i=0; i<l; i++)
    {
        int vec = vecini[nod][i];
        if(culoare[vec] == 0)
        {
            DFS(vec, color);
        }
    }
}

void Graf::problema_DFS()
{
    in>>nrNoduri>>nrMuchii;
    vecini.resize(nrNoduri+1);

    for(int i=0; i<nrMuchii; i++)
    {
        int x,y;
        in>>x>>y;
        vecini[x].push_back(y);
        vecini[y].push_back(x);
    }

    int nrCulori = 0;
    for(int nod=1; nod<=nrNoduri; nod++)
    {
        if(culoare[nod] == 0)
        {
            nrCulori++;
            DFS(nod, nrCulori);
        }
    }
    out << nrCulori;
}

void Graf::problema_BFS()
{
    in>>nrNoduri>>nrMuchii>>Start;

    vecini.resize(nrNoduri+1);
    queue<int> coada;
    bool vizitat[nrNoduri+1];
    int inaltime[nrNoduri+1];

    for(int i=1; i<=nrNoduri+1; i++)
    {
        vizitat[i] = false;
        inaltime[i] = 0;
    }
    for(int i=0; i<nrMuchii; i++)
    {
        int x,y;
        in>>x>>y;
        vecini[x].push_back(y);
    }
    coada.push(Start);
    vizitat[Start] = true;
    while(!coada.empty())
    {
        int nod_curent = coada.front();
        coada.pop();
        int l = vecini[nod_curent].size();
        for(int i=0; i<l; i++)
        {
            int x = vecini[nod_curent][i];
            if(!vizitat[x])
            {
                vizitat[x] = true;
                coada.push(x);
                inaltime[x] = inaltime[nod_curent] + 1;
            }
        }
    }
    for(int i=1; i<=nrNoduri; i++)
    {
        if(inaltime[i] == 0 && i != Start)
            out << -1 <<' ';
        else
            out<<inaltime[i]<<' ';
    }
}

int main()
{
    Graf G;
    in.close();
    out.close();
    return 0;
}

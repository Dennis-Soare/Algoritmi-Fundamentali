#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <stack>
#include <bits/stdc++.h>

using namespace std;

const int NMAX = numeric_limits<int>::max();

class Graf {
    int nrNoduri, nrMuchii;
    vector<vector<int>> vecini;

    void afis_lista_ad();
    vector<int> BFS(const int Start);
    vector<int> DFS(const int Start);
    vector<int> distante_BFS(const int Start);
    void parcurgere_DFS(int nod_curent, vector<int>& sol, vector<bool>& viz);
    vector<int> sortare_topologica(vector<int> grad_intern);
    void reuneste(int x, int y, vector<int>& tata, vector<int>& inaltime);
    int reprez(int x, vector<int>& tata);
    void Initializare(int u,vector<int>& r);
    int Reprez_Kruskal(int u, vector<int>& r);
    void Reuneste_Kruskal(int u, int v, vector<int>&r);

public:
    void problema_BFS();
    void problema_DFS();
    void problema_sortaret();
    void Havel_Hakimi();
    void problema_APM();
    void Paduri_de_multimi_disjuncte();
    void problema_Dijkstra();
    void problema_Bellman_Ford();
};

void Graf::afis_lista_ad()
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

vector<int> Graf::distante_BFS(const int Start)
{
    queue<int> coada;
    vector<int> inaltime;
    vector<bool> vizitat;
    inaltime.resize(nrNoduri+1);
    vizitat.resize(nrNoduri+1);
    for(int i=0; i<=nrNoduri; i++) {
        inaltime[i] = 0;
        vizitat[i] = false;
    }

    coada.push(Start);
    vizitat[Start] = true;

    while(!coada.empty()) {
        int nod_curent = coada.front();
        coada.pop();
        for(int i=0; i< int(vecini[nod_curent].size()); i++) {
            int x = vecini[nod_curent][i];
            if(!vizitat[x]) {
                vizitat[x] = true;
                coada.push(x);
                inaltime[x] = inaltime[nod_curent] + 1;
            }
        }
    }
    for(int i=1; i<=nrNoduri; i++) {
        if(inaltime[i] == 0 && i != Start)
            inaltime[i] = -1;
    }
    return inaltime;
}

vector<int> Graf::BFS(const int Start)
{
    vector<int>sol;
    queue<int> coada;
    vector<bool> vizitat;
    vizitat.resize(nrNoduri+1);
    for(int i=0; i<=nrNoduri; i++) {
        vizitat[i] = false;
    }

    coada.push(Start);
    vizitat[Start] = true;

    while(!coada.empty()) {
        int nod_curent = coada.front();
        coada.pop();
        sol.push_back(nod_curent);
        for(int i=0; i< int(vecini[nod_curent].size()); i++) {
            int x = vecini[nod_curent][i];
            if(!vizitat[x]) {
                vizitat[x] = true;
                coada.push(x);
            }
        }
    }
    return sol;
}

void Graf::problema_BFS()
{
    ifstream in("bfs.in");
    ofstream out("bfs.out");

    int Start;
    in>>nrNoduri>>nrMuchii>>Start;

    vecini.resize(nrNoduri+1);
    for(int i=0; i<nrMuchii; i++)
    {
        int x,y;
        in>>x>>y;
        vecini[x].push_back(y);
    }
    vector<int> distante = distante_BFS(Start);
    for(int i=1; i<=nrNoduri; i++) {
        out << distante[i]<< ' ';
    }

    in.close();
    out.close();
}

void Graf::parcurgere_DFS(int nod_curent, vector<int>& sol, vector<bool>& viz)
{
    viz[nod_curent] = true;
    sol.push_back(nod_curent);
    for(int i=0; i<int(vecini[nod_curent].size()); i++) {
        int nod = vecini[nod_curent][i];
        if(!viz[nod]) {
            parcurgere_DFS(nod, sol , viz);
        }
    }
}

vector<int> Graf::DFS(const int Start)
{
    vector<int> sol;
    vector<bool> viz;
    viz.resize(nrNoduri+1);
    for(int i=0; i<=nrNoduri; i++) {
        viz[i] = false;
    }
    int nod_curent = Start;
    parcurgere_DFS(nod_curent, sol, viz);
    return sol;
}

void Graf::problema_DFS()
{
    ifstream in("dfs.in");
    ofstream out("dfs.out");

    in>>nrNoduri>>nrMuchii;
    vecini.resize(nrNoduri+1);

    for(int i=0; i<nrMuchii; i++){
        int x,y;
        in>>x>>y;
        vecini[x].push_back(y);
        vecini[y].push_back(x);
    }

    vector<int> nr_comp_conexa;
    nr_comp_conexa.resize(nrNoduri+1);
    for(int i=1; i<=nrNoduri; i++) {
        nr_comp_conexa[i] = 0;
    }

    int k=0;
    for(int i=1; i<=nrNoduri; i++){
        if(nr_comp_conexa[i] == 0) {
            ++k;
            vector<int> noduri = DFS(i);
            for(int j=0; j<int(noduri.size()); j++) {
                nr_comp_conexa[noduri[j]] = k;
            }
        }
    }
    out<<k;
    in.close();
    out.close();
}

vector<int> Graf::sortare_topologica(vector<int> grad_intern)
{
    vector<int> sol;
    int contor = 0;
    while(contor < nrNoduri) {
        for(int nod = 1; nod <= nrNoduri; nod++) {
            if(grad_intern[nod] == 0) {
                contor++;
                sol.push_back(nod);
                for(int i=0; i<int(vecini[nod].size()); i++){
                    grad_intern[vecini[nod][i]]--;
                }
                grad_intern[nod]--;
            }
        }
    }
    return sol;
}

void Graf::problema_sortaret()
{
    ifstream in("sortaret.in");
    ofstream out("sortaret.out");

    in>>nrNoduri>>nrMuchii;
    vecini.resize(nrNoduri+1);
    vector<int> grad_intern;
    grad_intern.resize(nrNoduri+1);
    for(int i=0; i<nrNoduri; i++){
        grad_intern[i]=0;
    }
    for(int i=0; i<nrMuchii; i++){
        int x,y;
        in>>x>>y;
        vecini[x].push_back(y);
        grad_intern[y]++;
    }

    vector<int> sol = sortare_topologica(grad_intern);
    for(int i=0; i<nrNoduri; i++){
        out<<sol[i]<<' ';
    }
    in.close();
    out.close();
}

void Graf::Havel_Hakimi()
{
    cout<<"Introduceti numarul de noduri: ";
    cin>>nrNoduri;
    int grad[nrNoduri+1], suma=0;
    cout<<"Introduceti " << nrNoduri <<" numere: ";
    for(int i=0; i<nrNoduri; i++) {
        cin >> grad[i];
        suma += grad[i];
        if(grad[i] > (nrNoduri - 1)) {
            cout<<"NU";
            return;
        }
    }
    if(suma % 2 == 1){
        cout <<"NU";
        return;
    }
    sort(grad, grad+nrNoduri, greater<int>());
    while(grad[0] != 0){
        for(int i=1; i<=grad[0]; i++){
            grad[i]--;
            if(grad[i] < 0){
                cout<<"NU";
                return;
            }
        }
        grad[0] = 0;
        sort(grad, grad+nrNoduri, greater<int>());
    }
    cout<<"DA";
}

void Graf::reuneste(int x, int y, vector<int>& tata, vector<int>& inaltime)
{
    int rx = reprez(x, tata), ry = reprez(y, tata);
    if(inaltime[rx] > inaltime[ry]) {
        tata[ry] = rx;
    }
    else {
        tata[rx] = ry;
        if(inaltime[rx] == inaltime[ry])
            inaltime[ry]++;
    }
}

int Graf::reprez(int x, vector<int>& tata)
{
    if(tata[x] == 0)
        return x;
    tata[x] = reprez(tata[x], tata);
    return tata[x];
}

void Graf::Paduri_de_multimi_disjuncte()
{
    ifstream in("disjoint.in");
    ofstream out("disjoint.out");

    int nrOperatii;
    vector<int> tata, inaltime;
    in>>nrNoduri>>nrOperatii;
    tata.resize(nrNoduri+1);
    inaltime.resize(nrNoduri+1);
    for(int i=0; i<=nrNoduri;i++){
        tata[i] = inaltime[i] = 0;
    }

    for(int i=0;i<nrOperatii;i++) {
        int operatie,x,y;
        in >> operatie >> x >> y;
        if(operatie == 1) {
            reuneste(x, y, tata, inaltime);
        }
        else if(operatie == 2) {
            if(reprez(x, tata) == reprez(y, tata))
                out <<"DA\n";
            else out<<"NU\n";
        }
    }

    in.close();
    out.close();
}

void Graf::Initializare(int u, vector<int>& r) {
    r[u] = u;
}

int Graf::Reprez_Kruskal(int u, vector<int>&r) {
    return r[u];
}

void Graf::Reuneste_Kruskal(int u, int v, vector<int>&r) {
    int r1 = Reprez_Kruskal(u, r);
    int r2 = Reprez_Kruskal(v, r);
    for(int k=1; k<=nrNoduri; k++) {
        if(r[k] == r2)
            r[k] = r1;
    }
}

void Graf::problema_APM()
{
    ifstream in("apm.in");
    ofstream out("apm.out");


    in>>nrNoduri>>nrMuchii;
    vecini.resize(nrNoduri+1);
    vector<int> r;
    r.resize(nrNoduri);
    vector<vector<pair<int, int>>>drumuri;
    drumuri.resize(nrNoduri+1);

    for(int i=0; i<nrMuchii; i++){
        int a,b,c;
        in>>a>>b>>c;
        drumuri[a].push_back(make_pair(b,c));
    }
    for(int i=1;i<=nrNoduri;i++){
        cout<<i<<": ";
        for(int j=0; j<int(drumuri[i].size());j++){
            cout<<drumuri[i][j].first<<' '<<drumuri[i][j].second<<"; ";
        }
        cout<<endl;
    }
    //sorteaza muchiile
    for(int i=1; i<=nrNoduri; i++){
        Initializare(i, r);
    }

    in.close();
    out.close();
}

void Graf::problema_Dijkstra()
{
    ifstream in("dijkstra.in");
    ofstream out("dijkstra.out");

    in>>nrNoduri>>nrMuchii;
    vector<vector<pair<int, int>>>drumuri;
    vector<int> dist, tata;
    vector<bool> vizitat;
    drumuri.resize(nrNoduri+1);
    dist.resize(nrNoduri+1);
    tata.resize(nrNoduri+1);
    vizitat.resize(nrNoduri+1);
    for(int i=0; i<=nrNoduri;i++) {
        tata[i] = 0;
        dist[i] = NMAX;
        vizitat[i] = false;
    }

    for(int i=0; i<nrMuchii; i++) {
        int a,b,c;
        in>>a>>b>>c;
        drumuri[a].push_back(make_pair(b,c));
    }
    /*for(int i=1;i<=nrNoduri;i++){
        cout<<i<<": ";
        for(int j=0; j<int(drumuri[i].size());j++){
            cout<<drumuri[i][j].first<<' '<<drumuri[i][j].second<<"; ";
        }
        cout<<endl;
    }*/

    bool ok = false;
    int nod_curent = 1;
    dist[1] = 0;
    while(!ok) {
        ok = true;
        vizitat[nod_curent] = true;
        for(int i=0; i<int(drumuri[nod_curent].size());i++){
            int vecin = drumuri[nod_curent][i].first;
            int cost = drumuri[nod_curent][i].second;
            if(dist[nod_curent] + cost < dist[vecin]){
                dist[vecin] = dist[nod_curent] + cost;
                tata[vecin] = nod_curent;
            }
        }
        int et_minima = NMAX;
        for(int i=1; i<=nrNoduri; i++){
            if(!vizitat[i] && dist[i] < et_minima){
                nod_curent = i;
                ok = false;
            }
        }
    }
    for(int i=2; i<=nrNoduri; i++) {
        out<<dist[i]<<' ';
    }

    in.close();
    out.close();
}

void Graf::problema_Bellman_Ford()
{
    ifstream in("bellmanford.in");
    ofstream out("bellmanford.out");

    in>>nrNoduri>>nrMuchii;
    vector<vector<pair<int, int>>>drumuri;
    vector<int> dist, tata;
    drumuri.resize(nrNoduri+1);
    dist.resize(nrNoduri+1);
    tata.resize(nrNoduri+1);
    for(int i=0; i<=nrNoduri;i++) {
        tata[i] = 0;
        dist[i] = NMAX;
    }

    for(int i=0; i<nrMuchii; i++) {
        int a,b,c;
        in>>a>>b>>c;
        drumuri[a].push_back(make_pair(b,c));
    }
    /*for(int i=1;i<=nrNoduri;i++){
        cout<<i<<": ";
        for(int j=0; j<int(drumuri[i].size());j++){
            cout<<drumuri[i][j].first<<' '<<drumuri[i][j].second<<"; ";
        }
        cout<<endl;
    }*/

    queue<int> coada;
    dist[1] = 0;
    coada.push(1);
    while(!coada.empty()) {
        int nod_curent = coada.front();
        coada.pop();
        for(int i=0; i<int(drumuri[nod_curent].size());i++) {
            int vecin = drumuri[nod_curent][i].first;
            int cost = drumuri[nod_curent][i].second;
            if(dist[nod_curent] + cost < dist[vecin]){
                dist[vecin] = dist[nod_curent] + cost;
                tata[vecin] = nod_curent;
                coada.push(vecin);
            }
        }
    }
    for(int i=2; i<=nrNoduri; i++) {
        out<<dist[i]<<' ';
    }
    in.close();
    out.close();
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

int main()
{
    Graf G;
    //G.problema_BFS();
    //G.problema_DFS();
    //G.problema_sortaret();
    //G.Havel_Hakimi();
    //G.problema_APM();
    //G.Paduri_de_multimi_disjuncte();
    //G.problema_Dijkstra();
    //G.problema_Bellman_Ford();
    return 0;
}

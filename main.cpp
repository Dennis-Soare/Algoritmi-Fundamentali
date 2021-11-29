#include <bits/stdc++.h>

using namespace std;

const int NMAX = numeric_limits<int>::max();

class Graf {
    int nrNoduri, nrMuchii;
    vector<vector<int>> vecini;

    int diametru_arbore();
    vector<vector<int>> Roy_Floyd(vector<vector<int>> matrice);
    vector<int> BFS(const int Start);
    vector<int> DFS(const int Start);
    vector<int> distante_BFS(const int Start);      // calculeaza distanta de la start la fiecare nod
    void parcurgere_DFS(const int nod_curent, vector<int>& sol, vector<bool>& viz);     //functie ajutatoare pentru DFS
    vector<int> sortare_topologica(vector<int> grad_intern);
    vector<vector<int>> ctc();
    void tarjan(int nod_curent, vector<bool>& viz, vector<bool>& pe_stiva, vector<int>& low, stack<int>& stiva, vector<vector<int>>& sol);
    vector<int> Dijkstra(const vector<vector<pair<int, int>>> &lista_ad, const int sursa);
    vector<pair<int,int>> Prim(const vector<vector<pair<int, int>>> &lista_ad, int& cost_arbore, const int sursa);
    vector<int> Bellman_Ford(const vector<vector<pair<int, int>>>&, const int sursa, bool& ciclu_negativ);
    void reuneste(int x, int y, vector<int>& tata, vector<int>& inaltime);  // pentru Paduri_de_multimi_disjuncte
    int reprez(int x, vector<int>& tata);                                   // pentru Paduri_de_multimi_disjuncte
    void afis_lista_ad();

public:
    void problema_BFS();
    void problema_DFS();
    void problema_sortaret();
    void Havel_Hakimi();
    void problema_APM();
    void Paduri_de_multimi_disjuncte();
    void problema_Dijkstra();
    void problema_Bellman_Ford();
    void problema_Roy_Floyd();
    void problema_darb();
    void problema_ctc();
};


int Graf::diametru_arbore()
{
    vector<int> bfs = BFS(1);
    vector<int> dist = distante_BFS(bfs[nrNoduri-1]);       //calculam distantele nodurilor de la ultimul nod gasit in bfs

    int sol = 0;
    for(int i=1; i<=nrNoduri; i++){         //determinam distanta dintre ultimul nod din bfs si cel mai indepartat nod de el
        if(dist[i]+1 > sol)
            sol = dist[i]+1;
    }

    return sol;
}

void Graf::problema_darb()
{
    ifstream in("darb.in");
    ofstream out("darb.out");

    in>>nrNoduri;
    vecini.resize(nrNoduri+1);
    for(int i=1;i<nrNoduri;i++){
        int x,y;
        in>>x>>y;
        vecini[x].push_back(y);
        vecini[y].push_back(x);
    }

    int sol = diametru_arbore();
    out<<sol;

    in.close();
    out.close();
}

void Graf::problema_Roy_Floyd()
{
    ifstream in("royfloyd.in");
    ofstream out("royfloyd.out");

    in>>nrNoduri;
    vector<vector<int>> matrice(nrNoduri);

    for(int i=0; i<nrNoduri; i++){
        for(int j=0; j<nrNoduri; j++){
            int x;
            in>>x;
            matrice[i].push_back(x);
        }
    }

    vector<vector<int>> sol = Roy_Floyd(matrice);

    for(int i=0; i<nrNoduri; i++){
        for(int j=0; j<nrNoduri; j++){
                out<< sol[i][j]<<' ';
        }
        out<<'\n';
    }

    in.close();
    out.close();
}

vector<vector<int>> Graf::Roy_Floyd(vector<vector<int>> matrice)
{
    for(int i=0; i<nrNoduri; i++){
        for(int j=0; j<nrNoduri; j++){
                if(matrice[i][j] == 0 && i!=j)          //daca nu exista muchie intre i si j
                    matrice[i][j] = 1e9;
        }
    }

    for(int k=0; k<nrNoduri; k++){
        for(int i=0; i<nrNoduri; i++){
            for(int j=0; j<nrNoduri; j++){
                if (matrice[i][k] + matrice[k][j] < matrice[i][j])
                        matrice[i][j] = matrice[i][k] + matrice[k][j];
            }
        }
    }

    for(int i=0; i<nrNoduri; i++){
        for(int j=0; j<nrNoduri; j++){
            if(matrice[i][j] == 1e9)
                matrice[i][j] = 0;
        }
    }
    return matrice;
}

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
    vector<int> dist(nrNoduri+1, 0);
    vector<bool> vizitat(nrNoduri+1, false);

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
                dist[x] = dist[nod_curent] + 1;
            }
        }
    }
    for(int i=1; i<=nrNoduri; i++) {
        if(dist[i] == 0 && i != Start)
            dist[i] = -1;
    }
    return dist;
}

vector<int> Graf::BFS(const int Start)
{
    vector<int>sol;
    queue<int> coada;
    vector<bool> vizitat(nrNoduri+1, false);

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

void Graf::parcurgere_DFS(const int nod_curent, vector<int>& sol, vector<bool>& viz)
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
    vector<bool> viz(nrNoduri+1, false);
    parcurgere_DFS(Start, sol, viz);
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

    vector<int> nr_comp_conexa(nrNoduri+1, 0);

    int sol = 0;

    for(int i=1; i<=nrNoduri; i++){
        if(nr_comp_conexa[i] == 0){
            sol++;
            vector<int> dfs = DFS(i);
            for(int i=0; i<int(dfs.size()); i++){
                nr_comp_conexa[dfs[i]] = sol;
            }
        }
    }
    out<<sol;

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

vector<pair<int,int>> Graf::Prim(const vector<vector<pair<int, int>>> &lista_ad, int& cost_arbore, const int sursa)
{
    vector<bool> viz(nrNoduri+1, false);
    vector<int> tata(nrNoduri+1, -1), cost_min(nrNoduri+1, NMAX);
    priority_queue<pair<int,int>, vector<pair<int,int>>, greater<pair<int,int>>> q;

    q.push({0, sursa});
    cost_min[sursa] = 0;

    while(!q.empty()){
        int nod = q.top().second;
        int cost_muchie = q.top().first;
        q.pop();

        if(!viz[nod]){

            viz[nod] = true;
            cost_arbore += cost_muchie;

            for(auto vecin: lista_ad[nod]){
                int nod_vecin = vecin.first;
                int cost_vecin = vecin.second;

                if(viz[nod_vecin] == false && cost_min[nod_vecin] > cost_vecin){
                    cost_min[nod_vecin] = cost_vecin;
                    q.push({cost_vecin, nod_vecin});
                    tata[nod_vecin] = nod;
                }
            }
        }

    }
    vector<pair<int,int>> sol;
    for(int i=1; i<=nrNoduri; i++){
        sol.push_back({tata[i],i});
    }
    return sol;

}

void Graf::problema_APM()
{
    ifstream in("apm.in");
    ofstream out("apm.out");


    in>>nrNoduri>>nrMuchii;
    vector<vector<pair<int, int>>> lista_ad(nrNoduri+1);

    for(int i=1; i<=nrMuchii; i++){
        int x,y,cost;
        in>>x>>y>>cost;
        lista_ad[x].push_back({y, cost});
        lista_ad[y].push_back({x, cost});
    }

    /*for(int i=1; i<=nrNoduri;i++){
        cout<<i<<": ";
        for(int j=0; j<int(lista_ad[i].size());j++){
            cout<<lista_ad[i][j].first <<' ' <<lista_ad[i][j].second<<", ";
        }
        cout<<endl;
    }*/

    int cost_arbore = 0;
    vector<pair<int,int>> sol = Prim(lista_ad, cost_arbore, 1);

    out<<cost_arbore<<'\n'<<nrNoduri-1<<'\n';

    for(int i=1; i<=nrNoduri-1; i++){
        out<<sol[i].first<<' '<<sol[i].second<<'\n';
    }


    in.close();
    out.close();
}

vector<int> Graf::Dijkstra(const vector<vector<pair<int, int>>>& lista_ad, const int sursa)
{
    vector<int> dist(nrNoduri+1, NMAX);
    vector<bool> viz(nrNoduri+1, false);
    dist[sursa] = 0;
    priority_queue<pair<int,int>, vector<pair<int,int>>, greater<pair<int,int>>> q;
    q.push({0, sursa});

    while(!q.empty()){
        int nod = q.top().second;
        q.pop();

        if(!viz[nod]){
            viz[nod] = true;
            for(auto vecin: lista_ad[nod]){
                int nod_vecin = vecin.first;
                int cost_vecin = vecin.second;
                if(dist[nod] + cost_vecin < dist[nod_vecin]){
                    dist[nod_vecin] = dist[nod] + cost_vecin;
                    q.push({dist[nod_vecin], nod_vecin});
                }
            }
        }
    }
    return dist;
}

void Graf::problema_Dijkstra()
{
    ifstream in("dijkstra.in");
    ofstream out("dijkstra.out");

    in>>nrNoduri>>nrMuchii;
    vector<vector<pair<int, int>>> lista_ad(nrNoduri+1);

    for(int i=1; i<=nrMuchii; i++){
        int x,y,cost;
        in>>x>>y>>cost;
        lista_ad[x].push_back(make_pair(y, cost));
    }

    /*for(int i=1; i<=nrNoduri;i++){
        cout<<i<<": ";
        for(int j=0; j<int(lista_ad[i].size());j++){
            cout<<lista_ad[i][j].first <<' ' <<lista_ad[i][j].second<<", ";
        }
        cout<<endl;
    }*/
    vector<int> sol = Dijkstra(lista_ad, 1);
    for(int i=2; i<=nrNoduri; i++){
        if(sol[i]==NMAX)
            out<<"0 ";
        else
            out<<sol[i]<<' ';
    }

    in.close();
    out.close();
}

vector<int> Graf::Bellman_Ford(const vector<vector<pair<int, int>>> &lista_ad, const int sursa, bool& ciclu_negativ)
{
    vector<int> dist(nrNoduri+1, NMAX), viz(nrNoduri+1, 0);
    dist[sursa] = 0;
    ciclu_negativ = false;

    priority_queue<pair<int,int>, vector<pair<int,int>>, greater<pair<int,int>>> q;
    q.push({0, sursa});

    while(!q.empty()){
        int nod = q.top().second;
        q.pop();
        viz[nod]++;

        if (viz[nod] == nrNoduri){
            dist.resize(0);
            ciclu_negativ = true;
            break;
        }

        for(auto vecin: lista_ad[nod]){
            int nod_vecin = vecin.first;
            int cost = vecin.second;

            if(dist[nod] + cost < dist[nod_vecin]){         //relaxam muchia
                dist[nod_vecin] = dist[nod] + cost;
                q.push({dist[nod_vecin], nod_vecin});
            }
        }

    }
    return dist;
}

void Graf::problema_Bellman_Ford()
{
    ifstream in("bellmanford.in");
    ofstream out("bellmanford.out");

    in>>nrNoduri>>nrMuchii;
    vector<vector<pair<int, int>>> lista_ad(nrNoduri+1);

    for(int i=1; i<=nrMuchii; i++){
        int x,y,cost;
        in>>x>>y>>cost;
        lista_ad[x].push_back({y, cost});
    }

    /*for(int i=1; i<=nrNoduri;i++){
        cout<<i<<": ";
        for(int j=0; j<int(lista_ad[i].size());j++){
            cout<<lista_ad[i][j].first <<' ' <<lista_ad[i][j].second<<"; ";
        }
        cout<<endl;
    }*/

    bool ciclu_negativ;
    vector<int> sol = Bellman_Ford(lista_ad, 1, ciclu_negativ);

    if(ciclu_negativ){
        out<<"Ciclu negativ!";
    }
    else{
        for(int i=2; i<=nrNoduri; i++){
            out<<sol[i]<<' ';
        }
    }


    in.close();
    out.close();
}

void Graf::tarjan(int nod_curent, vector<bool>& viz, vector<bool>& pe_stiva, vector<int>& low, stack<int>& stiva, vector<vector<int>>& sol)
{
    viz[nod_curent] = true;
    stiva.push(nod_curent);
    pe_stiva[nod_curent] = true;

    for(int i=0; i<int(vecini[nod_curent].size()); i++){
        int vecin = vecini[nod_curent][i];
        if(!viz[vecin]){
            tarjan(vecin, viz, pe_stiva, low, stiva, sol);
        }
        if(pe_stiva[vecin]){
            low[vecin] = min(low[vecin], low[nod_curent]);
        }
    }

    /*if(low[nod_curent] == nod_curent){
        vector<int> componenta_conexa;
        int x = stiva.top();
        while(low[x] == nod_curent){
            componenta_conexa.push_back(x);
            pe_stiva[x] = false;
            stiva.pop();
            x = stiva.top();
        }
        sol.push_back(componenta_conexa);
    }*/
}

vector<vector<int>> Graf::ctc()
{
    stack<int> stiva;
    vector<bool> viz(nrNoduri+1, false), pe_stiva(nrNoduri+1, false);
    vector<int> low(nrNoduri+1);
    vector<vector<int>> sol;

    for(int i=1; i<=nrNoduri; i++){
        low[i] = i;
    }

    tarjan(1, viz, pe_stiva, low, stiva, sol);
    /*for(int i=0; i<int(sol.size()); i++){
        for(int j=0; j<int(sol[i].size()); j++){
            cout<<sol[i][j]<<' ';
        }
        cout<<'\n';
    }*/
    return sol;
}

void Graf::problema_ctc()
{
    ifstream in("ctc.in");
    ofstream out("ctc.out");

    in>>nrNoduri>>nrMuchii;
    vecini.resize(nrNoduri+1);
    for(int i=1; i<=nrMuchii; i++){
        int x,y;
        in>>x>>y;
        vecini[x].push_back(y);
    }
    afis_lista_ad();
    vector<vector<int>> sol;
    ctc();

    in.close();
    out.close();
}

int main()
{
    Graf G;
    //G.problema_BFS();
    //G.problema_DFS();
    //G.problema_ctc();
    //G.problema_sortaret();
    //G.Havel_Hakimi();
    G.problema_APM();
    //G.Paduri_de_multimi_disjuncte();
    //G.problema_Dijkstra();
    //G.problema_Bellman_Ford();
    //G.problema_Roy_Floyd();
    //G.problema_darb();

    return 0;
}

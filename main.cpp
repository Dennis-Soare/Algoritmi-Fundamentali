#include <bits/stdc++.h>

using namespace std;

const int NMAX = numeric_limits<int>::max();


struct muchie {
    int sursa, destinatie, capacitate;
};

class Graf {
    int nr_noduri, nr_muchii;
    vector<vector<muchie>> vecini;


    bool construieste_lant_nesat_BF(const int sursa, const int dest, vector<int>& tata,const vector<vector<int>>& cap, vector<vector<int>>& flux);
    void parcurgere_dfs(const int nod_curent, vector<int>& sol, vector<bool>& viz);     //functie ajutatoare pentru dfs
    vector<vector<int>> ctc();
    void tarjan(int nod_curent, vector<bool>& viz, vector<bool>& pe_stiva, vector<int>& low, stack<int>& stiva, vector<vector<int>>& sol);
    void reuneste(int x, int y, vector<int>& tata, vector<int>& inaltime);  // pentru paduri_de_multimi_disjuncte
    int reprez(int x, vector<int>& tata);                                   // pentru paduri_de_multimi_disjuncte

public:
    Graf(){};
    Graf(int nr_noduri);
    void adauga_muchie(int, int, int);
    void afis_lista_ad();
    ~Graf();

    vector<int> bfs(const int sursa);
    vector<int> dfs(const int sursa);
    vector<int> distante_bfs(const int sursa);      // calculeaza distanta de la sursa la fiecare nod
    int diametru_arbore();
    vector<int> sortare_topologica();
    bool havel_hakimi(int grad[]);
    vector<pair<int,int>> prim(int& cost_arbore);
    vector<int> dijkstra(const int sursa);
    vector<int> bellman_ford(const int sursa, bool& ciclu_negativ);
    int flux_maxim(const int sursa, const int dest);
    vector<vector<int>> roy_floyd(vector<vector<int>> matrice);
    vector<int> ciclu_eulerian(vector<vector<pair<int,int>>>&);
    void paduri_de_multimi_disjuncte();

};

Graf::Graf(int nr_noduri)
{
    this -> nr_noduri = nr_noduri;
    vecini.resize(nr_noduri + 1);
    this -> nr_muchii = 0;
}

Graf::~Graf()
{
    vecini.clear();
}

void Graf::adauga_muchie(int sursa, int destinatie, int capacitate = 1)
{
    nr_muchii++;
    muchie m;
    m.sursa = sursa;
    m.destinatie = destinatie;
    m.capacitate = capacitate;
    vecini[sursa].push_back(m);
}

void Graf::afis_lista_ad()
{
    for(int i=1; i<=nr_noduri; i++) {
        cout << i << ": ";
        for(int j=0; j<int(vecini[i].size()); j++) {
            cout << vecini[i][j].destinatie << ' '<< vecini[i][j].capacitate << "; ";
        }
        cout <<'\n';
    }
}


vector<int> Graf::bfs(const int sursa)                  // O(n+m)
{
    vector<int> sol;
    queue<int> coada;
    vector<bool> vizitat(nr_noduri+1, false);

    coada.push(sursa);
    vizitat[sursa] = true;

    while(!coada.empty()) {
        int nod_curent = coada.front();
        coada.pop();
        sol.push_back(nod_curent);

        for(int i=0; i< int(vecini[nod_curent].size()); i++) {
            int x = vecini[nod_curent][i].destinatie;
            if(!vizitat[x]) {
                vizitat[x] = true;
                coada.push(x);
            }
        }
    }
    return sol;
}

vector<int> Graf::distante_bfs(const int sursa)         // O(n+m)
{
    queue<int> coada;
    vector<int> dist(nr_noduri+1, 0);
    vector<bool> vizitat(nr_noduri+1, false);

    coada.push(sursa);
    vizitat[sursa] = true;

    while(!coada.empty()) {
        int nod_curent = coada.front();
        coada.pop();

        for(int i=0; i< int(vecini[nod_curent].size()); i++) {
            int x = vecini[nod_curent][i].destinatie;

            if(!vizitat[x]) {
                vizitat[x] = true;
                coada.push(x);
                dist[x] = dist[nod_curent] + 1;
            }
        }
    }
    for(int i=1; i<=nr_noduri; i++) {
        if(dist[i] == 0 && i != sursa)
            dist[i] = -1;
    }
    return dist;
}

void Graf::parcurgere_dfs(const int nod_curent, vector<int>& sol, vector<bool>& viz)
{
    viz[nod_curent] = true;
    sol.push_back(nod_curent);
    for(int i=0; i<int(vecini[nod_curent].size()); i++) {
        int nod = vecini[nod_curent][i].destinatie;
        if(!viz[nod]) {
            parcurgere_dfs(nod, sol , viz);
        }
    }
}

vector<int> Graf::dfs(const int sursa)                  // O(n+m)
{
    vector<int> sol;
    vector<bool> viz(nr_noduri+1, false);
    parcurgere_dfs(sursa, sol, viz);
    return sol;
}

int Graf::diametru_arbore()                             // O(n+m)
{
    vector<int> elemente_bfs = bfs(1);
    vector<int> dist = distante_bfs(elemente_bfs[nr_noduri-1]);       //calculam distantele nodurilor de la ultimul nod gasit in bfs

    int sol = 0;
    for(int i=1; i<=nr_noduri; i++){         //determinam distanta dintre ultimul nod din bfs si cel mai indepartat nod de el
        if(dist[i]+1 > sol)
            sol = dist[i]+1;
    }

    return sol;
}

vector<int> Graf::sortare_topologica()
{
    vector<int> grad_intern(nr_noduri + 1, 0);

    for(int i=1; i<=nr_noduri; i++) {
        for(int j=0; j<int(vecini[i].size()); j++) {
            grad_intern[vecini[i][j].destinatie]++;
        }
    }

    vector<int> sol;
    int contor = 0;
    while(contor < nr_noduri) {
        for(int nod = 1; nod <= nr_noduri; nod++) {
            if(grad_intern[nod] == 0) {         //selectam nodurile cu gradul intern 0
                contor++;
                sol.push_back(nod);
                for(int i=0; i<int(vecini[nod].size()); i++){       //updatam gradul intern al vecinilor dupa ce am scos nodul selectat
                    grad_intern[vecini[nod][i].destinatie]--;
                }
                grad_intern[nod]--;
            }
        }
    }
    return sol;
}

bool Graf::havel_hakimi(int grad[])
{
    int suma=0;
    for(int i=0; i<nr_noduri; i++) {
        suma += grad[i];
        if(grad[i] > (nr_noduri - 1)) {
            return false;
        }
    }

    if(suma % 2 == 1){
        return false;
    }

    sort(grad, grad+nr_noduri, greater<int>());

    while(grad[0] != 0){
        for(int i=1; i<=grad[0]; i++){
            grad[i]--;
            if(grad[i] < 0){
                return false;
            }
        }
        grad[0] = 0;
        sort(grad, grad+nr_noduri, greater<int>());
    }
    return true;
}

vector<int> Graf::dijkstra(const int sursa)             // O(m * log(n))
{
    vector<int> dist(nr_noduri+1, NMAX);
    vector<bool> viz(nr_noduri+1, false);
    priority_queue<pair<int,int>, vector<pair<int,int>>, greater<pair<int,int>>> q;

    q.push({0, sursa});
    dist[sursa] = 0;

    while(!q.empty()){
        int nod = q.top().second;
        q.pop();

        if(!viz[nod]){
            viz[nod] = true;
            for(auto vecin: vecini[nod]){             //parcurgem toti vecinii nodului curent
                int nod_vecin = vecin.destinatie;
                int cost_vecin = vecin.capacitate;

                if(dist[nod] + cost_vecin < dist[nod_vecin]){       //actualizam distantele
                    dist[nod_vecin] = dist[nod] + cost_vecin;
                    q.push({dist[nod_vecin], nod_vecin});
                }
            }
        }
    }
    return dist;
}

vector<pair<int,int>> Graf::prim(int& cost_arbore)      // O(m * log(n))
{
    vector<bool> viz(nr_noduri+1, false);
    vector<int> tata(nr_noduri+1, -1), cost_min(nr_noduri+1, NMAX);
    priority_queue<pair<int,int>, vector<pair<int,int>>, greater<pair<int,int>>> q;

    q.push({0, 1});
    cost_min[1] = 0;

    while(!q.empty()){
        int nod = q.top().second;
        int cost_muchie = q.top().first;
        q.pop();

        if(!viz[nod]){              // alegem muchia cea mai mica care pleaca din componenta conexa si se termina in afara ei
                                    // nodul proaspat unit devine parte din componenta conexa
            viz[nod] = true;
            cost_arbore += cost_muchie;

            for(auto vecin: vecini[nod]){         // parcurgem muchiile nodului nou
                int nod_vecin = vecin.destinatie;
                int cost_vecin = vecin.capacitate;

                if(!viz[nod_vecin] && cost_min[nod_vecin] > cost_vecin){  // daca nu am vizitat nodul si muchia gasita are un cost mai bun
                    cost_min[nod_vecin] = cost_vecin;                     // updatam costul minim si introducem costul si nodul nou in heap
                    q.push({cost_vecin, nod_vecin});
                    tata[nod_vecin] = nod;
                }
            }
        }
    }

    vector<pair<int,int>> sol;
    for(int i=1; i<=nr_noduri; i++){
        sol.push_back({tata[i],i});
    }
    return sol;

}

vector<int> Graf::bellman_ford(const int sursa, bool& ciclu_negativ)        // O(n*m)
{
    ciclu_negativ = false;
    vector<int> dist(nr_noduri+1, NMAX), viz(nr_noduri+1, 0);
    priority_queue<pair<int,int>, vector<pair<int,int>>, greater<pair<int,int>>> q;

    q.push({0, sursa});
    dist[sursa] = 0;

    while(!q.empty()){
        int nod = q.top().second;
        q.pop();
        viz[nod]++;

        if (viz[nod] == nr_noduri){                  // daca se mai fac actualizari la pasul n, inseamna ca exista un ciclu negativ in graf
            dist.resize(0);
            ciclu_negativ = true;
            break;
        }

        for(auto vecin: vecini[nod]){             //parcurgem toti vecinii nodului curent
            int nod_vecin = vecin.destinatie;
            int cost = vecin.capacitate;

            if(dist[nod] + cost < dist[nod_vecin]){         //relaxam muchia
                dist[nod_vecin] = dist[nod] + cost;
                q.push({dist[nod_vecin], nod_vecin});
            }
        }
    }
    return dist;
}

bool Graf::construieste_lant_nesat_BF(const int sursa, const int dest, vector<int>& tata, const vector<vector<int>>& cap, vector<vector<int>>& flux)
{
    queue<int> coada;           // pentru bfs
    bool am_ajuns = false;      // bool care verifica daca bfs a ajuns la nodul destinatie

    for(int i=1; i<=nr_noduri; i++){
        tata[i] = 0;
    }

    tata[sursa] = -1;           // pentru a marca nodul sursa vizitat
    coada.push(sursa);

    while(!coada.empty()){
        int i = coada.front();
        coada.pop();

        if(i == dest){
            am_ajuns = true;
            continue;
        }

        for(int j=0; j<int(vecini[i].size()); j++){
            int vecin = vecini[i][j].destinatie;
            if((cap[i][vecin] - flux[i][vecin] > 0) && tata[vecin] == 0){
                coada.push(vecin);
                tata[vecin] = i;
            }
        }
    }
    return am_ajuns;
}

int Graf::flux_maxim(const int sursa, const int dest)
{
    vector<int> tata(nr_noduri+1, 0);
    vector<vector<int>> flux(nr_noduri+1, vector<int>(nr_noduri+1, 0));
    vector<vector<int>> cap(nr_noduri+1, vector<int>(nr_noduri+1, 0));
    int flux_max = 0;

    for(int i=1; i<=nr_noduri; i++) {
        for(int j=0; j<int(vecini[i].size()); j++) {
            int vecin = vecini[i][j].destinatie;
            cap[i][vecin] = vecini[i][j].capacitate;
        }
    }

    // Algoritm Edmonds-Karp
    // O(n * m^2)

   while(construieste_lant_nesat_BF(sursa, dest, tata, cap, flux)){
        for(int i=0; i<int(vecini[dest].size()); i++){           // Revizuim toate drumurile gasite de la sursa la destinatie dupa o parcurgere BF
            int nod = vecini[dest][i].destinatie;
            int u = dest;
            int update_flux = NMAX;

            if(tata[nod] == 0){
                continue;
            }

            tata[dest] = nod;
            while(u != sursa){
                update_flux = min(update_flux, cap[tata[u]][u] - flux[tata[u]][u]);
                u = tata[u];
            }

            u = dest;
            while(u != sursa){              // Revizuieste flux lant
                flux[tata[u]][u] += update_flux;
                flux[u][tata[u]] -= update_flux;
                u = tata[u];
            }

            flux_max += update_flux;
        }

    }

    return flux_max;
}

vector<vector<int>> Graf::roy_floyd(vector<vector<int>> matrice)        // O(n^3)
{
    for(int i=0; i<nr_noduri; i++){
        for(int j=0; j<nr_noduri; j++){
                if(matrice[i][j] == 0 && i!=j)          //daca nu exista muchie intre i si j
                    matrice[i][j] = 1e9;
        }
    }

    for(int k=0; k<nr_noduri; k++){
        for(int i=0; i<nr_noduri; i++){
            for(int j=0; j<nr_noduri; j++){
                if (matrice[i][k] + matrice[k][j] < matrice[i][j])
                        matrice[i][j] = matrice[i][k] + matrice[k][j];
            }
        }
    }

    for(int i=0; i<nr_noduri; i++){
        for(int j=0; j<nr_noduri; j++){
            if(matrice[i][j] == 1e9)
                matrice[i][j] = 0;
        }
    }
    return matrice;
}

vector<int> Graf::ciclu_eulerian(vector<vector<pair<int,int>>>& lista_ad)
{
    vector<int> ciclu;
    vector<bool> viz_muchie(nr_muchii+1, false);
    stack<int> noduri_nesaturate;

    for(int i=1; i<=nr_noduri; i++) {
        if(int(lista_ad[i].size()) % 2 == 1){      // daca un nod are gradul impar graful nu este eulerian
            ciclu.push_back(-1);
            return ciclu;
        }
    }

    noduri_nesaturate.push(1);

    while(!noduri_nesaturate.empty()) {
        int nod = noduri_nesaturate.top();

        if(lista_ad[nod].empty()) {
            noduri_nesaturate.pop();
            ciclu.push_back(nod);
        }
        else {
            int vecin = lista_ad[nod].back().first;
            int nr_muchie = lista_ad[nod].back().second;
            lista_ad[nod].pop_back();

            if(!viz_muchie[nr_muchie]) {
                viz_muchie[nr_muchie] = true;
                noduri_nesaturate.push(vecin);
            }
        }
    }

    ciclu.pop_back();
    return ciclu;
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

void Graf::paduri_de_multimi_disjuncte()
{
    ifstream in("disjoint.in");
    ofstream out("disjoint.out");

    int nrOperatii;
    in>>nr_noduri>>nrOperatii;

    vector<int> tata(nr_noduri+1, 0), inaltime(nr_noduri+1, 0);

    // fiecare multime este memorata sub forma unui arbore, avand ca reprezentant radacina

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


void problema_bfs()
{
    ifstream in("bfs.in");
    ofstream out("bfs.out");

    int n, m, s;
    in >> n >> m >> s;

    Graf g(n);
    for(int i=0; i<m; i++) {
        int x, y;
        in >> x >> y;
        g.adauga_muchie(x,y);
    }

    vector<int> sol = g.distante_bfs(s);

    for(int i=1; i<=n; i++) {
        out << sol[i] <<' ';
    }

    in.close();
    out.close();
}

void problema_dfs()
{
    ifstream in("dfs.in");
    ofstream out("dfs.out");

    int n, m;
    in>> n >> m;

    Graf g(n);
    for(int i=0; i<m; i++){
        int x,y;
        in>>x>>y;
        g.adauga_muchie(x,y);
    }


    vector<bool> viz(n+1, false);
    int nr_comp_conexe = 0;

    for(int i=1; i<=n; i++) {
        if(!viz[i]) {
            nr_comp_conexe++;
            vector<int> dfs = g.bfs(i);
            for(auto x: dfs) {
                viz[x] = true;
            }

        }
    }
    out << nr_comp_conexe;

    in.close();
    out.close();
}

void problema_darb()
{
    ifstream in("darb.in");
    ofstream out("darb.out");

    int nr_noduri;

    in>>nr_noduri;
    Graf g(nr_noduri);

    for(int i=1;i<nr_noduri;i++){
        int x,y;
        in>>x>>y;

        g.adauga_muchie(x,y);
        g.adauga_muchie(y,x);
    }


    int sol = g.diametru_arbore();
    out<<sol;


    in.close();
    out.close();
}

void problema_sortaret()
{
    ifstream in("sortaret.in");
    ofstream out("sortaret.out");

    int nr_noduri, nr_muchii;
    in>>nr_noduri>>nr_muchii;
    Graf g(nr_noduri);

    for(int i=0; i<nr_muchii; i++){
        int x,y;
        in>>x>>y;
        g.adauga_muchie(x,y);
    }

    vector<int> sol = g.sortare_topologica();
    for(auto x: sol) {
        out << x << ' ';
    }

    in.close();
    out.close();
}

void problema_havel_hakimi()
{
    int nr_noduri;

    cout<<"Introduceti numarul de noduri: ";
    cin>>nr_noduri;
    Graf g(nr_noduri);

    int grad[nr_noduri+1];
    cout<<"Introduceti " << nr_noduri <<" numere: ";

    for(int i=0; i<nr_noduri; i++) {
        cin >> grad[i];
    }

    bool HH = g.havel_hakimi(grad);

    if(HH)
        cout<<"DA";
    else cout<<"NU";
}

void problema_dijkstra()
{
    ifstream in("dijkstra.in");
    ofstream out("dijkstra.out");

    int nr_noduri, nr_muchii;
    in>>nr_noduri>>nr_muchii;
    Graf g(nr_noduri);

    for(int i=1; i<=nr_muchii; i++){
        int x,y,cost;
        in>>x>>y>>cost;
        g.adauga_muchie(x,y,cost);
    }

    vector<int> sol = g.dijkstra(1);

    for(int i=2; i<=nr_noduri; i++){
        if(sol[i]==NMAX)
            out<<"0 ";
        else
            out<<sol[i]<<' ';
    }

    in.close();
    out.close();
}

void problema_apm()
{
    ifstream in("apm.in");
    ofstream out("apm.out");

    int nr_muchii, nr_noduri;
    in>>nr_noduri>>nr_muchii;
    Graf g(nr_noduri);

    for(int i=1; i<=nr_muchii; i++){
        int x,y,cost;
        in>>x>>y>>cost;
        g.adauga_muchie(x,y,cost);
        g.adauga_muchie(y,x,cost);
    }


    int cost_arbore = 0;
    vector<pair<int,int>> sol = g.prim(cost_arbore);

    out<<cost_arbore<<'\n'<<nr_noduri-1<<'\n';

    for(int i=1; i<=nr_noduri-1; i++){
        out<<sol[i].first<<' '<<sol[i].second<<'\n';
    }

    in.close();
    out.close();
}

void problema_bellman_ford()
{
    ifstream in("bellmanford.in");
    ofstream out("bellmanford.out");

    int nr_noduri, nr_muchii;
    in>>nr_noduri>>nr_muchii;
    Graf g(nr_noduri);

    for(int i=1; i<=nr_muchii; i++){
        int x,y,cost;
        in>>x>>y>>cost;
        g.adauga_muchie(x,y,cost);
    }

    bool ciclu_negativ;
    vector<int> sol = g.bellman_ford(1, ciclu_negativ);

    if(ciclu_negativ){
        out<<"Ciclu negativ!";
    }
    else{
        for(int i=2; i<=nr_noduri; i++){
            out<<sol[i]<<' ';
        }
    }

    in.close();
    out.close();
}

void problema_flux_maxim()
{
    ifstream in("maxflow.in");
    ofstream out("maxflow.out");

    int nr_muchii, nr_noduri;
    in>>nr_noduri>>nr_muchii;
    Graf g(nr_noduri);

    for(int i=1; i<=nr_muchii; i++){
            int x,y,z;
            in>>x>>y>>z;
            g.adauga_muchie(x,y,z);
            g.adauga_muchie(y,x,0);
    }

    int maxflow = g.flux_maxim(1, nr_noduri);
    out<<maxflow;

    in.close();
    out.close();
}

void problema_roy_floyd()
{
    ifstream in("royfloyd.in");
    ofstream out("royfloyd.out");

    int nr_noduri;
    in>>nr_noduri;
    Graf g(nr_noduri);
    vector<vector<int>> matrice(nr_noduri);

    for(int i=0; i<nr_noduri; i++){
        for(int j=0; j<nr_noduri; j++){
            int x;
            in>>x;
            matrice[i].push_back(x);
        }
    }

    vector<vector<int>> sol = g.roy_floyd(matrice);

    for(int i=0; i<nr_noduri; i++){
        for(int j=0; j<nr_noduri; j++){
                out<< sol[i][j]<<' ';
        }
        out<<'\n';
    }

    in.close();
    out.close();
}

void problema_ciclu_eulerian()
{
    ifstream in("ciclueuler.in");
    ofstream out("ciclueuler.out");

    int nr_noduri, nr_muchii;
    in>>nr_noduri>>nr_muchii;
    Graf g(0);
    vector<vector<pair<int,int>>> lista_ad(nr_noduri+1);

    for(int i=1; i<=nr_muchii; i++){
            int x,y;
            in>>x>>y;
            lista_ad[x].push_back({y, i});
            lista_ad[y].push_back({x, i});
    }

    vector<int> sol = g.ciclu_eulerian(lista_ad);

    for(auto el: sol)
        out << el << ' ';

    in.close();
    out.close();
}

int main()
{
    //problema_bfs();
    //problema_dfs();
    //problema_darb();
    //problema_sortaret();
    //problema_havel_hakimi();
    //problema_dijkstra();
    //problema_apm();
    //problema_bellman_ford();
    //problema_flux_maxim();
    //problema_roy_floyd();
    //problema_ctc();
    //problema_ciclu_eulerian();
    //problema_cuplaj();


    return 0;
}

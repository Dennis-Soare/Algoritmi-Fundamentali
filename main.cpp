#include <bits/stdc++.h>

using namespace std;

const int NMAX = numeric_limits<int>::max();

class Graf {
    int nr_noduri, nr_muchii;
    vector<vector<int>> vecini;


    bool construieste_lant_nesat_BF(const int sursa, const int dest, vector<int>& tata,const vector<vector<int>>& cap, vector<vector<int>>& flux);
    vector<int> distante_bfs(const int Start);      // calculeaza distanta de la start la fiecare nod
    void parcurgere_dfs(const int nod_curent, vector<int>& sol, vector<bool>& viz);     //functie ajutatoare pentru dfs
    vector<vector<int>> ctc();
    void tarjan(int nod_curent, vector<bool>& viz, vector<bool>& pe_stiva, vector<int>& low, stack<int>& stiva, vector<vector<int>>& sol);
    void reuneste(int x, int y, vector<int>& tata, vector<int>& inaltime);  // pentru paduri_de_multimi_disjuncte
    int reprez(int x, vector<int>& tata);                                   // pentru paduri_de_multimi_disjuncte
    void afis_lista_ad();

public:
    vector<int> bfs(const int Start);
    vector<int> dfs(const int Start);
    vector<int> sortare_topologica(vector<int> grad_intern);
    bool havel_hakimi(int grad[]);
    vector<int> dijkstra(const vector<vector<pair<int, int>>> &lista_ad, const int sursa);
    vector<pair<int,int>> prim(const vector<vector<pair<int, int>>> &lista_ad, int& cost_arbore, const int sursa);
    vector<int> bellman_ford(const vector<vector<pair<int, int>>>&, const int sursa, bool& ciclu_negativ);
    vector<vector<int>> roy_floyd(vector<vector<int>> matrice);
    int diametru_arbore();
    int flux_maxim(const vector<vector<int>>& cap, const int sursa, const int dest);
    vector<int> ciclu_eulerian(vector<vector<pair<int,int>>>& lista_ad, const int sursa);


    void problema_bfs();
    void problema_dfs();
    void problema_sortaret();
    void problema_havel_hakimi();
    void problema_apm();
    void paduri_de_multimi_disjuncte();
    void problema_dijkstra();
    void problema_bellman_ford();
    void problema_roy_floyd();
    void problema_darb();
    void problema_ctc();
    void problema_flux_maxim();
    void problema_ciclu_eulerian();
    void problema_cuplaj();
};



vector<int> Graf::ciclu_eulerian(vector<vector<pair<int,int>>>& lista_ad, const int sursa)
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

    noduri_nesaturate.push(sursa);

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

void Graf::problema_ciclu_eulerian()
{
    ifstream in("ciclueuler.in");
    ofstream out("ciclueuler.out");

    in>>nr_noduri>>nr_muchii;
    vector<vector<pair<int,int>>> lista_ad(nr_noduri+1);

    for(int i=1; i<=nr_muchii; i++){
            int x,y;
            in>>x>>y;
            lista_ad[x].push_back({y, i});
            lista_ad[y].push_back({x, i});
    }

    vector<int> sol = ciclu_eulerian(lista_ad, 1);

    for(auto el: sol)
        out << el << ' ';

    in.close();
    out.close();
}

void Graf::problema_cuplaj()
{
    ifstream in("cuplaj.in");
    ofstream out("cuplaj.out");

    int card_l, card_r, e;
    in >> card_l >> card_r >> e;

    nr_noduri = card_l + card_r + 2;         // cardinalele multimilor + nodul sursa + nodul destinatie
    nr_muchii = card_l + card_l + e;
    vecini.resize(nr_noduri + 1);
    vector<vector<int>> cap(nr_noduri + 1, vector<int>(nr_noduri + 1, 0));       //matrice cu capacitati

    int nod_sursa = card_l + card_r + 1;
    int nod_dest = card_l + card_r + 2;

    for(int i=1; i<=card_l; i++) {
        cap[nod_sursa][i] = 1;
        vecini[nod_sursa].push_back(i);
        vecini[i].push_back(nod_sursa);
    }
    for(int i=card_l+1; i <= card_l + card_r; i++) {
        cap[i][nod_dest] = 1;
        vecini[i].push_back(nod_dest);
        vecini[nod_dest].push_back(i);
    }

    for(int i=1; i<=e; i++) {
        int u, v;
        in >> u >> v;
        vecini[u].push_back(card_l + v);
        vecini[card_l + v].push_back(u);
        cap[u][card_l + v] = 1;
    }

    /*afis_lista_ad();
    cout<<endl;

    for(int i=1; i<=nr_noduri; i++) {
        for(int j=1; j<=nr_noduri; j++) {
            cout << cap[i][j] <<' ';
        }
        cout << '\n';
    }
    cout <<endl;*/

    int x = flux_maxim(cap, nod_sursa, nod_dest);
    out <<x;

    in.close();
    out.close();
}

bool Graf::construieste_lant_nesat_BF(const int sursa, const int dest, vector<int>& tata,const vector<vector<int>>& cap, vector<vector<int>>& flux)
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

        for(int &j: vecini[i]){
            if((cap[i][j] - flux[i][j] > 0) && tata[j] == 0){
                coada.push(j);
                tata[j] = i;
            }
        }
    }

    return am_ajuns;
}

int Graf::flux_maxim(const vector<vector<int>>& cap, const int sursa, const int dest)
{
    vector<int> tata(nr_noduri+1, 0);
    vector<vector<int>> flux(nr_noduri+1, vector<int>(nr_noduri+1, 0));
    int flux_max = 0;

    // Algoritm Edmonds-Karp
    // O(n * m^2)

   while(construieste_lant_nesat_BF(sursa, dest, tata, cap, flux)){
        for(auto &nod: vecini[dest]){           // Revizuim toate drumurile gasite de la sursa la destinatie dupa o parcurgere BF
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

void Graf::problema_flux_maxim()
{
    ifstream in("maxflow.in");
    ofstream out("maxflow.out");

    in>>nr_noduri>>nr_muchii;
    vector<vector<int>> cap(nr_noduri+1, vector<int>(nr_noduri+1, 0));       //matrice cu capacitati

    vecini.resize(nr_noduri+1);

    for(int i=1; i<=nr_muchii; i++){
            int x,y,z;
            in>>x>>y>>z;
            vecini[x].push_back(y);
            vecini[y].push_back(x);
            cap[x][y] = z;
    }

       cout<<endl;

    for(int i=1; i<=nr_noduri; i++) {
        for(int j=1; j<=nr_noduri; j++) {
            cout << cap[i][j] <<' ';
        }
        cout << '\n';
    }

    int maxflow = flux_maxim(cap, 1, nr_noduri);
    cout<<maxflow;

    in.close();
    out.close();
}

int Graf::diametru_arbore()
{
    vector<int> bfs = bfs(1);
    vector<int> dist = distante_bfs(bfs[nr_noduri-1]);       //calculam distantele nodurilor de la ultimul nod gasit in bfs

    int sol = 0;
    for(int i=1; i<=nr_noduri; i++){         //determinam distanta dintre ultimul nod din bfs si cel mai indepartat nod de el
        if(dist[i]+1 > sol)
            sol = dist[i]+1;
    }

    return sol;
}

void Graf::problema_darb()
{
    ifstream in("darb.in");
    ofstream out("darb.out");

    in>>nr_noduri;
    vecini.resize(nr_noduri+1);
    for(int i=1;i<nr_noduri;i++){
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

void Graf::problema_roy_floyd()
{
    ifstream in("royfloyd.in");
    ofstream out("royfloyd.out");

    in>>nr_noduri;
    vector<vector<int>> matrice(nr_noduri);

    for(int i=0; i<nr_noduri; i++){
        for(int j=0; j<nr_noduri; j++){
            int x;
            in>>x;
            matrice[i].push_back(x);
        }
    }

    vector<vector<int>> sol = roy_floyd(matrice);

    for(int i=0; i<nr_noduri; i++){
        for(int j=0; j<nr_noduri; j++){
                out<< sol[i][j]<<' ';
        }
        out<<'\n';
    }

    in.close();
    out.close();
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

void Graf::afis_lista_ad()
{
    for(int i=1; i<=nr_noduri; i++)
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

vector<int> Graf::distante_bfs(const int Start)         // O(n+m)
{
    queue<int> coada;
    vector<int> dist(nr_noduri+1, 0);
    vector<bool> vizitat(nr_noduri+1, false);

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
    for(int i=1; i<=nr_noduri; i++) {
        if(dist[i] == 0 && i != Start)
            dist[i] = -1;
    }
    return dist;
}

vector<int> Graf::bfs(const int Start)                  // O(n+m)
{
    vector<int>sol;
    queue<int> coada;
    vector<bool> vizitat(nr_noduri+1, false);

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

void Graf::problema_bfs()
{
    ifstream in("bfs.in");
    ofstream out("bfs.out");

    int Start;
    in>>nr_noduri>>nr_muchii>>Start;

    vecini.resize(nr_noduri+1);
    for(int i=0; i<nr_muchii; i++)
    {
        int x,y;
        in>>x>>y;
        vecini[x].push_back(y);
    }
    vector<int> distante = distante_bfs(Start);
    for(int i=1; i<=nr_noduri; i++) {
        out << distante[i]<< ' ';
    }

    in.close();
    out.close();
}

void Graf::parcurgere_dfs(const int nod_curent, vector<int>& sol, vector<bool>& viz)
{
    viz[nod_curent] = true;
    sol.push_back(nod_curent);
    for(int i=0; i<int(vecini[nod_curent].size()); i++) {
        int nod = vecini[nod_curent][i];
        if(!viz[nod]) {
            parcurgere_dfs(nod, sol , viz);
        }
    }
}

vector<int> Graf::dfs(const int Start)                  // O(n+m)
{
    vector<int> sol;
    vector<bool> viz(nr_noduri+1, false);
    parcurgere_dfs(Start, sol, viz);
    return sol;
}

void Graf::problema_dfs()
{
    ifstream in("dfs.in");
    ofstream out("dfs.out");

    in>>nr_noduri>>nr_muchii;
    vecini.resize(nr_noduri+1);

    for(int i=0; i<nr_muchii; i++){
        int x,y;
        in>>x>>y;
        vecini[x].push_back(y);
        vecini[y].push_back(x);
    }

    vector<int> nr_comp_conexa(nr_noduri+1, 0);

    int sol = 0;

    for(int i=1; i<=nr_noduri; i++){
        if(nr_comp_conexa[i] == 0){
            sol++;
            vector<int> dfs = dfs(i);
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
    while(contor < nr_noduri) {
        for(int nod = 1; nod <= nr_noduri; nod++) {
            if(grad_intern[nod] == 0) {         //selectam nodurile cu gradul intern 0
                contor++;
                sol.push_back(nod);
                for(int i=0; i<int(vecini[nod].size()); i++){       //updatam gradul intern al vecinilor dupa ce am scos nodul selectat
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

    in>>nr_noduri>>nr_muchii;

    vecini.resize(nr_noduri+1);
    vector<int> grad_intern(nr_noduri+1, 0);

    for(int i=0; i<nr_muchii; i++){
        int x,y;
        in>>x>>y;
        vecini[x].push_back(y);
        grad_intern[y]++;
    }

    vector<int> sol = sortare_topologica(grad_intern);
    for(int i=0; i<nr_noduri; i++){
        out<<sol[i]<<' ';
    }
    in.close();
    out.close();
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

void Graf::problema_havel_hakimi()
{
    cout<<"Introduceti numarul de noduri: ";
    cin>>nr_noduri;

    int grad[nr_noduri+1];
    cout<<"Introduceti " << nr_noduri <<" numere: ";

    for(int i=0; i<nr_noduri; i++) {
        cin >> grad[i];
    }

    bool HH = havel_hakimi(grad);

    if(HH)
        cout<<"DA";
    else cout<<"NU";
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

vector<pair<int,int>> Graf::prim(const vector<vector<pair<int, int>>> &lista_ad, int& cost_arbore, const int sursa)     // O(m * log(n))
{
    vector<bool> viz(nr_noduri+1, false);
    vector<int> tata(nr_noduri+1, -1), cost_min(nr_noduri+1, NMAX);
    priority_queue<pair<int,int>, vector<pair<int,int>>, greater<pair<int,int>>> q;

    q.push({0, sursa});
    cost_min[sursa] = 0;

    while(!q.empty()){
        int nod = q.top().second;
        int cost_muchie = q.top().first;
        q.pop();

        if(!viz[nod]){              // alegem muchia cea mai mica care pleaca din componenta conexa si se termina in afara ei
                                    // nodul proaspat unit devine parte din componenta conexa
            viz[nod] = true;
            cost_arbore += cost_muchie;

            for(auto vecin: lista_ad[nod]){         // parcurgem muchiile nodului nou
                int nod_vecin = vecin.first;
                int cost_vecin = vecin.second;

                if(!viz[nod_vecin] && cost_min[nod_vecin] > cost_vecin){       //daca nu am vizitat nodul si muchia gasita are un cost mai bun
                    cost_min[nod_vecin] = cost_vecin;                          //updatam costul minim si introducem costul si nodul nou in coada cu prioritati
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

void Graf::problema_apm()
{
    ifstream in("apm.in");
    ofstream out("apm.out");


    in>>nr_noduri>>nr_muchii;
    vector<vector<pair<int, int>>> lista_ad(nr_noduri+1);

    for(int i=1; i<=nr_muchii; i++){
        int x,y,cost;
        in>>x>>y>>cost;
        lista_ad[x].push_back({y, cost});
        lista_ad[y].push_back({x, cost});
    }

    /*for(int i=1; i<=nr_noduri;i++){
        cout<<i<<": ";
        for(int j=0; j<int(lista_ad[i].size());j++){
            cout<<lista_ad[i][j].first <<' ' <<lista_ad[i][j].second<<", ";
        }
        cout<<endl;
    }*/

    int cost_arbore = 0;
    vector<pair<int,int>> sol = prim(lista_ad, cost_arbore, 1);

    out<<cost_arbore<<'\n'<<nr_noduri-1<<'\n';

    for(int i=1; i<=nr_noduri-1; i++){
        out<<sol[i].first<<' '<<sol[i].second<<'\n';
    }


    in.close();
    out.close();
}

vector<int> Graf::dijkstra(const vector<vector<pair<int, int>>>& lista_ad, const int sursa)         // O(m * log(n))
{
    vector<int> dist(nr_noduri+1, NMAX);
    vector<bool> viz(nr_noduri+1, false);
    dist[sursa] = 0;
    priority_queue<pair<int,int>, vector<pair<int,int>>, greater<pair<int,int>>> q;
    q.push({0, sursa});

    while(!q.empty()){
        int nod = q.top().second;
        q.pop();

        if(!viz[nod]){
            viz[nod] = true;
            for(auto vecin: lista_ad[nod]){             //parcurgem toti vecinii nodului curent
                int nod_vecin = vecin.first;
                int cost_vecin = vecin.second;
                if(dist[nod] + cost_vecin < dist[nod_vecin]){       //actualizam distantele
                    dist[nod_vecin] = dist[nod] + cost_vecin;
                    q.push({dist[nod_vecin], nod_vecin});
                }
            }
        }
    }
    return dist;
}

void Graf::problema_dijkstra()
{
    ifstream in("dijkstra.in");
    ofstream out("dijkstra.out");

    in>>nr_noduri>>nr_muchii;
    vector<vector<pair<int, int>>> lista_ad(nr_noduri+1);

    for(int i=1; i<=nr_muchii; i++){
        int x,y,cost;
        in>>x>>y>>cost;
        lista_ad[x].push_back(make_pair(y, cost));
    }

    /*for(int i=1; i<=nr_noduri;i++){
        cout<<i<<": ";
        for(int j=0; j<int(lista_ad[i].size());j++){
            cout<<lista_ad[i][j].first <<' ' <<lista_ad[i][j].second<<", ";
        }
        cout<<endl;
    }*/
    vector<int> sol = dijkstra(lista_ad, 1);
    for(int i=2; i<=nr_noduri; i++){
        if(sol[i]==NMAX)
            out<<"0 ";
        else
            out<<sol[i]<<' ';
    }

    in.close();
    out.close();
}

vector<int> Graf::bellman_ford(const vector<vector<pair<int, int>>> &lista_ad, const int sursa, bool& ciclu_negativ)        // O(n*m)
{
    vector<int> dist(nr_noduri+1, NMAX), viz(nr_noduri+1, 0);
    dist[sursa] = 0;
    ciclu_negativ = false;

    priority_queue<pair<int,int>, vector<pair<int,int>>, greater<pair<int,int>>> q;
    q.push({0, sursa});

    while(!q.empty()){
        int nod = q.top().second;
        q.pop();
        viz[nod]++;

        if (viz[nod] == nr_noduri){                  // daca se mai fac actualizari la pasul n, inseamna ca exista un ciclu negativ in graf
            dist.resize(0);
            ciclu_negativ = true;
            break;
        }

        for(auto vecin: lista_ad[nod]){             //parcurgem toti vecinii nodului curent
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

void Graf::problema_bellman_ford()
{
    ifstream in("bellmanford.in");
    ofstream out("bellmanford.out");

    in>>nr_noduri>>nr_muchii;
    vector<vector<pair<int, int>>> lista_ad(nr_noduri+1);

    for(int i=1; i<=nr_muchii; i++){
        int x,y,cost;
        in>>x>>y>>cost;
        lista_ad[x].push_back({y, cost});
    }

    /*for(int i=1; i<=nr_noduri;i++){
        cout<<i<<": ";
        for(int j=0; j<int(lista_ad[i].size());j++){
            cout<<lista_ad[i][j].first <<' ' <<lista_ad[i][j].second<<"; ";
        }
        cout<<endl;
    }*/

    bool ciclu_negativ;
    vector<int> sol = bellman_ford(lista_ad, 1, ciclu_negativ);

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
    vector<bool> viz(nr_noduri+1, false), pe_stiva(nr_noduri+1, false);
    vector<int> low(nr_noduri+1);
    vector<vector<int>> sol;

    for(int i=1; i<=nr_noduri; i++){
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

    in>>nr_noduri>>nr_muchii;
    vecini.resize(nr_noduri+1);
    for(int i=1; i<=nr_muchii; i++){
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
    //G.problema_bfs();
    //G.problema_dfs();
    //G.problema_ctc();
    //G.problema_sortaret();
    //G.problema_havel_hakimi();
    //G.problema_apm();
    //G.paduri_de_multimi_disjuncte();
    //G.problema_dijkstra();
    //G.problema_bellman_ford();
    //G.problema_roy_floyd();
    //G.problema_darb();
    //G.problema_flux_maxim();
    //G.problema_ciclu_eulerian();
    G.problema_cuplaj();

    return 0;
}

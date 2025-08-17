#include<bits/stdc++.h>
using namespace std;


#define int long long 
#define N 1000001

vector<int> g[N];
vector<int> dep(N),par(N);
int up[N][21];
int n;

void dfs(int v,int p) {
    par[v]=p;
    for (auto u:g[v]) {
        if (u!=p) {
            dep[u]=dep[v]+1;
            dfs(u,v);
        }
    }
}

void preprocessLCA() {
    for (int v=1;v<=n;v++) up[v][0]=par[v]; 
    for (int i=1;i<=20;i++) {
        for (int v=1;v<=n;v++) up[v][i]=up[up[v][i-1]][i-1];
    }
}

int lca(int u,int v) {
    if (dep[u]<dep[v]) swap(u,v); 
    int d=dep[u]-dep[v];
    for (int i=0;i<=20;i++) {
        if (d&(1LL<<i)) u=up[u][i];
    }
    if (u==v) return u; 
    for (int i=20;i>=0;i--) {
        if (up[u][i]!=up[v][i]) {
            u=up[u][i]; 
            v=up[v][i];
        }
    }
    return par[u];
}

void soln() {   
    cin>>n;

    for (int i=0;i<n-1;i++) {
        int u,v; cin>>u>>v; 
        g[u].push_back(v); 
        g[v].push_back(u);
    }

    dfs(1,0); 
    preprocessLCA();

    int a,b; cin>>a>>b; 
    cout<<lca(a,b)<<endl;

}   



signed main() {
    ios_base::sync_with_stdio(0);
    cin.tie(0);
    #ifndef ONLINE_JUDGE
        freopen("input.txt","r",stdin);
        freopen("output.txt","w",stdout);
    #endif
    
 

    int tc=1; 
    while (tc--) {
        soln();
    }
    
    
 
    return 0;

}   
 

 

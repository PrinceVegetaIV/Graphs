#include <bits/stdc++.h>
using namespace std;

#define PB push_back
#define MAX 100000 //maximum node
vector<int> graph[MAX];
vector<int> cost[MAX];
int dis[MAX];
bool vis[MAX];

void printGraph(vector<int> Graph[],int nod) {
	for (int i=0;i<nod;i++) {
		cout<<i<<"--> ";
		for (auto x:graph[i]) cout<<x<<" ";
		cout<<endl;
	}
}

void DFS(vector<int> graph[],int source) {
	vis[source]=1;
	for (auto x:graph[source]) {
		if (!vis[x]) DFS(graph,x);
	}
}

void solve() {
	int nod,edg; cin>>nod>>edg;
	for (int i=0;i<edg;i++) {
		int x,y; cin>>x>>y;
		graph[x].PB(y);
		graph[y].PB(x);
	}
	printGraph(graph,nod);
	DFS(graph,0);
	for (int i=0;i<nod;i++) cout<<vis[i]<<" ";
}

 

int main() {
	/*freopen("input.txt","r",stdin);
	freopen("output.txt","w",stdout);*/	

	solve();


	return 0;
}

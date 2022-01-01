#include <bits/stdc++.h>
using namespace std;

#define PB push_back
#define MAX 105 //maximum node
vector <int> graph[MAX];
bool vis[MAX];
vector<int> result; //topologically sorted nodes (from ulta directon)

void printGraph(vector<int> graph[],int nod) {
	for (int i=1;i<=nod;i++) {
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
	result.PB(source); //no child of source left unexplored
}

void solve() {
	int nod,edg; cin>>nod>>edg;
	for (int i=0;i<edg;i++) {
		int x,y; cin>>x>>y;
		graph[x].PB(y);
		//directed graph
	}
	
	for (int i=1;i<=nod;i++) {
		if (!vis[i]) DFS(graph,i);
	}
	
	reverse(result.begin(),result.end()); //now sorted topologically
	
	for (auto x:result) cout<<x<<" ";
}

int main() {
	/*freopen("input.txt","r",stdin);
	freopen("output.txt","w",stdout);*/	

	solve();


	return 0;
}

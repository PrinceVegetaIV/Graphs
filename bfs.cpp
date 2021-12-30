#include <bits/stdc++.h>
using namespace std;

#define PB push_back
#define MAX 100000 //maximum node
vector<int> graph[MAX];
int dis[MAX];
int vis[MAX];

void printGraph(vector<int> graph[],int nod) {
	for (int i=0;i<=nod;i++) {
		cout<<i<<"--> ";
		for (auto x:graph[i]) cout<<x<<" ";
		cout<<endl;
	}
}

void BFS(vector<int> graph[],int source) {
	queue<int> Q;
	Q.push(source);
	vis[source]++;
	dis[source]=0;
	while (!Q.empty()) {
		int ff=Q.front(); 
		Q.pop();
		for (auto x:graph[ff]) {
			if (!vis[x]) {
				Q.push(x);
				dis[x]=dis[ff]+1;
				vis[x]++;
			}
		}
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
	BFS(graph,1);
	for (int i=1;i<=nod;i++) cout<<dis[i]<<" ";
}

 

int main() {
	/*freopen("input.txt","r",stdin);
	freopen("output.txt","w",stdout);*/	

	solve();


	return 0;
}

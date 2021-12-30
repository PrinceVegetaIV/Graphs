#include <bits/stdc++.h>
using namespace std;

#define PB push_back
#define MAX 100000 //maximum node
vector<int> graph[MAX];
int dis[MAX];
bool vis[MAX];

void printGraph(vector<int> graph[],int nod) {
	for (int i=0;i<nod;i++) {
		cout<<i<<"--> ";
		for (auto x:graph[i]) cout<<x<<" ";
		cout<<endl;
	}
} 

void multiBFS(vector<int> graph[],vector<int> sources) {
	queue<int> Q;
	for (auto source:sources) {
		Q.push(source);
		vis[source]=1;
		dis[source]=0;
	}
	while (!Q.empty()) {
		int ff=Q.front();
		Q.pop();
		for (auto x:graph[ff]) {
			if (!vis[x]) {
				vis[x]=1;
				dis[x]=dis[ff]+1;
				Q.push(x);
			}
		}
	}
}

void solve() {
	int nod, edg; cin>>nod>>edg;
	for (int i=0;i<edg;i++) {
		int x,y; cin>>x>>y;
		graph[x].PB(y);
		graph[y].PB(x);
	}
	int numberofSources; cin>>numberofSources;
	vector<int> sources;
	for (int i=0;i<numberofSources;i++) {
		int source; cin>>source;
		sources.PB(source);
	}
	
	multiBFS(graph,sources);
	for (int i=1;i<=nod;i++) {
		if (dis[i]==0) continue; //this is source
		cout<<"node "<<i<<"--> "<<dis[i]<<endl;;
	}
}

int main() {
	/*freopen("input.txt","r",stdin);
	freopen("output.txt","w",stdout);*/	

	solve();


	return 0;
}

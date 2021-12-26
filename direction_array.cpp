#include <bits/stdc++.h>
using namespace std;

#define PB push_back
#define MAX 100000 //maximum node
#define N 5
#define valid(x,y) x>=0 && y>=0 && x<N && y<N && !vis[x][y]
int dx[]={-1,1,0,0,-1,-1,1,1};
int dy[]={0,0,-1,1,-1,1,-1,1}; //king's move
int kx[]={1,1,2,2,-1,-1,-2,-2};
int ky[]={2,-2,1,-1,2,-2,1,-1}; //knight's move
vector<int> graph[MAX];
vector<int> cost[MAX];
int A[N][N];
bool vis[N][N];


void printGraph(vector<int> graph[],int nod) {
	for (int i=0;i<nod;i++) {
		cout<<i<<"--> ";
		for (auto x:graph[i]) cout<<x<<" ";
		cout<<endl;
	}
}

void DFS(int x,int y) {
	vis[x][y]=1;
	for (int i=0;i<8;i++) {
		int nx=x+kx[i];
		int ny=y+ky[i];
		if (valid(nx,ny)) {
			vis[nx][ny]=1;
			DFS(nx,ny);
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
	DFS(2,2);
	for(int i=0;i<N;i++) {
		for (int j=0;j<N;j++) cout<<vis[i][j]<<" ";
		cout<<endl;
	}
}

int main() {
	/*freopen("input.txt","r",stdin);
	freopen("output.txt","w",stdout);*/	

	solve();


	return 0;
}

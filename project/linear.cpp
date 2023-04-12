#include<bits/stdc++.h>
using namespace std;
const int N =10000000;
const int M =30000;
bool vis[N+5];
int p[N+5],cnt;
int main(int argc, char *argv[]){
    int m =atoi(argv[1]);
    for(int i=2;i<=N;++i){
        if(!vis[i])p[++cnt]=i;
        for(int j=1;j<=cnt;++j){
            if(1ll*p[j]*i>N)break;
            vis[p[j]*i]=1;
            if(i%p[j]==0)break;
        }
    }
    for(int i=1;i<=cnt&&p[i]<=m;++i){
        printf("%d\n",p[i]);
    }
}
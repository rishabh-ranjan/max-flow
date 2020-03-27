/*/ Author: kal013 /*/
//#pragma GCC optimize ("O3")
#include <bits/stdc++.h>
#include <cstdio>
using namespace std;
#define int long long
#define ll long long
typedef vector<int> vi;
typedef pair<int,int> pii;
typedef pair<pii,int> ppi;
typedef vector<pii> vpi;
typedef priority_queue<ppi> max_heap_ppi;
#define For(i,n) for(int i=0;i<(int)n;i++)
#define Rev(i,n) for(int i=n-1;i>=0;i--)
#define Rep(i,n) for(int i=1;i<=n;++i)
#define F first
#define S second
#define pb push_back
#define mp make_pair
#define d0(x) cout<<(x)<<" "
#define d1(x) cout<<(x)<<endl
#define d2(x,y) cout<<(x)<<" "<<(y)<<endl
#define d3(x,y,z) cout<<(x)<<" "<<(y)<<" "<<(z)<<endl
#define d4(a,b,c,d) cout<<(a)<<" "<<(b)<<" "<<(c)<<" "<<(d)<<endl
#define d5(a,b,c,d,e) cout<<(a)<<" "<<(b)<<" "<<(c)<<" "<<(d)<<" "<<(e)<<endl
#define d6(a,b,c,d,e,f) cout<<(a)<<" "<<(b)<<" "<<(c)<<" "<<(d)<<" "<<(e)<<" "<<(f)<<endl
#define fio ios_base::sync_with_stdio(false);cin.tie(NULL);cout.tie(NULL)
#define all(v) (v).begin(),(v).end()
#define sz(v) (v).size()

template<class T> ostream& operator<<(ostream &os, vector<T> V) {
    os << "[ ";
    for(auto v : V) os << v << " ";
    return os << "]";
}
template<class T> ostream& operator<<(ostream &os, set<T> S){
    os << "{ ";
    for(auto s:S) os<<s<<" ";
    return os<<"}";
}
template<class L, class R> ostream& operator<<(ostream &os, pair<L,R> P) {
    return os << "(" << P.first << "," << P.second << ")";
}
template<class L, class R> ostream& operator<<(ostream &os, map<L,R> M) {
    os << "{ ";
    for(auto m:M) os<<"("<<m.F<<":"<<m.S<<") ";
    return os<<"}";
}

#define TRACE
#ifdef TRACE
    #define trace(...) __f(#__VA_ARGS__, __VA_ARGS__)
    template <typename Arg1>
    void __f(const char* name, Arg1&& arg1){
        cerr << name << " : " << arg1 << endl;
    }
    template <typename Arg1, typename... Args>
    void __f(const char* names, Arg1&& arg1, Args&&... args){
        const char* comma = strchr(names + 1, ',');
        cerr.write(names, comma - names) << " : " << arg1<<" | ";
        __f(comma+1, args...);
    }
#else
    #define trace(...)
#endif
const int maxn=1e6+100;
const int MaxN=1e5+100;
const int mod=1e9+7;
#ifdef int
const int INF=1e16;
#else
const int INF=1e9;
//#include <ext/pb_ds/assoc_container.hpp> 
//#include <ext/pb_ds/tree_policy.hpp> 
//using namespace __gnu_pbds;
//#define ordered_set tree<int, null_type,less<int>, rb_tree_tag,tree_order_statistics_node_update>
// find_by_order(k)  returns iterator to kth element starting from 0;
// order_of_key(k) returns count of elements strictly smaller than k;
// erase,insert same as normal set
#endif

inline int fast_expo(int base,int power,int modulo=mod){
    base%=modulo;
    if (base<0) base+=modulo;
    ll x=base,cnt=power,ans=1;
    while(cnt){
        if (cnt&1) ans=(ans*x)%modulo;
        x=(x*x)%modulo;
        cnt>>=1;
    }
    // int tmp=ans;
    return ans;
}
inline int inv(int base,int modulo=mod){
    return fast_expo(base,modulo-2,modulo);
}


/*/-----------------------------Code begins----------------------------------/*/


// Code for max flow copied from Far_behind repo.
const ll inf = (1e9);
template<typename T>
struct edge {ll x, y;	T cap, flow;};
template<typename T>
struct DinicFlow {
    // *** change inf accordingly *****
    T eps;
    vector <edge<T>> e;
    vector <ll> cur, d;
    vector < vector <ll> > adj;
    vector <bool> vis;
    ll n, source, sink;
    DinicFlow() {}
    DinicFlow(ll v,T eps) {
        n = v;
        this->eps=eps;
        cur = vector <ll> (n + 1);
        d = vector <ll> (n + 1);
        adj = vector < vector <ll> > (n + 1);
        vis= vector< bool> (n+1);
        For(i,n+1) vis[i]=false;
    }
    void addEdge(ll from, ll to, T cap,bool dir=true) {
        edge<T> e1 = {from, to, cap, (T)0};
        edge<T> e2 = {to, from, (dir)?(T)0:cap, (T)0};
        adj[from].push_back(e.size()); e.push_back(e1);
        adj[to].push_back(e.size()); e.push_back(e2);
    }
    ll bfs() {
        queue <ll> q;
        for(ll i = 0; i <= n; ++i) d[i] = -1;
        q.push(source); d[source] = 0;
        while(!q.empty() and d[sink] < 0) {
            ll x = q.front(); q.pop();
            for(ll i = 0; i < (ll)adj[x].size(); ++i) {
                ll id = adj[x][i], y = e[id].y;
                if(d[y] < 0 and e[id].flow < e[id].cap) {
                    q.push(y); d[y] = d[x] + 1;
                }
            }
        }
        return d[sink] >= 0;
    }
    T dfs(ll x, T flow) {
        if(flow<=eps) return 0;
        if(x == sink) return flow;
        for(;cur[x] < (ll)adj[x].size(); ++cur[x]) {
            ll id = adj[x][cur[x]], y = e[id].y;
            if(d[y] != d[x] + 1) continue;
            T pushed = dfs(y, min(flow, e[id].cap - e[id].flow));
            if(pushed>0) {
                e[id].flow += pushed;
                e[id ^ 1].flow -= pushed;
                return pushed;
            }
        }
        return (T)0;
    }
    T maxFlow(ll src, ll snk) {
        this->source = src; this->sink = snk;
        T flow = 0;
        while(bfs()) {
            for(ll i = 0; i <= n; ++i) cur[i] = 0;
            while(T pushed = dfs(source, inf)) {
                flow += pushed;
                cerr<<pushed;
                break;
            }

            break;
        }
        return flow;
    }
};

map<pair<int,int>,double> capacities;

void solve(){
	long double eps;
	int n,m,s,t;
	cin>>eps>>n>>m>>s>>t;
	DinicFlow<int> U(n+1,0);
	For(i,m){
		int from, to,cap;
		cin>>from>>to>>cap;
        capacities[{from,to}]+=cap;
        capacities[{to,from}]+=cap;
        U.addEdge(from,to,cap,false);
	}
    double f=U.maxFlow(s,t);
    cout<<f<<endl;
    int mn;
    cin>>mn;
    vector<double> Fs(n+1);

    For(i,mn){
        int from,to;
        double flow;
        cin>>from>>to>>flow;
        Fs[from]+=flow;
        Fs[to]-=flow;
        assert((capacities[{from,to}]>=abs(flow)));
    }

    for(int i=1;i<=n;++i){
        if (i!=s && i!=t) assert(Fs[i]==0);
    }
    assert(abs(Fs[s]+Fs[t])<eps);

    double mxf;
    cin>>mxf;

    assert(f*(1-eps)<=mxf && mxf<=f);
    assert(f*(1-eps)<=abs(Fs[s]) && abs(Fs[s])<=f);
    assert(f*(1-eps)<=abs(Fs[t]) && abs(Fs[t])<=f);

}
signed main(){
    // Use "set_name".max_load_factor(0.25);"set_name".reserve(512); with unordered set
    // Uncomment codechef for large input files. Doesn't work on codeforces.
   
    fio;
    int t=1;
    // cin>>t;
    while(t--) {
        solve();
    }
}
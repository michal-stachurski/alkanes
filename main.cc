#include <bits/stdc++.h>
using namespace std;

using ld = long double;
#define EPS (1e-6)
#define MIN_DIST (0.0f)
#define MAX_CNT (0x1ffff)
struct pt {
  ld x, y;
  pt() {}
  pt(ld _x, ld _y) : x(_x), y(_y) {}
  bool operator==(const pt& p) const { return (abs(x-p.x)<EPS && abs(y-p.y)<EPS); }
  bool operator!=(const pt& p) const { return !(*this == p); }
  pt operator+(const pt& p) const { return pt(x + p.x, y + p.y); }
  pt operator-(const pt& p) const { return pt(x - p.x, y - p.y); }
  ld cross(const pt& p) const { return x * p.y - y * p.x; }
  ld cross(const pt& a, const pt& b) const {
    return (a - *this).cross(b - *this);
  }
};

ld dist(pt a, pt b) {
  return sqrt((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y));
}

struct line {
  pt a, b;
};

bool zero(ld x) { return abs(x) < EPS; }

int sgn(ld x) {
  if (zero(x)) return 0;
  else if (x < 0) return -1;
  else return 1;
}

bool leq(ld x, ld y) { return x < y || abs(x-y) < EPS; }

bool inter1(ld a, ld b, ld c, ld d) {
  if (a > b) swap(a, b);
  if (c > d) swap(c, d);
  return leq(max(a, c),min(b, d));
}

bool check_inter(const pt& a, const pt& b, const pt& c, const pt& d) {
  if (zero(c.cross(a, d)) && zero(c.cross(b, d)))
    return inter1(a.x, b.x, c.x, d.x) && inter1(a.y, b.y, c.y, d.y);
  return sgn(a.cross(b, c)) != sgn(a.cross(b, d)) &&
         sgn(c.cross(d, a)) != sgn(c.cross(d, b));
}


int randint(int l, int r) {
  static random_device rd;
  static mt19937 gen(rd());
  return uniform_int_distribution<>(l, r)(gen);
}

constexpr int MAX_R = 0x100000;
vector<pt> roots;
void construct_root() {
  for (ld k = 0.0; k < (ld)MAX_R; k = k + 1.0) {
    ld phi = (2*M_PI_2l*k) / (ld)MAX_R;
    roots.emplace_back(cos(phi), sin(phi));
  }
}
pt randpt(pt p) {
  pt r = roots[randint(0,MAX_R-1)];
  return p+r;
}

bool mid(line x, line y) {
  return 
    x.a == y.a || x.a == y.b || 
    x.b == y.a || x.b == y.b;
}

bool ok(vector<line>& l, line new_line) {
  for (line& z : l) {
    if (!mid(z,new_line) && check_inter(z.a,z.b, new_line.a,new_line.b)) {
      return false;
    }
  }
  return true;
}

bool za_blisko(vector<pt>& p, pt a) {
  for (pt x : p) {
    if (dist(x,a) < MIN_DIST) {
      return true;
    }
  }
  return false;
}

FILE *file;

struct tree {
  vector<vector<int>> adj;
  size_t encoding;
  int n;

  tree() { n = 1; adj = {{}}; encoding = f(); }
  tree(int v, const tree& t) {
    n = t.n+1; adj = t.adj;
    adj.resize(n);
    insert_edge(n-1, v);
    encoding = f();
  }
  void insert_edge(int u, int v) {
    adj[u].push_back(v);
    adj[v].push_back(u);
  }

  string encode(int v, int p = -1) const {
    vector<string> enc;
    for (int u : adj[v]) {
      if (u != p) {
        enc.push_back(encode(u, v));
      }
    }
    if (adj[v].size() == 1) return "10";
    sort(enc.begin(), enc.end());
    string res;
    for (string& e : enc) res += e;
    return "1" + res + "0";
  }
  string encode() const {
    if (n == 1) return "10";
    if (n == 2) return "1100";
    vector<int> c = center();
    if (c.size() == 1) {
      return encode(c.front());
    } else {
      return min(
        encode(c.front()),
        encode(c.back())
      );
    }
    return "2";
  }
  size_t f() const {
    return stoll(encode(), nullptr, 2);
  }

  vector<tree> gen() const {
    vector<tree> res;
    for (int v = 0; v < n; v++) {
      if (adj[v].size() < 4) {
        res.emplace_back(v, *this);
      }
    }
    return res;
  }
  vector<int> center() const {
    vector<int> centroid;
    vector<int> sz(n);
    function<void (int, int)> dfs = [&](int u, int p) {
      sz[u] = 1;
      bool is_centroid = true;
      for (auto v : adj[u])
        if (v != p) {
          dfs(v, u);
          sz[u] += sz[v];
          if (sz[v] > n / 2) is_centroid = false;
        }
      if (n - sz[u] > n / 2) is_centroid = false;
      if (is_centroid) centroid.push_back(u);
    };
    dfs(0, -1);
    return centroid;
  }
  bool operator==(const tree& t) const {
    return encoding == t.encoding;
  }
  bool operator<(const tree& t) const {
    return encoding < t.encoding;
  }

  void read() {
    cin >> n;
    adj.resize(n);
    for (int i = 1; i < n; i++) {
      int u, v;
      cin >> u >> v;
      insert_edge(u, v);
    }
  }
  void print() {
    for (int v = 0; v < n; v++) {
      cout << v << ": ";
      for (int u : adj[v]) {
        cout << u << " ";
      }
      cout << "\n";
    }
  }


  void dfs(vector<line>& l, vector<pt>& cord, int v, int p = -1) const {
    for (int u : adj[v]) {
      if (u != p) {
        pt cand = randpt(cord[v]);
        int cnt = 0;
        while (!ok(l,{cord[v],cand}) || dist(cand,cord[p]) < M_SQRT2 ||
          za_blisko(cord,cand)) {
          if (++cnt >= MAX_CNT) { break; }
          cand = randpt(cord[v]);
        }
        cord[u] = cand;
        l.push_back({cord[u],cord[v]});
        // dfs(l, cord, u, v);
      }
    }

    for (int u : adj[v]) {
      if (u != p){
        dfs(l,cord,u,v);
      }
    }
  }
  void draw() const {
    fprintf(file, "\\begin{tikzpicture}\n");
    int c = center().front();
    vector<line> l;
    vector<pt> cord(n, pt(0,0));
    // cord[c] = pt(0,0);
    dfs(l,cord,c);
    for (int v = 0; v < n; v++) {
      for (int u : adj[v]) {
        if (u < v) continue;
        fprintf(file, "\t\\draw (%Lf,%Lf) -- (%Lf,%Lf);\n", cord[v].x, cord[v].y, cord[u].x, cord[u].y);
      }
    }
    fprintf(file, "\\end{tikzpicture}\n\n");
  }
};

const int LIM = 20;
set<tree> trees[LIM+1];

int main() {

  file = fopen("trees.tex", "w");

  tree one;
  trees[1].insert(one);
  for (int n = 2; n <= LIM; n++) {
    for (const tree& t : trees[n-1]) {
      for (const tree& s : t.gen()) {
        trees[n].insert(s);
      }
    }
    cout << n << ": " << trees[n].size() << '\n';
  }
  construct_root();

  int i = 0;
  for (const tree& t : trees[LIM]) {
    t.draw();
    cerr << ++i << '\r';
  }
  fclose(file);
}
// 109°28'
#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

struct Cell {
    int score = 0;
    int matches = 0;
    int aligned = 0;
    int qcons = 0;
    int tcons = 0;
};

static bool better(const Cell &a, const Cell &b) {
    if (a.score != b.score) return a.score > b.score;
    if (a.matches != b.matches) return a.matches > b.matches;
    if (a.aligned != b.aligned) return a.aligned < b.aligned;
    return (a.qcons + a.tcons) > (b.qcons + b.tcons);
}

static Cell add_diag(Cell c, int delta, bool match) {
    c.score += delta;
    c.aligned += 1;
    c.qcons += 1;
    c.tcons += 1;
    if (match) c.matches += 1;
    if (c.score <= 0) return Cell{};
    return c;
}

static Cell add_gap(Cell c, int delta, bool consume_query) {
    c.score += delta;
    c.aligned += 1;
    if (consume_query) c.qcons += 1;
    else c.tcons += 1;
    if (c.score <= 0) return Cell{};
    return c;
}

static const std::string AA = "ARNDCQEGHILKMFPSTWYV";
static const int B62[20][20] = {
    { 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0},
    {-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3},
    {-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3},
    {-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3},
    { 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1},
    {-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2},
    {-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2},
    { 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3},
    {-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3},
    {-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3},
    {-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1},
    {-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2},
    {-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1},
    {-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1},
    {-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2},
    { 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2},
    { 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0},
    {-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3},
    {-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1},
    { 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4}
};

static int aa_index(char c) {
    auto p = AA.find(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
    return p == std::string::npos ? -1 : static_cast<int>(p);
}

static int subst(char a, char b) {
    int i = aa_index(a), j = aa_index(b);
    if (i < 0 || j < 0) return a == b ? 1 : -1;
    return B62[i][j];
}

static Cell local_align(const std::string &q, const std::string &t) {
    const int gap_open = -10;
    const int gap_extend = -1;
    const size_t n = t.size();
    std::vector<Cell> pm(n + 1), px(n + 1), py(n + 1);
    std::vector<Cell> cm(n + 1), cx(n + 1), cy(n + 1);
    Cell best;
    for (size_t i = 1; i <= q.size(); ++i) {
        cm[0] = cx[0] = cy[0] = Cell{};
        for (size_t j = 1; j <= n; ++j) {
            Cell mbase = pm[j - 1];
            if (better(px[j - 1], mbase)) mbase = px[j - 1];
            if (better(py[j - 1], mbase)) mbase = py[j - 1];
            cm[j] = add_diag(mbase, subst(q[i - 1], t[j - 1]), q[i - 1] == t[j - 1]);

            Cell x1 = add_gap(pm[j], gap_open, true);
            Cell x2 = add_gap(px[j], gap_extend, true);
            Cell x3 = add_gap(py[j], gap_open, true);
            cx[j] = x1;
            if (better(x2, cx[j])) cx[j] = x2;
            if (better(x3, cx[j])) cx[j] = x3;

            Cell y1 = add_gap(cm[j - 1], gap_open, false);
            Cell y2 = add_gap(cy[j - 1], gap_extend, false);
            Cell y3 = add_gap(cx[j - 1], gap_open, false);
            cy[j] = y1;
            if (better(y2, cy[j])) cy[j] = y2;
            if (better(y3, cy[j])) cy[j] = y3;

            if (better(cm[j], best)) best = cm[j];
            if (better(cx[j], best)) best = cx[j];
            if (better(cy[j], best)) best = cy[j];
        }
        pm.swap(cm); px.swap(cx); py.swap(cy);
    }
    return best;
}

static std::vector<std::string> split_tab(const std::string &line) {
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string x;
    while (std::getline(ss, x, '\t')) out.push_back(x);
    return out;
}

int main(int argc, char **argv) {
    if (argc != 3) {
        std::cerr << "usage: stage2_5_local_aligner INPUT.tsv OUTPUT.tsv\n";
        return 2;
    }
    std::ifstream in(argv[1]);
    std::ofstream out(argv[2]);
    if (!in || !out) {
        std::cerr << "failed to open input or output\n";
        return 3;
    }
    std::string line;
    std::getline(in, line);
    out << "pair_id\tseed_accession\tcandidate_accession\talignment_score\tmatches\taligned_columns"
           "\tseed_residues_aligned\tcandidate_residues_aligned\tsequence_identity"
           "\tseed_coverage\tcandidate_coverage\n";
    out << std::fixed << std::setprecision(8);
    while (std::getline(in, line)) {
        auto f = split_tab(line);
        if (f.size() < 5) continue;
        Cell b = local_align(f[3], f[4]);
        double ident = b.aligned ? static_cast<double>(b.matches) / b.aligned : 0.0;
        double qcov = f[3].empty() ? 0.0 : static_cast<double>(b.qcons) / f[3].size();
        double tcov = f[4].empty() ? 0.0 : static_cast<double>(b.tcons) / f[4].size();
        out << f[0] << '\t' << f[1] << '\t' << f[2] << '\t' << b.score << '\t'
            << b.matches << '\t' << b.aligned << '\t' << b.qcons << '\t' << b.tcons << '\t'
            << ident << '\t' << qcov << '\t' << tcov << '\n';
    }
    return 0;
}

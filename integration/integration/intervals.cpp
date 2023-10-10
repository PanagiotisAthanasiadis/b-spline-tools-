#include <iostream>
#include <vector>
#include <set>

using namespace std;

int main() {
    vector<double> v{0.0, 0.0, 0.0, 1.5, 2.0, 3.0, 4.0, 4.0, 5.0, 5.0, 5.0};
    set<pair<double, double>> intervals;

    for (int i = 0; i < v.size() - 1; i++) {
        double start = v[i];
        double end = v[i+1];

        intervals.insert(make_pair(start, end));
    }

    cout << "Unique Intervals:" << endl;

    for (auto interval : intervals) {
        cout << interval.first << "-" << interval.second << endl;
    }

    return 0;
}
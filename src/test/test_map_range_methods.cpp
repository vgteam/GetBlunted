#include "utility.hpp"

#include <iostream>
#include <vector>
#include <map>

using bluntifier::less_than;
using bluntifier::greater_than;
using bluntifier::less_than_or_equal;
using bluntifier::greater_than_or_equal;

using std::cout;
using std::vector;
using std::map;


int main(){
    vector<size_t> keys = {1,2,3,5};

    map<size_t, size_t> m;

    for (auto k: keys){
        m[k] = k;
    }

    for (auto& item: m){
        cout << item.first << " " << item.second << '\n';
    }
    cout << '\n';

    for (auto k: keys){
        cout << "less_than: " << k << '\n';
        vector <map<size_t,size_t>::iterator> result;

        less_than(m, k, result);

        for (auto& r: result){
            cout << r->first << " " << r->second << '\n';
        }
        cout << '\n';
    }

    for (auto k: keys){
        cout << "greater_than: " << k << '\n';
        vector <map<size_t,size_t>::iterator> result;

        greater_than(m, k, result);

        for (auto& r: result){
            cout << r->first << " " << r->second << '\n';
        }
        cout << '\n';
    }


    for (auto k: keys){
        cout << "less_than_or_equal: " << k << '\n';
        vector <map<size_t,size_t>::iterator> result;

        less_than_or_equal(m, k, result);

        for (auto& r: result){
            cout << r->first << " " << r->second << '\n';
        }
        cout << '\n';
    }


    for (auto k: keys){
        cout << "greater_than_or_equal: " << k << '\n';
        vector <map<size_t,size_t>::iterator> result;

        greater_than_or_equal(m, k, result);

        for (auto& r: result){
            cout << r->first << " " << r->second << '\n';
        }
        cout << '\n';
    }

    {
        cout << "less_than: " << 0 << '\n';
        vector <map<size_t,size_t>::iterator> result;

        less_than(m, 0, result);

        for (auto& r: result){
            cout << r->first << " " << r->second << '\n';
        }
        cout << '\n';
    }


    {
        cout << "less_than: " << 4 << '\n';
        vector <map<size_t,size_t>::iterator> result;

        less_than(m, 4, result);

        for (auto& r: result){
            cout << r->first << " " << r->second << '\n';
        }
        cout << '\n';
    }

    {
        cout << "less_than: " << 6 << '\n';
        vector <map<size_t,size_t>::iterator> result;

        less_than(m, 6, result);

        for (auto& r: result){
            cout << r->first << " " << r->second << '\n';
        }
        cout << '\n';
    }

    {
        cout << "greater_than: " << 0 << '\n';
        vector <map<size_t,size_t>::iterator> result;

        greater_than(m, 0, result);

        for (auto& r: result){
            cout << r->first << " " << r->second << '\n';
        }
        cout << '\n';
    }


    {
        cout << "greater_than: " << 4 << '\n';
        vector <map<size_t,size_t>::iterator> result;

        greater_than(m, 4, result);

        for (auto& r: result){
            cout << r->first << " " << r->second << '\n';
        }
        cout << '\n';
    }

    {
        cout << "greater_than: " << 6 << '\n';
        vector <map<size_t,size_t>::iterator> result;

        greater_than(m, 6, result);

        for (auto& r: result){
            cout << r->first << " " << r->second << '\n';
        }
        cout << '\n';
    }


    {
        cout << "less_than_or_equal: " << 0 << '\n';
        vector <map<size_t,size_t>::iterator> result;

        less_than_or_equal(m, 0, result);

        for (auto& r: result){
            cout << r->first << " " << r->second << '\n';
        }
        cout << '\n';
    }


    {
        cout << "less_than_or_equal: " << 4 << '\n';
        vector <map<size_t,size_t>::iterator> result;

        less_than_or_equal(m, 4, result);

        for (auto& r: result){
            cout << r->first << " " << r->second << '\n';
        }
        cout << '\n';
    }

    {
        cout << "less_than_or_equal: " << 6 << '\n';
        vector <map<size_t,size_t>::iterator> result;

        less_than_or_equal(m, 6, result);

        for (auto& r: result){
            cout << r->first << " " << r->second << '\n';
        }
        cout << '\n';
    }


    {
        cout << "greater_than_or_equal: " << 0 << '\n';
        vector <map<size_t,size_t>::iterator> result;

        greater_than_or_equal(m, 0, result);

        for (auto& r: result){
            cout << r->first << " " << r->second << '\n';
        }
        cout << '\n';
    }


    {
        cout << "greater_than_or_equal: " << 4 << '\n';
        vector <map<size_t,size_t>::iterator> result;

        greater_than_or_equal(m, 4, result);

        for (auto& r: result){
            cout << r->first << " " << r->second << '\n';
        }
        cout << '\n';
    }

    {
        cout << "greater_than_or_equal: " << 6 << '\n';
        vector <map<size_t,size_t>::iterator> result;

        greater_than_or_equal(m, 6, result);

        for (auto& r: result){
            cout << r->first << " " << r->second << '\n';
        }
        cout << '\n';
    }

    return 0;
}

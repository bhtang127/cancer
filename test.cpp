#include <iostream>
#include <vector>
#include <random>

using namespace std;

int main(){
    mt19937 eng{random_device{}()};
    vector<double> p(4);
    p[0] = 0.1;
    p[1] = 0.2;
    p[2] = 0.3;
    p[3] = 0.4;
    int choice = discrete_distribution<int>(p.begin(),p.end())(eng);
    cout<<choice<<endl;
    return 0;
}
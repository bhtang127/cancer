#include <iostream>
#include <vector>
#include <random>

using namespace std;

int main(){
    mt19937 eng{random_device{}()};
    vector<double> p(4);
    vector<int> count(4);
    p[0] = 0.1;
    p[1] = 0.2;
    p[2] = 0.3;
    p[3] = 0.4;
    for(int i=0;i<100;i++){
        count [ discrete_distribution<int>(p.begin(),p.end())(eng) ] ++;
    }
    cout<<count[0]<<" "<<count[1]<<" "<<count[2]<<" "<<count[3]<<" "<<endl;
    return 0;
}
#include <iostream>
#include <cstdlib>
#include <set>
#include <vector>

using namespace std;

vector< vector<unsigned int> > sampleID(int cycles, vector<int> n_samples){
    vector< vector<unsigned int> > retval(n_samples.size());
    set<unsigned int> samples;
    vector<int>::iterator it;
    
    int i = 0;
    for(it = n_samples.begin(); it != n_samples.end(); it++, i++){
        if(*it == 1){
            // cout<<"Here for 1"<<endl;
            retval[i].push_back((unsigned int) rand() % (1 << cycles));
        }

        else{
            samples.clear();
            // cout << "Here for 2,3" << endl;
            while(samples.size() < *it){
                samples.insert( (unsigned int) rand() % (1 << cycles) );
            }
            // cout<<"is here?"<<endl;
            retval[i].assign( samples.begin(), samples.end() );
        }
    }
    return retval;
}

unsigned int* translateUID(unsigned int ID, int cycles){
    unsigned int* UIDs;
    unsigned int uid1, uid2;
    if(cycles == 0){
        UIDs = (unsigned int*) malloc(2 * sizeof(unsigned int));
        UIDs[0] = 0;
        UIDs[1] = 0;
        return UIDs;
    }
    else if(ID % 2 == 1){
        uid1 = translateUID((unsigned int) (ID/2), cycles-1)[1];
        UIDs = (unsigned int*) malloc(2 * sizeof(unsigned int));
        UIDs[0] = uid1;
        UIDs[1] = uid1;
        return UIDs;
    }
    else{
        uid1 = translateUID(ID/2, cycles-1)[0];
        uid2 = ID / 2 + 1 << (cycles-1);
        UIDs = (unsigned int*) malloc(2 * sizeof(unsigned int));
        UIDs[0] = uid1;
        UIDs[1] = uid2;
        return UIDs;
    }

} 

unsigned int* translateWBC(unsigned int ID, int cycles){
    unsigned int* WBCs;
    unsigned int w1, w2;
    if(ID == 0){
        WBCs = (unsigned int*) malloc(2 * sizeof(unsigned int));
        WBCs[0] = 0;
        WBCs[1] = 1;
        return WBCs;
    }
    else if(ID < 1<<(cycles-1)){
        WBCs = (unsigned int*) malloc(2 * sizeof(unsigned int));
        WBCs[0] = 1;
        WBCs[1] = 1;
        return WBCs;
    }
    else{
        return translateWBC( ID - 1<<(cycles-1), cycles-1 );
    }
} 

vector<unsigned int>* translate(vector< vector<unsigned int> > IDs, int cycles, const char tag){
    vector<unsigned int>* translation;
    translation = (vector<unsigned int>*) malloc(2 * sizeof(vector<unsigned int>));
    vector<unsigned int>::iterator it;
    unsigned int *ids;
    if(tag == 'U') {
        for(int i=0; i < IDs.size(); i++){
            for(it=IDs[i].begin(); it!=IDs[i].end(); it++){
                ids = translateUID(*it, cycles);
                translation[0].push_back(ids[0]);
                translation[1].push_back(ids[1]);
            }
        }
    }
    else {
        for(int i=0; i < IDs.size(); i++){
            for(it=IDs[i].begin(); it!=IDs[i].end(); it++){
                ids = translateWBC(*it, cycles);
                translation[0].push_back(ids[0]);
                translation[1].push_back(ids[1]);
            }
        }
    }
    return translation;
}

vector< vector<unsigned int> >

vector<int> expand(vector<int> n_samples){
    vector<int> indexes;
    for(int i=0; i < n_samples.size(); ++i){
        if(n_samples[i] > 0){
            for(int j=0; j < n_samples[i]; ++j){
            indexes.push_back(i);
            }
        }
    }
    return indexes;
}



int main(){
    int c = 30, iter = 15655642;
    vector< vector<unsigned int> > l;
    vector<int> ns;
    for(int i=0; i < iter; i++){
        ns.push_back((int) rand()%2);
    }
    l = sampleID(c, ns);
    translate(l, c, 'U');
    translate(l, c, 'W');    
    expand(ns);
    cout << l.size() << endl;
    return 0;
}
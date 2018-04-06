#include <iostream>
#include <cstdlib>
#include <set>
#include <vector>

using namespace std;

vector<unsigned int> sampleID(int cycles, vector<int> n_samples){
    vector<unsigned int> retval(n_samples);
    set<unsigned int> samples;
    
    if(n_samples == 1){
        retval.push_back((unsigned int) rand() % (1 << cycles));
        return retval;
    }

    else{
        while(samples.size() < n_samples){
            samples.insert( (unsigned int) rand() % (1 << cycles) );
        }
        copy( samples.begin(), samples.end(), retval.begin() );
        return retval;
    }
}

def translate_aux(id, cycles):
    if cycles == 0:
        return [0, 0]
    if id % 2 == 1:
        uid = translate_aux( int(id / 2), cycles - 1 )[1]
        return [uid, uid]
    else:
        uid1 = translate_aux( int(id / 2), cycles - 1 )[0]
        uid2 = id / 2 + 2 ** (cycles - 1)
        return [uid1, uid2]

def translate(samples_id, cycles):
    if samples_id is None:
        return None
    n = len(samples_id)
    UID = np.array([ translate_aux(id, cycles) for id in samples_id ], dtype = np.int32)
    return UID

vector<unsigned int> 

vector<vector<unsigned int>> translateUID(vector<unsigned int> IDs, int cycles, const char tag){
    vector<vector<unsigned int>> translation;

}

int main(){
    int c = 30, n_samples = 2, iter = 15655642;
    vector<unsigned int> l;
    for(int i=0; i < iter; i++){
        l = sampleID(c, n_samples);
    }
    cout << l.size() << endl;
    return 0;
}
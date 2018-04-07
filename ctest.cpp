#include <iostream>
#include <cstdlib>
#include <set>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iterator>

using namespace std;

vector<unsigned int> sampleID(int cycles, int n_samples){

    vector<unsigned int> retval(n_samples);
    set<unsigned int> samples;
    
    if(n_samples == 1){
        // cout<<"Here for 1"<<endl;
        retval[0] = ((unsigned int) rand() % (1 << cycles));
    }
    else{
        samples.clear();
        // cout << "Here for 2,3" << endl;
        while(samples.size() < n_samples){
            samples.insert( (unsigned int) rand() % (1 << cycles) );
        }
        // cout<<"is here?"<<endl;
        retval.assign( samples.begin(), samples.end() );
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

unsigned int factorial(unsigned int n){
 	unsigned int retval = 1;
 	for (int i = n; i > 1; --i)
 		retval *= i;
 	return retval;
}

double poisson_pmf(int k, double mu){
    if(k == 0) return exp(-mu);
    return exp(-mu) * pow(mu, k) / factorial(k);
}

int random_choice(vector<long double> p_val){
    default_random_engine generator;
    uniform_real_distribution<long double> distribution(0.0,1.0);
    int choice = 0;
    long double cumm, p = distribution(generator);
    distribution.reset();
    for(cumm=p_val[0]; cumm < p; choice++){
        cumm += p_val[choice];
    }
    return choice;
}

vector< set<unsigned int> > find_all_ancestors(vector<unsigned int> ids, int cycles){
    vector< set<unsigned int> > ancestors(cycles);
    for(int i=0; i<cycles; i++){
        for(int j=0; j<ids.size(); j++){
            ancestors[i].insert((unsigned int) (ids[j] / (1<<i)));
        }
    }
    return ancestors;
}

vector< vector<unsigned int> > cpp_poisson_mutation(vector<int> nsamples_per_molecule, int cycles, 
                                                int bases_per_amplicon, double error_rate, const char tag){

    int n_molecule = nsamples_per_molecule.size();
    unsigned int aid1, aid2, *ID, ns_max;
    ns_max = *max_element(nsamples_per_molecule.begin(), nsamples_per_molecule.end());
    double control_mean = ns_max * (cycles-floor(log2(ns_max))+2) * bases_per_amplicon * error_rate;
    double true_mean;
    int cut_off = 0, ns, randnum, choice, aux, *mutation_count, layer, index;
    double zero_control_mean;
    long double psum;
    vector<long double> control_pval(2), choice_pval, poisson_pval, adjust_pval;
    vector<unsigned int> layers, aids, state;
    vector<unsigned int>::iterator it;
    set<unsigned int>::iterator st;

    int sum=0, ancsum=0;
    for(int i=0; i < n_molecule; i++) sum += nsamples_per_molecule[i];
    vector< vector<unsigned int> > retval(sum), mutation_loc, mutation_loc_aux;
    vector<unsigned int> ids;
    vector< set<unsigned int> > ancestors;

    while( n_molecule * poisson_pmf(cut_off, control_mean) > 1e-5 ) ++cut_off;
    cut_off = max(cut_off, 2);
    cout<<"cut_off: "<<cut_off<<" iteration: "<<n_molecule<<endl;
    
    if(tag == 'U'){
        for(int i=0,loc=0; i<n_molecule; i++){
            ns = nsamples_per_molecule[i];
            if(i % (n_molecule/10) == 0){
                cout<<"Processing: "<<i*100.0/n_molecule<<"%"<<endl;
            }
            zero_control_mean = ns * (cycles-floor(log2(ns))+2) * bases_per_amplicon * error_rate;
            control_pval[0] = poisson_pmf(0, zero_control_mean); 
            control_pval[1] = 1 - control_pval[0];
            if(random_choice(control_pval) == 0){
                ids = sampleID(cycles, ns);
                for(int j=0; j<ids.size(); j++){
                    ID = translateUID( ids[j], cycles );
                    retval[loc].push_back(ID[0]);
                    retval[loc].push_back(ID[1]);
                    retval[loc].push_back(0);
                    retval[loc].push_back(0);
                    loc++;       
                }
            }
            else{
                ids = sampleID(cycles, ns);
                ancestors = find_all_ancestors(ids, cycles);
                ancsum = 0;
                for(int k=0; k < ancestors.size(); k++){
                    ancsum += ancestors[k].size();
                }
                true_mean = ancsum * bases_per_amplicon * error_rate;

                psum = 0.0;
                poisson_pval.clear();
                for(int k=0; k < cut_off; k++){
                    poisson_pval.push_back(poisson_pmf(k, true_mean));
                    psum += poisson_pval[k];
                }
                poisson_pval.push_back(1 - psum);

                adjust_pval.clear();
                adjust_pval.push_back((poisson_pval[0] - control_pval[0]) / control_pval[1]);
                for(int k=1; k < poisson_pval.size(); k++){
                    adjust_pval.push_back(poisson_pval[k] / control_pval[1]);
                }

                randnum = random_choice(adjust_pval);
                if(randnum == 0){
                    for(int j=0; j<ids.size(); j++){
                        ID = translateUID( ids[j], cycles );
                        retval[loc].push_back(ID[0]);
                        retval[loc].push_back(ID[1]);
                        retval[loc].push_back(0);
                        retval[loc].push_back(0);
                        loc++;       
                    }
                }
                else{
                    choice_pval.clear();
                    for(int k=0; k < ancestors.size(); k++){
                        choice_pval.push_back(ancestors[k].size() / (long double) ancsum);
                    }
                    layers.clear();
                    aids.clear();
                    for(int k=0; k < randnum; k++){
                        choice = random_choice(choice_pval);
                        layers.push_back(choice);
                        aux = rand() % ancestors[choice].size();
                        st = ancestors[choice].begin();
                        for(int l=0; l<aux; l++) st++;
                        aids.push_back(*st);
                    }
                    
                    mutation_count = (int *)malloc(ns*2 * sizeof(int));
                    for(int k=0; k < 2*ns; k++) mutation_count[k] = 0;

                    mutation_loc.clear();
                    mutation_loc.resize(randnum * ns);
                    mutation_loc_aux.clear();
                    mutation_loc_aux.resize(randnum * ns);
                    for(int k=0; k < randnum; k++){
                        mutation_loc[k].push_back(aids[k]);
                        mutation_loc[k].push_back(layers[k]);
                        mutation_loc[k].push_back(1);
                    }
                    for(int k=0; k < cycles; k++){
                        mutation_loc_aux.clear();
                        mutation_loc_aux.resize(randnum * ns);
                        aux = 0;
                        for(int m=0; m < mutation_loc.size(); m++){
                            if(mutation_loc[m].size() == 0) continue;
                            if(mutation_loc[m][1] == 0) {
                                mutation_loc_aux[aux].push_back(mutation_loc[m][0]);
                                mutation_loc_aux[aux].push_back(mutation_loc[m][1]);
                                mutation_loc_aux[aux].push_back(mutation_loc[m][2]);
                                aux++;
                            }
                            else if(mutation_loc[m][2] == 1){
                                aid1 = mutation_loc[m][0] * 2 + mutation_loc[m][0] % 2 == 0;
                                layer = mutation_loc[m][1];
                                if(ancestors[layer-1].count(aid1) == 0) continue;
                                else{
                                    mutation_loc_aux[aux].push_back(aid1);
                                    mutation_loc_aux[aux].push_back(mutation_loc[m][1]-1);
                                    mutation_loc_aux[aux].push_back(2);
                                    aux++;
                                }
                            }
                            else{
                                layer = mutation_loc[m][1];
                                aid1 = mutation_loc[m][1] * 2;
                                aid2 = mutation_loc[m][1] * 2 + 1;
                                if(ancestors[layer-1].count(aid1) > 0){
                                    mutation_loc_aux[aux].push_back(aid1);
                                    mutation_loc_aux[aux].push_back(mutation_loc[m][1]-1);
                                    mutation_loc_aux[aux].push_back(mutation_loc[m][2]);
                                    aux++;
                                }
                                if(ancestors[layer-1].count(aid2) > 0){
                                    mutation_loc_aux[aux].push_back(aid2);
                                    mutation_loc_aux[aux].push_back(mutation_loc[m][1]-1);
                                    mutation_loc_aux[aux].push_back(mutation_loc[m][2]);
                                    aux++;
                                }
                            }
                        }
                        mutation_loc = mutation_loc_aux;
                    }
                    for(int k=0; k < mutation_loc.size(); k++){
                        if(mutation_loc[k].size() == 0) continue;
                        else{
                            if(mutation_loc[k][2] == 1){
                                index = distance(ids.begin(), find(ids.begin(), ids.end(), mutation_loc[k][0]));
                                mutation_count[2 * index + (mutation_loc[k][0] % 2 == 0)] += 1;
                            }
                            else{
                                index = distance(ids.begin(), find(ids.begin(), ids.end(), mutation_loc[k][0]));
                                mutation_count[2 * index] += 1;
                                mutation_count[2 * index + 1] += 1;                                                                
                            }
                        }
                    }
                    for(int k=0; k < ns; k++){
                        ID = translateUID( ids[k], cycles );
                        retval[loc].push_back(ID[0]);
                        retval[loc].push_back(ID[1]);
                        retval[loc].push_back(mutation_count[2 * k]);
                        retval[loc].push_back(mutation_count[2 * k + 1]);
                        loc++;
                    }
                    free(mutation_count);
                }
            }
        }
    }
    else{
        for(int i=0,loc=0; i < n_molecule; i++){
            ns = nsamples_per_molecule[i];
            if(i % (n_molecule/10) == 0){
                cout<<"Processing: "<<i*100.0/n_molecule<<"%"<<endl;
            }
            zero_control_mean = ns * (cycles-floor(log2(ns))+2) * bases_per_amplicon * error_rate;
            control_pval[0] = poisson_pmf(0, zero_control_mean); 
            control_pval[1] = 1 - control_pval[0];
            if(random_choice(control_pval) == 0){
                ids = sampleID(cycles, ns);
                for(int j=0; j < ids.size(); j++){
                    ID = translateWBC( ids[j], cycles );
                    retval[loc].push_back(ID[0]);
                    retval[loc].push_back(ID[1]);
                    retval[loc].push_back( ids[j] < (1<<(cycles-1)) );
                    retval[loc].push_back(0);
                    retval[loc].push_back(0);
                    loc++;       
                }
            }
            else{
                ids = sampleID(cycles, ns);
                ancestors = find_all_ancestors(ids, cycles);
                ancsum = 0;
                for(int k=0; k < ancestors.size(); k++){
                    ancsum += ancestors[k].size();
                }
                true_mean = ancsum * bases_per_amplicon * error_rate;

                psum = 0.0;
                poisson_pval.clear();
                for(int k=0; k < cut_off; k++){
                    poisson_pval.push_back(poisson_pmf(k, true_mean));
                    psum += poisson_pval[k];
                }
                poisson_pval.push_back(1 - psum);

                adjust_pval.clear();
                adjust_pval.push_back((poisson_pval[0] - control_pval[0]) / control_pval[1]);
                for(int k=1; k < poisson_pval.size(); k++){
                    adjust_pval.push_back(poisson_pval[k] / control_pval[1]);
                }

                randnum = random_choice(adjust_pval);
                if(randnum == 0){
                    for(int j=0; j<ids.size(); j++){
                        ID = translateWBC( ids[j], cycles );
                        retval[loc].push_back(ID[0]);
                        retval[loc].push_back(ID[1]);
                        retval[loc].push_back( ids[j] < (1<<(cycles-1)) );                        
                        retval[loc].push_back(0);
                        retval[loc].push_back(0);
                        loc++;       
                    }
                }
                else{
                    choice_pval.clear();
                    for(int k=0; k < ancestors.size(); k++){
                        choice_pval.push_back(ancestors[k].size() / (long double) ancsum);
                    }
                    layers.clear();
                    aids.clear();
                    for(int k=0; k < randnum; k++){
                        choice = random_choice(choice_pval);
                        layers.push_back(choice);
                        aux = rand() % ancestors[choice].size();
                        st = ancestors[choice].begin();
                        for(int l=0; l<aux; l++) st++;
                        aids.push_back(*st);
                    }
                    
                    mutation_count = (int *)malloc(ns*2 * sizeof(int));
                    for(int k=0; k < 2*ns; k++) mutation_count[k] = 0;

                    mutation_loc.clear();
                    mutation_loc.resize(randnum * ns);
                    mutation_loc_aux.clear();
                    mutation_loc_aux.resize(randnum * ns);
                    for(int k=0; k < randnum; k++){
                        mutation_loc[k].push_back(aids[k]);
                        mutation_loc[k].push_back(layers[k]);
                        mutation_loc[k].push_back(1);
                    }
                    for(int k=0; k < cycles; k++){
                        mutation_loc_aux.clear();
                        mutation_loc_aux.resize(randnum * ns);
                        aux = 0;
                        for(int m=0; m < mutation_loc.size(); m++){
                            if(mutation_loc[m].size() == 0) continue;
                            if(mutation_loc[m][1] == 0) {
                                mutation_loc_aux[aux].push_back(mutation_loc[m][0]);
                                mutation_loc_aux[aux].push_back(mutation_loc[m][1]);
                                mutation_loc_aux[aux].push_back(mutation_loc[m][2]);
                                aux++;
                            }
                            else if(mutation_loc[m][2] == 1){
                                aid1 = mutation_loc[m][0] * 2 + mutation_loc[m][0] % 2 == 0;
                                layer = mutation_loc[m][1];
                                if(ancestors[layer-1].count(aid1) == 0) continue;
                                else{
                                    mutation_loc_aux[aux].push_back(aid1);
                                    mutation_loc_aux[aux].push_back(mutation_loc[m][1]-1);
                                    mutation_loc_aux[aux].push_back(2);
                                    aux++;
                                }
                            }
                            else{
                                layer = mutation_loc[m][1];
                                aid1 = mutation_loc[m][1] * 2;
                                aid2 = mutation_loc[m][1] * 2 + 1;
                                if(ancestors[layer-1].count(aid1) > 0){
                                    mutation_loc_aux[aux].push_back(aid1);
                                    mutation_loc_aux[aux].push_back(mutation_loc[m][1]-1);
                                    mutation_loc_aux[aux].push_back(mutation_loc[m][2]);
                                    aux++;
                                }
                                if(ancestors[layer-1].count(aid2) > 0){
                                    mutation_loc_aux[aux].push_back(aid2);
                                    mutation_loc_aux[aux].push_back(mutation_loc[m][1]-1);
                                    mutation_loc_aux[aux].push_back(mutation_loc[m][2]);
                                    aux++;
                                }
                            }
                        }
                        mutation_loc = mutation_loc_aux;
                    }
                    cout<<"Aggregate"<<endl;
                    for(int k=0; k < mutation_loc.size(); k++){
                        if(mutation_loc[k].size() == 0) continue;
                        else{
                            if(mutation_loc[k][2] == 1){
                                index = distance(ids.begin(), find(ids.begin(), ids.end(), mutation_loc[k][0]));
                                mutation_count[2 * index + (mutation_loc[k][0] % 2 == 0)] += 1;
                            }
                            else{
                                index = distance(ids.begin(), find(ids.begin(), ids.end(), mutation_loc[k][0]));
                                mutation_count[2 * index] += 1;
                                mutation_count[2 * index + 1] += 1;                                                                
                            }
                        }
                    }
                    for(int k=0; k < ns; k++){
                        ID = translateWBC( ids[k], cycles );
                        retval[loc].push_back(ID[0]);
                        retval[loc].push_back(ID[1]);
                        retval[loc].push_back( ids[k] < (1<<(cycles-1)) );                                                
                        retval[loc].push_back(mutation_count[2 * k]);
                        retval[loc].push_back(mutation_count[2 * k + 1]);
                        loc++;
                    }
                    free(mutation_count);
                }
            }
        }
    }

    return retval;
}

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
    int c = 30, iter = 1000000, s=1;
    vector< vector<unsigned int> > l;
    vector<int> ns;
    vector<unsigned int> ws(2);
    set<int> ss;
    set<int>::iterator st;
    vector< vector<unsigned int> > res;
    for(int i=0; i < iter; i++){
        ns.push_back((int) 1 + rand()%3);
    }
    res = cpp_poisson_mutation(ns,30,33,1e-6,'U');
    cout<< res.size() << " " << res[0].size() << endl;
    return 0;
}

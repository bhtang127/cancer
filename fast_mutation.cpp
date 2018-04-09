#include "fast_mutation.h"
#include <iostream>
#include <ctime>

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

double poisson_pmf(int k, double lambda){
    return pow(M_E, k * log(lambda) - lambda - lgamma(k + 1.0));
}

int random_choice(std::vector<long double> p_val, long double r){
    long double cumm = p_val[0];
    int i=0;
    for(; cumm < r;)
        cumm += p_val[++i];
    return i;
}

std::vector< std::set<unsigned int> > find_all_ancestors(std::vector<unsigned int> ids, int cycles){
    std::vector< std::set<unsigned int> > ancestors(cycles);
    for(int i=0; i<cycles; i++){
        for(int j=0; j<ids.size(); j++){
            ancestors[i].insert((unsigned int) (ids[j] / (1<<i)));
        }
    }
    return ancestors;
}

void inplace_sampling(int* samples, int len, int total_samples, int restriction){
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<int> dist(0, len-1);
    int cumm = 0, randi;
    while(cumm < total_samples){
        randi = dist(mt);
        if(restriction != 0 && samples[randi] > restriction) continue;
        else{
            samples[randi]++;
            cumm++;
        }
    }
    return;
}

using uint32 = unsigned int;

class Random {
public:
    Random() = default;
    Random(std::mt19937::result_type seed) : eng(seed) {}
    uint32 DrawInt(uint32 Max);
    int Choice(std::vector<long double> p_val);
    long double DrawDouble(double max);

private:        
    std::mt19937 eng{std::random_device{}()};
};

uint32 Random::DrawInt(uint32 Max)
{
    return std::uniform_int_distribution<uint32>{0, Max - 1}(eng);
}

int Random::Choice(std::vector<long double> p_val){
    return std::discrete_distribution<int>{p_val.begin(), p_val.end()}(eng);
}

long double Random::DrawDouble(double max){
    return   std::uniform_real_distribution<long double>{0.0, max}(eng);
}

int find_index(int cycles, std::set<uint32> indexset, uint32 randi){
    uint32 index = 0, count = 0;
    for(; index < 1<<cycles; index++){
        if(indexset.count(index) != 0) count++;
        if(count == (randi+1)) return index; 
    }
}

std::vector<uint32> sampleID(int cycles, int n_samples, Random rd){
    std::vector<uint32> retval(n_samples);
    std::set<uint32> maintain;
    uint32 randi, ind;
    if(n_samples == 1){
        retval[0] = rd.DrawInt(1<<cycles);
    }
    else if(cycles < 15){
        std::set<uint32> indexes;
        for(int j=0; j < n_samples; j++){
            randi = rd.DrawInt(1<<cycles - j);
            ind = find_index(cycles, indexes, randi);
            indexes.insert(ind);
            retval[j] = ind;
        }
        std::sort(retval.begin(), retval.end());
    }
    else{
        while(maintain.size() < n_samples){
            maintain.insert( rd.DrawInt(1<<cycles) );
        }
        retval.assign( maintain.begin(), maintain.end() );
    }
    return retval;
}

void inplace_poisson_mutation(unsigned int* retval, int rows, int columns,
                              unsigned int* nsamples_per_molecule, int n_molecule, int cycles, 
                              int bases_per_amplicon, double error_rate){


    unsigned int aid1, aid2, *ID, ns_max = 0;
    for(int i=0; i < n_molecule; ++i){
        if(nsamples_per_molecule[i] > ns_max){
            ns_max = nsamples_per_molecule[i];
        }
    }

    int total_mutation_counts = 0;
    double control_mean = ns_max * (cycles-floor(log2(ns_max))+2) * bases_per_amplicon * error_rate;
    double true_mean;
    int cut_off = 0, ns, randnum, choice, aux, *mutation_count, layer, index;
    double zero_control_mean;
    long double psum;
    std::vector<long double> control_pval(2), choice_pval, poisson_pval, adjust_pval;
    std::vector<unsigned int> layers, aids, state;
    std::set<unsigned int> idset;
    std::vector<unsigned int>::iterator it;
    std::set<unsigned int>::iterator st;

    std::cout<<"New Version 2: "<<std::endl;
    std::cout<<"control_mean: "<<control_mean<<std::endl;

    int sum=0, ancsum=0;
    for(int i=0; i < n_molecule; i++) sum += nsamples_per_molecule[i];
    std::vector< std::vector<unsigned int> > mutation_loc, mutation_loc_aux;
    std::vector<unsigned int> ids;
    std::vector< std::set<unsigned int> > ancestors;
    
    while( n_molecule * poisson_pmf(cut_off, control_mean) > 1e-5 ) ++cut_off;
    cut_off = std::max(cut_off, 2);
    std::cout<<"cut_off: "<<cut_off<<" iteration: "<<n_molecule<<std::endl;
    
    Random rd;
    srand( time(NULL) );

    if(columns == 4){
        for(int i=0,loc=0; i<n_molecule; i++){
            ns = nsamples_per_molecule[i];
            if(i % (n_molecule/10) == 0){
                std::cout<<"Processing: "<<i*100.0/n_molecule<<"%"<<std::endl;
            }
            zero_control_mean = ns * (cycles-floor(log2(ns))+2) * bases_per_amplicon * error_rate;
            control_pval[0] = poisson_pmf(0, zero_control_mean); 
            control_pval[1] = 1 - control_pval[0];
            if(rd.Choice(control_pval) == 0){
                ids = sampleID(cycles, ns, rd);
                for(int j=0; j<ids.size(); j++){
                    ID = translateUID( ids[j], cycles );
                    retval[columns * loc] = ID[0];
                    retval[columns * loc + 1] = ID[1];
                    loc++;       
                }
            }
            else{
                ids = sampleID(cycles, ns, rd);
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

                // for(int k=0; k < adjust_pval.size(); k++){
                //     std::cout << adjust_pval[k] << " ";
                // }
                // std::cout << std::endl;

                randnum = random_choice(adjust_pval, rd.DrawDouble(1));
                if(randnum == 0){
                    for(int j=0; j<ids.size(); j++){
                        ID = translateUID( ids[j], cycles );
                        retval[columns * loc] = ID[0];
                        retval[columns * loc + 1] = ID[1];
                        loc++;       
                    }
                }
                else{
                    total_mutation_counts += randnum;
                    // std::cout<<i<<" ";
                    // if(randnum > 1) std::cout<<"here a multiple mutations of "<<randnum<<" at iteration "<<i<<std::endl;
                    choice_pval.clear();
                    for(int k=0; k < ancestors.size(); k++){
                        choice_pval.push_back(ancestors[k].size() / (long double) ancsum);
                    }
                    layers.clear();
                    aids.clear();
                    for(int k=0; k < randnum; k++){
                        choice = random_choice(choice_pval, rd.DrawDouble(1));
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
                                aid1 = mutation_loc[m][0] * 2 + (mutation_loc[m][0] % 2 == 0);
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
                                index = distance(ids.begin(), std::find(ids.begin(), ids.end(), mutation_loc[k][0]));
                                mutation_count[2 * index + (mutation_loc[k][0] % 2 == 0)] += 1;
                            }
                            else{
                                index = distance(ids.begin(), std::find(ids.begin(), ids.end(), mutation_loc[k][0]));
                                mutation_count[2 * index] += 1;
                                mutation_count[2 * index + 1] += 1;                                                                
                            }
                        }
                    }
                    for(int k=0; k < ns; k++){
                        ID = translateUID( ids[k], cycles );
                        retval[columns * loc] = ID[0];
                        retval[columns * loc + 1] = ID[1];
                        retval[columns * loc + 2] += mutation_count[2 * k];
                        retval[columns * loc + 3] += mutation_count[2 * k + 1];
                        loc++;
                    }
                    free(mutation_count);
                }
            }
        }
    }
    else{
        for(int i=0,loc=0; i < n_molecule; i++){
            // std::cout<<"Start Here"<<std::endl;            
            ns = nsamples_per_molecule[i];
            if(i % (n_molecule/10) == 0){
                std::cout<<"Processing: "<<i*100.0/n_molecule<<"%"<<std::endl;
            }
            zero_control_mean = ns * (cycles-floor(log2(ns))+2) * bases_per_amplicon * error_rate;
            control_pval[0] = poisson_pmf(0, zero_control_mean); 
            control_pval[1] = 1 - control_pval[0];
            // std::cout<<"Here"<<std::endl;
            if(rd.Choice(control_pval) == 0){
                // std::cout<<ns<<" "<<i<<std::endl;             
                ids = sampleID(cycles, ns, rd);
                // std::cout<<"pass here"<<std::endl;                             
                for(int j=0; j < ids.size(); j++){
                    ID = translateWBC( ids[j], cycles );
                    retval[columns * loc] = ID[0];
                    retval[columns * loc + 1] = ID[1];
                    retval[columns * loc + 2] = ( ids[j] < (1<<(cycles-1)) );
                    loc++;       
                }
            }            
            else{
                // std::cout<<"Else Here"<<std::endl;            
                ids = sampleID(cycles, ns, rd);
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

                randnum = random_choice(adjust_pval, rd.DrawDouble(1));
                if(randnum == 0){
                    for(int j=0; j<ids.size(); j++){
                        ID = translateWBC( ids[j], cycles );
                        retval[columns * loc] = ID[0];
                        retval[columns * loc + 1] = ID[1];
                        retval[columns * loc + 2] = ( ids[j] < (1<<(cycles-1)) );                       
                        loc++;       
                    }
                }
                else{
                    total_mutation_counts += randnum;

                    choice_pval.clear();
                    for(int k=0; k < ancestors.size(); k++){
                        choice_pval.push_back(ancestors[k].size() / (long double) ancsum);
                    }
                    layers.clear();
                    aids.clear();
                    for(int k=0; k < randnum; k++){
                        choice = random_choice(choice_pval, rd.DrawDouble(1));
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
                                aid1 = mutation_loc[m][0] * 2 + (mutation_loc[m][0] % 2 == 0);
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
                                index = distance(ids.begin(), std::find(ids.begin(), ids.end(), mutation_loc[k][0]));
                                mutation_count[2 * index + (mutation_loc[k][0] % 2 == 0)] += 1;
                            }
                            else{
                                index = distance(ids.begin(), std::find(ids.begin(), ids.end(), mutation_loc[k][0]));
                                mutation_count[2 * index] += 1;
                                mutation_count[2 * index + 1] += 1;                                                                
                            }
                        }
                    }
                    for(int k=0; k < ns; k++){
                        ID = translateWBC( ids[k], cycles );
                        retval[columns * loc] = ID[0];
                        retval[columns * loc + 1] = ID[1];
                        retval[columns * loc + 2] = ( ids[k] < (1<<(cycles-1)) );                                                
                        retval[columns * loc + 3] += mutation_count[2 * k];
                        retval[columns * loc + 4] += mutation_count[2 * k + 1];
                        loc++;
                    }
                    free(mutation_count);
                }
            }
        }
    }
    std::cout<<"total mutation count: "<<total_mutation_counts<<std::endl;
    return;
}

void inplace_expand(int* indexes, int olen, int* n_samples, int ilen){
    for(int i=0, j=0; i < ilen; ++i){
        if(n_samples[i] > 0){
            for(int k=0; k < n_samples[i]; ++k, ++j){
                indexes[j] = i;
            }
        }
    }
    return;
}


int main(){
    int c = 30, iter = 10000, iter2=20000000, s=1;
    std::vector< std::vector<unsigned int> > l;
    std::vector<unsigned int> ns, ns2;
    std::vector<unsigned int> ws(2);
    std::set<int> ss;
    std::set<int>::iterator st;
    std::vector< std::vector<unsigned int> > res;
    int sum = 0, sum1 = 0,  r;
    unsigned int *retval;
    for(int i=0; i < iter2; i++){
        r = (int) 1 + rand()%3;
        ns.push_back(r);
        sum += r;
    }
    for(int i=0; i < iter; i++){
        ns2.push_back(2000);
        sum1 += 2000;
    }
    retval = (unsigned int*)malloc( sum1 * 4 * sizeof(unsigned int) );
    inplace_poisson_mutation(retval, sum1, 4, &ns2[0], iter, 30,33,1e-6);
    retval = (unsigned int*)malloc( sum * 5 * sizeof(unsigned int) );
    inplace_poisson_mutation(retval, sum, 5, &ns[0], iter2, 5,33,1e-6);

    // ns.push_back(0);
    // ns.push_back(1);
    // ns.push_back(2);
    // ns.push_back(3);
    // ns.erase(ns.begin()+2);
    // std::cout<< ns[0] <<" "<<ns[1]<<" "<<ns[2]<<" "<<std::endl;
    return 0;
}

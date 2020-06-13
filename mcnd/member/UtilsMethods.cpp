//
//  UtilsMethods.cpp
//  
//
//  Created by Rui Shibasaki on 30/11/2018.
//
//

#include "UtilsMethods.hpp"
#include <bitset>

//===================================================================================
//===================================================================================
// Read/ Translate methods
//===================================================================================
//===================================================================================


void MCND_read_data(std::string fname, Data & data) {
    
    std::ifstream file;
    file.open(fname.c_str());
    if (!file.is_open()) {
        std::cout<<"Failure to open ufl datafile: %s\n "<<fname;
        abort();
    }
    
    int& ndemands = data.ndemands;
    int& narcs = data.narcs;
    int& nnodes = data.nnodes;
    
    
    std::string s;
    std::istringstream ss;
    int style=0;
    getline(file,s);
    if(s == "MULTIGEN.DAT:" ||s == " MULTIGEN.DAT:"){
        style=1;
        getline(file,s);
    }
    ss.str(s);
    
    // read number of locations and number of customers
    
    ss>>nnodes;
    ss>>narcs;
    ss>>ndemands;
    
    ss.clear();
    
    data.d_k.resize(ndemands);
    data.arcs.resize(narcs);
    data.grid.resize(nnodes*nnodes,-1);
    data.kgrid.resize(nnodes*nnodes,0.0);
    if(style == 1){
        double cost=0;
        for(int i=0; i<narcs; ++i){
            getline(file,s);
            ss.str(s);
            ss>>data.arcs[i].i;
            ss>>data.arcs[i].j;
            ss>>cost;
            ss>>data.arcs[i].capa;
            ss>>data.arcs[i].f;
            data.grid[(data.arcs[i].i-1)*nnodes+(data.arcs[i].j-1)]=i;
            ss.clear();
            data.arcs[i].c.resize(ndemands, cost);
            data.arcs[i].b.resize(ndemands, data.arcs[i].capa);
        }
        
        for(int i=0; i<ndemands; ++i){
            getline(file,s);
            ss.str(s);
            ss>> data.d_k[i].O >> data.d_k[i].D >> data.d_k[i].quantity;
            data.kgrid[(data.d_k[i].O-1)*nnodes+(data.d_k[i].D-1)] += data.d_k[i].quantity ;

            ss.clear();
        }
        //calculate variable bounds
        for(int k=0;k<ndemands;++k){
            for(int e=0; e<narcs; ++e){
                if(data.d_k[k].quantity < data.arcs[e].b[k])
                data.arcs[e].b[k] = data.d_k[k].quantity;
            }
        }
    }else if(style == 0){
        int pk=0;
        for(int i=0; i<narcs; ++i){
            getline(file,s);
            ss.str(s);
            ss>>data.arcs[i].i;
            ss>>data.arcs[i].j;
            ss>>data.arcs[i].f;
            ss>>data.arcs[i].capa; //std::cout<<data.arcs[i].capa<<std::endl;
            ss>>pk;
            data.grid[(data.arcs[i].i-1)*nnodes+(data.arcs[i].j-1)]=i;
            ss.clear();
            data.arcs[i].c.resize(ndemands, 0.0);
            data.arcs[i].b.resize(ndemands, 0.0);
            for(int k=0;k<ndemands;++k){
                getline(file,s);
                ss.str(s);
                ss>>pk;
                ss>>data.arcs[i].c[pk-1];
                ss>>data.arcs[i].b[pk-1];
                ss.clear();
            }
        }
        // demands origin b+, destination b-
        
        int node; double b;
        for(int i=0; i<2*ndemands; ++i){
            getline(file,s);
            ss.str(s);
            ss>> pk >> node >> b;
            if(b<0) data.d_k[pk-1].D = node;
            else{data.d_k[pk-1].O = node; data.d_k[pk-1].quantity = b;}
            ss.clear();
        }
        for(int i=0; i<ndemands; ++i){
            data.kgrid[(data.d_k[i].O-1)*nnodes+(data.d_k[i].D-1)]+= data.d_k[i].quantity;
        }

    }
    file.close();
    
}

//----------------------------------------------------------------------------------


int read_init_sol(std::string fname, std::string instance, double * xy, double& val) {
	std::ifstream file;
    file.open(fname.c_str());
    if (!file.is_open()) {
        std::cout<<"Failure to open datafile: %s\n "<<fname;
        abort();
    }
    
    std::string s;
    std::string instance_inline;
    std::istringstream ss;
    val=0;
    int cont=0;
    

    while (getline(file,s)) {
        ss.str(s);
        ss>>instance_inline;
        if(instance_inline == instance){
        	ss>>val;
        	std::cout<<std::setprecision(15)<<instance_inline<<" "<<val<<std::endl;
        	break;
        }
        ss.clear();
    }
    file.close();
    if(val==0) return -1;
    else if(xy==0) return 0;
    instance_inline = "../results/heursolutions/sol"+instance_inline.substr(instance_inline.find("/")+1); 
    std::cout<<"file: "<<instance_inline<<std::endl;
    file.open(instance_inline.c_str()); 
    if (!file.is_open()) {
        //std::cout<<"Failure to open datafile: %s\n "<<fname;
        //abort();
        std::cout<<"no solution vector"<<std::endl;
        return 0;
    }
    
    std::string aux;
	int sz, arc;
    ss.clear();
    getline(file,s);
    ss.str(s);
    
    ss>>aux;
	ss>>sz;
	ss.clear();
	//std::cout<<aux<<" "<<sz<<std::endl;
    for(;sz--;){
    	getline(file,s);
    	ss.str(s);
    	ss>>aux;
    	ss>>arc;
    	xy[arc] =1.0;
    	//std::cout<<arc<<std::endl;
    	ss.clear();
    }
    file.close();
    return 1;
}

//----------------------------------------------------------------------------------

void read_init_dualsol(std::string fname, const Data & data, double * dual){
	std::ifstream file;
    file.open(fname.c_str());
    if (!file.is_open()) {
        std::cout<<"Failure to open datafile: %s\n "<<fname;
        abort();
    }
    std::string s;
    std::string aux;
    std::istringstream ss;
     int iaux;
    for(int k=0;k<data.ndemands;++k){
    	for(int i=0;i<data.nnodes;++i){
    		getline(file,s);
    		ss.str(s);
    		ss>>aux;
    		ss>>iaux;
    		ss>>iaux;
    		ss>>dual[k*data.nnodes+i];
    		ss.clear();
    	}
    }
    int cont=0;
    while (getline(file,s)) {
        ss.str(s);
        ss>>aux;
    	ss>>iaux;
    	ss>>dual[data.ndemands*data.nnodes+cont];
    	++cont;
        ss.clear();
    }
    file.close();
    
}

//----------------------------------------------------------------------------------

bool getScalar(int id,const std::vector<PairF> & mapset, double & gam1, double & gam2){
    if(mapset[id].fst < 0 &&  mapset[id].fst != mapset[id].snd){
        std::cout<<"!!!!!! PROBLEMS !!!! getScalar"<<std::endl;
    }
    if(mapset[id].fst < 0 || mapset[id].snd < 0)
        return false;
    
    gam1 = mapset[id].fst;
    gam2 = mapset[id].snd;
    
    return true;
}

//===================================================================================
//===================================================================================
// Flow methods
//===================================================================================
//===================================================================================

void remove_i(int p, std::list<int> &N){
    
    std::list<int>::iterator it = N.begin();
    std::advance(it,p);
    N.erase(it);
}

//----------------------------------------------------------------------------------

Pair1 arg_min(const std::list<int> & nodes, const std::vector<double> & costs){
    Pair1 k;
    k.node = nodes.front();
    k.pos = 0;
    std::list<int>::const_iterator it = nodes.begin();
    for (int i=1; i<nodes.size(); ++i){
        if(costs[k.node-1]>costs[*(++it)-1]){
            k.node = *(it);
            k.pos = i;
        }
    }
    return k;
    
}

//----------------------------------------------------------------------------------


double dijkstra(int source, int sink, const std::vector<int> &grid,
                std::vector<Pair1> & preced, const std::vector<double>& cij){
    
    int nnodes= (int)preced.size();
    std::vector<double> potentials(nnodes,0.0);
    std::list<int> nodes(nnodes);
    std::list<int>::iterator it;
    Pair1 i;
    
    it = nodes.begin();
    for (int i=0; i<nnodes; ++i) {
        *(it++) = i+1;
        potentials[i] = __DBL_MAX__;
        preced[i].node = 0;
        preced[i].pos = -1;
    }
    
    potentials[source-1] = 0.0;
    while (nodes.size()>=1) {
        i=arg_min(nodes, potentials);
        remove_i(i.pos, nodes);
        if(i.node == sink) return potentials[sink-1];
        for (int j=nnodes;j--;) {
            if(grid[(i.node-1)*nnodes+j]<0) continue;
            int a = grid[(i.node-1)*nnodes+j];
            if (potentials[j]-potentials[i.node-1]-cij[a]>1e-30){
                potentials[j]=potentials[i.node-1]+cij[a];
                preced[j].node=i.node;
                preced[j].pos=a;
            }
        }
    }
    
    double ret = potentials[sink-1];
    potentials.clear();
    return ret;
}

//===================================================================================
//===================================================================================
// Knapsack dp
//===================================================================================
//===================================================================================


int MaxKnapsack(int n, int W, const int * wt, const  int *val, int * kws){
    int i, w;
    int W1 = (W+1);
    if(!kws) kws = new int [(n+1)*W1];
    
    // Build table K[][] in bottom up manner
    for (i = 0; i <= n; i++){
        for (w = 0; w <= W; w++){
            if (i==0 || w==0)
            kws[i*W1+w] = 0;
            else if (wt[i-1] <= w)
            kws[i*W1+w] = max( val[i-1] + kws[(i-1)*W1+ w-wt[i-1]],  kws[(i-1)*W1 + w] );
            else
            kws[i*W1+w] = kws[(i-1)*W1 + w];
        }
    }
    return kws[n*W1+W];
}

//===================================================================================
//===================================================================================
// Math methods
//===================================================================================
//===================================================================================

double fmin(double a, double b){ return (a < b)? a : b;  }

//----------------------------------------------------------------------------------

unsigned int min(unsigned int a, unsigned int b) { return (a < b)? a : b; }

//----------------------------------------------------------------------------------

double max(double a, double b) { return (a > b)? a : b; }

//----------------------------------------------------------------------------------

int max(int a, int b) { return (a > b)? a : b; }

//----------------------------------------------------------------------------------

int factorial(unsigned int n,  unsigned int stop){
    if (n == 0 || n==stop)
        return 1;
    return n * factorial(n - 1,stop);
}

//----------------------------------------------------------------------------------

int combination(unsigned int n, unsigned int p){
    if(p>n) return 0;
    int stop=n-p;
    int div=p;
    if(p>(n-p)){
        div = stop;
        stop = p;
    }
    //std::cout<<stop<<" "<<div<<std::endl;
    return factorial(n,stop)/factorial(div,0);
}

//----------------------------------------------------------------------------------


void  setBit(unsigned int *A,  int iA , int bit){
    unsigned int flag = 1;
    A[iA] |= (flag) << (bit);  // Set the bit at the k-th position in A[i]

}

//----------------------------------------------------------------------------------

void  clearBit(unsigned int *A,  int iA, int bit){
    A[iA] &= ~(1 << (bit));
}

//----------------------------------------------------------------------------------


int testBit(unsigned int *A,  int iA , int bit){
    return ( (A[iA] & (1 << (bit) )) != 0 ) ;
}



//----------------------------------------------------------------------------------

int countSetBits(unsigned int n){
    unsigned int count = 0;
    while (n){
        n &= (n-1) ;
        count++;
    }
    return count;
}

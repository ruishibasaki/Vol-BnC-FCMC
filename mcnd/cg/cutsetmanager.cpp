//
//  cutsetmanager.cpp
//  
//
//  Created by Rui Shibasaki on 20/05/2019.
//

#include "cutsetmanager.hpp"

#include <algorithm>

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//  initializing methods
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

void
CutSetManager::initialize(const Data * d) {
    data=d;
    sets.initialize(data->nnodes);
    //sets.print();
}

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//  main methods
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------


int
CutSetManager::cutset_generation_main(const double * ystar, const double * y, bool prepro){
    int total=0;
    std::vector<int> s(data->nnodes,0);
    for(int i=0;i<data->nnodes;++i){
        total+= compute_cutset(s, ystar, y, i, true, prepro);
        //if(total>=max)break;
        total+= compute_cutset(s, ystar, y, i, false, prepro);
        //if(total>=max)break;
        //break;
    }
    
    s.clear();
    return total;
}

//-------------------------------------------------------------------------------------------

int
CutSetManager::compute_cutset(std::vector<int>& s,  const double * ystar, const double * y, int initial, bool stosb, bool prepro){
    std::deque<int> ss_;
    std::deque<int> s_s;
    std::deque<int> ss_k;
    std::deque<int> s_sk;
    std::vector<int> saux(data->nnodes,0);
    double dss, uss, ds_s, us_s;
    int added, found;
    CutSet * newSS_;
    
    added =0;
    if(prepro) found = lin_separate_cutset(initial, ystar, saux,stosb, 0.05);
    else found = separate_cutset(initial, y, saux, stosb);
    for(int c=found;c>=1;--c){
        for(int n=data->nnodes; n--;) if(saux[n]<=c && saux[n]>0) s[n] = 1;
        if(!stosb){
            for(int n=data->nnodes; n--;) s[n] = (s[n]) ? 0 : 1;
        }
        newSS_ = sets.collect(s);
        if(newSS_!=0){
            make_cutset( s,  ss_,ss_k,  uss, dss, s_s,s_sk, us_s, ds_s);
            if(dss>0 || ds_s>0){
                added += sets.addCutSet( newSS_, ss_, ss_k,  uss,  dss,  s_s,  s_sk,  us_s,  ds_s);
            }else delete newSS_;
        }
        ss_.clear();
        s_s.clear();
        ss_k.clear();
        s_sk.clear();
        s.assign(data->nnodes,0);
        dss = uss = ds_s = us_s = 0;
    }
    return added;
}



//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//  Cutset methods
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------


void
CutSetManager::make_cutset(const std::vector<int> & s, std::deque<int> & ss_, std::deque<int> & ss_k, double &uss, double &dss, std::deque<int> & s_s, std::deque<int> & s_sk, double &us_s, double &ds_s){
    int arc;
    double k;
    dss = uss = us_s = ds_s = 0;
    //std::cout<<"make cutset: ";
    std::vector<double>css(data->narcs,0.0);
    std::vector<double>cs_s(data->narcs,0.0);
    for(int i=data->nnodes;i--;){
        if(s[i]==0) continue;
        //std::cout<<i+1<<" : "<<std::endl;
        for(int j=data->nnodes;j--;){
            arc = data->grid[i*data->nnodes+j];
            if(arc>=0 && s[j]==0){
                ss_.push_back(arc);
                uss+=data->arcs[arc].capa;
                css[arc]=1.0;
                //std::cout<<arc<<"("<<data->arcs[arc].i<<"-"<<data->arcs[arc].j<<"): "<<std::endl;
            }
           /* k = data->kgrid[i*data->nnodes+j];
            if(s[j]==0){
                //std::cout<<k<<"("<<data->d_k[k].O<<"-"<<data->d_k[k].D<<")"<<" qnt: "<<data->d_k[k].quantity<<std::endl;
                dss+=k;
            }*/
            arc = data->grid[j*data->nnodes+i];
            if(arc>=0 && s[j]==0){
                s_s.push_back(arc);
                us_s+=data->arcs[arc].capa;
                cs_s[arc]=1.0;
                //std::cout<<"("<<data->arcs[arc].i<<"-"<<data->arcs[arc].j<<"): "<<arc<<std::endl;
            }
           /* k = data->kgrid[j*data->nnodes+i];
            if(s[j]==0){
                ds_s+=k;
            }*/
        }
        //std::cout<<std::endl;
    }
    int o,d;
    double phi1, phi2;
    dss=0; ds_s=0;
    std::vector<Pair1> preced(data->nnodes, Pair1());
    for(int k =data->ndemands;k--;){
        o = data->d_k[k].O;
        d = data->d_k[k].D;
        phi1= dijkstra(o, d, data->grid, preced, css);
        dss+=phi1*data->d_k[k].quantity;
        if(phi1>0)ss_k.push_back(k);
        phi2= dijkstra(o, d, data->grid, preced, cs_s);
        ds_s+=phi2*data->d_k[k].quantity;
        if(phi2>0)s_sk.push_back(k);
        /*if(s[o-1]==1 && s[d-1]==0){
            phi1= dijkstra(o, d, data->grid, preced, css);
            //dss+=phi*data->d_k[k].quantity;
            if(phi1!=1.0) std::cout<<std::setprecision(15)<<"phi1+: "<<phi1<<" "<<o<<" "<<d<<": "<<data->d_k[k].quantity<<std::endl;
        }
        if(s[o-1]==1 && s[d-1]==1){
            phi1= dijkstra(o, d, data->grid, preced, css);
            //dss+=phi*data->d_k[k].quantity;
            //ds_s+=phi*data->d_k[k].quantity;
            if(phi1!=0.0) std::cout<<std::setprecision(15)<<"phi1o: "<<phi1<<" "<<o<<" "<<d<<": "<<data->d_k[k].quantity<<std::endl;
        }
        if(s[o-1]==0 && s[d-1]==1){
            //phi= dijkstra(o, d, data->grid, preced, cs_s);
            //ds_s+=phi*data->d_k[k].quantity;
            if(phi2!=1.0) std::cout<<std::setprecision(15)<<"phi2+: "<<phi2<<" "<<o<<" "<<d<<": "<<data->d_k[k].quantity<<std::endl;
        }
        if(s[o-1]==0 && s[d-1]==0){
            //phi= dijkstra(o, d, data->grid, preced, cs_s);
            //ds_s+=phi*data->d_k[k].quantity;
            //dss+=phi*data->d_k[k].quantity;
            if(phi2!=0.0) std::cout<<std::setprecision(15)<<"phi2o: "<<phi2<<" "<<o<<" "<<d<<": "<<data->d_k[k].quantity<<std::endl;
        }*/
    }
    preced.clear();
    css.clear();
    cs_s.clear();
    //std::cout<<"dss: "<<dss<<" ds_s:"<<ds_s<<std::endl;
}

//-------------------------------------------------------------------------------------------

int
CutSetManager::separate_cutset(int initial, const double * y, std::vector<int>& s, bool stosb){
    
    double  uOne, dss;
    double uPlus, dssPlus, delta;
    int arc,k, num;
    int jToS=-1;

    uOne=0; dss=0;
    s[initial]=1;
    for(int j=data->nnodes;j--;){
        if(stosb) arc = data->grid[initial*data->nnodes+j];
        else arc = data->grid[j*data->nnodes+initial];
        if(arc>=0){
            if(y[arc]==1)uOne+=data->arcs[arc].capa;
        }
        if(stosb) k = data->kgrid[initial*data->nnodes+j];
        else k = data->kgrid[j*data->nnodes+initial];
        dss+=k;
    }
    //std::cout<<"U1 "<<uOne<<" Dss: "<<dss<<" initial: "<<initial+1<<std::endl;
    num =1;
    delta = (uOne - dss);
    while(1){
        uPlus=0; dssPlus=0;
        jToS = get_jToS(s, y,  uPlus, dssPlus, stosb);
        if(jToS>=0){
            //std::cout<<"add "<<jToS+1<<" plus:"<<dssPlus-uPlus<<std::endl;
            uOne+=uPlus;
            dss+=dssPlus;
            delta = (uOne - dss);
            s[jToS]=num;
            if(delta<0) ++num;
        }else{ break;}
    }
    //std::cout<<"U1 "<<uOne<<" Dss: "<<dss<<" delta: "<<-delta<<std::endl;
    return num;
    
}

//-------------------------------------------------------------------------------------------

int
CutSetManager::get_jToS(const std::vector<int> & s, const double * y, double &uPlus, double&dssPlus, bool stosb){
    int jToS=-1;
    int arc;
    double k;
    double deltaPlus =0;double uPlusj;double dssPlusj;
    for(int j=data->nnodes;j--;){
        if(s[j]>0) continue;
        uPlusj=0;dssPlusj=0;
        for(int i=data->nnodes;i--;){
            if(s[i]==1){
                if(stosb){
                    k = data->kgrid[i*data->nnodes+j];
                    arc = data->grid[i*data->nnodes+j];
                }else{
                    k = data->kgrid[j*data->nnodes+i];
                    arc = data->grid[j*data->nnodes+i];
                }
                if(arc>=0)if(y[arc]==1){ uPlusj-=data->arcs[arc].capa;}
                dssPlusj-=k;
            }
        }
        if(uPlusj==0) continue;
        for(int i=data->nnodes;i--;){
            if(s[i]==0){
                if(stosb){
                    k = data->kgrid[j*data->nnodes+i];
                    arc = data->grid[j*data->nnodes+i];
                }else{
                    k = data->kgrid[i*data->nnodes+j];
                    arc = data->grid[i*data->nnodes+j];
                }
                if(arc>=0)if(y[arc]==1)uPlusj+=data->arcs[arc].capa;
                dssPlusj+=k;
            }
        }
        double deltaj =  dssPlusj - uPlusj;
        if(deltaj>=0 && deltaj>=deltaPlus){
            uPlus = uPlusj;
            dssPlus = dssPlusj;
            deltaPlus = deltaj;
            jToS = j;
        }
    }
    return jToS;
}


//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//Preprocessing methods
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

int
CutSetManager::lin_separate_cutset(int initial, const double * ystar, std::vector<int>& s, bool stosb, double tag){
    
    double yPlus;
    int  num;
    int jToS=-1;
    
    s[initial]=1;
    
    //std::cout<<" initial: "<<initial+1<<std::endl;
    num =1;
    while(1){
        yPlus=0;
        jToS = get_lin_jToS(s, ystar,  yPlus, stosb, tag);
        if(jToS>=0){
            //std::cout<<"add "<<jToS+1<<" plus:"<<yPlus<<std::endl;
            s[jToS]=num;
        }else{ break;}
    }
    //std::cout<<" num: "<<num<<std::endl;
    return num;
}

//-------------------------------------------------------------------------------------------

int
CutSetManager::get_lin_jToS(const std::vector<int> & s, const double * ystar, double &uPlus, bool stosb, double tag){
    int jToS=-1;
    int arc;
    double maxyss;
    double yPlus, ctag;
    uPlus=0;
    for(int j=data->nnodes;j--;){
        if(s[j]>0) continue;
        maxyss =0;
        for(int i=data->nnodes;i--;){
            if(s[i]==1){
                if(stosb){ arc = data->grid[i*data->nnodes+j];}
                else{arc = data->grid[j*data->nnodes+i];}
                if(arc>=0){
                    ctag = std::fabs(tag-ystar[arc]);
                    if(ctag > maxyss)
                        maxyss = ctag;
                }
            }
        }
        if(maxyss==0) continue;
        yPlus = -maxyss;
        maxyss =0;
        for(int i=data->nnodes;i--;){
            if(s[i]==0){
                if(stosb){arc = data->grid[j*data->nnodes+i];}
                else{arc = data->grid[i*data->nnodes+j];}
                if(arc>=0){
                    ctag = std::fabs(tag-ystar[arc]);
                    if(ctag > maxyss)
                        maxyss = ctag;
                }
            }
        }
        yPlus += maxyss;
        //std::cout<<"try "<<j+1<<" yss+: "<<yPlus<<std::endl;
        if(yPlus<=uPlus){
            uPlus = yPlus;
            jToS = j;
        }
    }
    return jToS;
}


#define UUana_cxx
#include "UUana.h"
#include  <vector>
const double PI = acos(-1.0);

using namespace std;

UUana::UUana(string filelist, int fr, int tr):from(fr),to(tr){

    TChain* chain  = new TChain("tt","");

    char fname[400];
    ifstream lis(filelist.c_str());
    int cnt=0;
    while(!lis.eof()){
        string filename;
        lis >> filename;
        sprintf(fname,"%s",filename.c_str());

        if(!filename.empty()) {
            if(cnt<to && cnt>=from){
                cout << fname << endl;   chain    ->Add(fname);
            }   
        }   
        cnt++;
        if(cnt>1000000) {cout<<"Too Many Files"<<endl;break;}
    }   
    fill_tree = false;
    Init(chain);
}

void UUana::InitHist(){

    sprintf(name,"random_%d_%d.root", from, to); fout = new TFile(name,"recreate");

    sprintf(name,"h_eta");     h_eta     = new TH1D(name,"",24000, -12, 12);
    sprintf(name,"h_pt");      h_pt 	 = new TH1D(name,"",1000, 0, 50);
    sprintf(name,"h_phi");     h_phi 	 = new TH1D(name,"",20000, -PI, PI);
    sprintf(name,"h_bimp");    h_bimp 	 = new TH1D(name,"",200, 0, 20);
    sprintf(name,"h_mul");     h_mul 	 = new TH1D(name,"",2000, 0, 20000 );
    sprintf(name,"h_multotal");     h_multotal 	 = new TH1D(name,"",2000, 0, 20000 );
    sprintf(name,"h_mulraw");     h_mulraw 	 = new TH1D(name,"",2000, 0, 20000 );
    sprintf(name,"h_npart");   h_npart 	 = new TH1D(name,"",500, 0-0.5, 500-0.5 );
    sprintf(name,"h_phiPT");   h_phiPT 	 = new TH2D(name,"",200, 0, 2*PI, 200, 0, 2*PI);
    sprintf(name,"h_thetaPT"); h_thetaPT = new TH2D(name,"",200, 0, 2*PI, 200, 0, 2*PI);
    sprintf(name,"h_etaphi");  h_etaphi  = new TH2D(name,"",60, -3, 3, 200, -PI, PI);
    sprintf(name,"h_spectator");  h_spectator  = new TH2D(name,"",500, 0-0.5, 500-0.5, 500, 0-0.5, 500-0.5);
    

    for(int ih=0; ih<NHAR; ih++){ sprintf(name,"h_ecc_%d", ih); h_ecc[ih] = new TH1D(name,"",11000, 0, 1.1); }

    for(int fb=0; fb<3; fb++){
        for(int ih=0; ih<NHAR-1; ih++){

            sprintf(name,"h_qn_fb_%d_%d", fb, ih);
            h_qn[fb][ih] = new TH1D(name,"",1000, -0.5, 1.1);
            sprintf(name,"q_{%d}^{obs}", ih+2); h_qn[fb][ih] ->GetXaxis()->SetTitle(name);
            sprintf(name,"h_psin_fb%d_%d", fb, ih);
            h_psin[fb][ih] = new TH1D(name,"",1000, -PI, PI);
            sprintf(name,"#Psi_{%d}^{obs}", ih+2); h_psin[fb][ih] ->GetXaxis()->SetTitle(name);
            sprintf(name,"h_qnxy_fb%d_%d", fb, ih);
            h_qnxy[fb][ih] = new TH2D(name,"",1000, -0.5, 0.5, 1000, -0.5, 0.5);
            sprintf(name,"q_{%d,x}^{obs}", ih+2); h_qnxy[fb][ih] ->GetXaxis()->SetTitle(name);
            sprintf(name,"q_{%d,y}^{obs}", ih+2); h_qnxy[fb][ih] ->GetYaxis()->SetTitle(name);

        }    
    }

    for(int ih=0; ih<NHAR-1; ih++){
        sprintf(name,"h_qndxy_%d", ih);
        h_qndxy[ih] = new TH2D(name,"",1000, -0.5, 0.5, 1000, -0.5, 0.5);
        h_qndxy[ih] ->Sumw2();
        sprintf(name,"q_{%d,x}^{B,obs}-q_{%d,x}^{F,obs}", ih+2, ih+2);
        h_qndxy[ih] ->GetXaxis()->SetTitle(name);
        sprintf(name,"q_{%d,y}^{B,obs}-q_{%d,y}^{F,obs}", ih+2, ih+2);
        h_qndxy[ih] ->GetYaxis()->SetTitle(name);
    }


    for(int ih=0; ih<NHAR-1; ih++){
        sprintf(name,"h_vns%d", ih+1);
        h_vns[ih] = new TH1D(name,"", 1000, -0.5, 1);
        sprintf(name,"h_vnc%d", ih+1);
        h_vnc[ih] = new TH1D(name,"", 1000, -0.5, 1);
    }

    for(int ih=0; ih<NHAR-1; ih++){

        sprintf(name,"hecc_psi_%d", ih);
        hecc_psi[ih] = new TH2D(name,"", 300, 0, 2*PI, 300, -PI, PI);
        sprintf(name,"#Phi_{%d}^{PP}", ih+2);   hecc_psi[ih] ->GetXaxis()->SetTitle(name);
        sprintf(name,"#Psi_{%d}^{EP,raw}", ih+2);   hecc_psi[ih] ->GetYaxis()->SetTitle(name);

        sprintf(name,"hecc_vn_ih%d", ih);
        hecc_vn[ih] = new TH2D(name,"", 500, -0.5, 1, 500,-0.5,1);
        sprintf(name,"ecc[%d]", ih+2);   hecc_vn[ih] ->GetXaxis()->SetTitle(name);
        sprintf(name,"v_{%d}", ih+2);   hecc_vn[ih] ->GetYaxis()->SetTitle(name);
    }

    outtree = new TTree("outtree","tree for debug");
    //    outtree->Branch("phibin",phibin, "phibin[200]/S");
    //    outtree->Branch("etabin",etabin, "etabin[240]/S");
    outtree->Branch("evtno",&evtno, "evtno/I");
    outtree->Branch("eff_trk",&ntrk, "eff_trk/I");
    
    outtree->Branch("thetaP",&thetaP, "thetaP/F");
    outtree->Branch("phiP",&phiP, "phiP/F");
    outtree->Branch("thetaT",&thetaT, "thetaT/F");
    outtree->Branch("phiT",&phiT, "phiT/F");
    
    outtree->Branch("eff_npart",&eff_npart, "eff_npart/I");
    outtree->Branch("npartP",&npartP, "npartP/I");
    outtree->Branch("npartT",&npartT, "npartT/I");
    
    
    outtree->Branch("mrawvnc",mrawvnc, "mrawvnc[6]/F");
    outtree->Branch("mrawvns",mrawvns, "mrawvns[6]/F");
    
    outtree->Branch("ecc",ecc, "ecc[6]/F");
    outtree->Branch("eccP",eccP, "eccP[6]/F");
    outtree->Branch("eccT",eccT, "eccT[6]/F");
    
    outtree->Branch("eccangP",eccangP, "eccangP[6]/F");
    outtree->Branch("eccangT",eccangT, "eccangT[6]/F");
    outtree->Branch("eccang",eccang, "eccang[6]/F");
    
    
    outtree->Branch("qnx",qnx, "qnx[3][6]/F");
    outtree->Branch("qny",qny, "qny[3][6]/F");
    outtree->Branch("qn",qn, "qn[3][6]/F");

    
    

}

bool UUana::calcEcc(){

    //eccentricity is calculated in recentered coordinate system

    std::vector<float> x_p;
    std::vector<float> y_p;
    std::vector<float> x_t;
    std::vector<float> y_t;
    std::vector<float> x_pt;
    std::vector<float> y_pt;
    
    std::vector<float> xp_neutron;
    std::vector<float> yp_neutron;
    std::vector<float> xt_neutron;
    std::vector<float> yt_neutron;
    std::vector<float> x_neutron;
    std::vector<float> y_neutron;
    
    

    float x_mean[3] = {0};
    float y_mean[3] = {0};
    int nparticipant[3] = {0};
    for(int i=0; i<npart; i++){
        
//        if(st_part[i]<3) continue; //inelastic collisions
        if(id_part[i]==1 ){
            x_p.push_back( x_part[i] );
            y_p.push_back( y_part[i] );
            
            x_mean[0] += x_part[i];
            y_mean[0] += y_part[i];
            nparticipant[0]++;
            
            if(pid_part[i]== 2112){
                xp_neutron.push_back(x_part[i]);
                yp_neutron.push_back(y_part[i]);
            }
        }else{
            x_t.push_back( x_part[i] );
            y_t.push_back( y_part[i] );
            
            x_mean[1] += x_part[i];
            y_mean[1] += y_part[i];
            nparticipant[1]++;
            
            if(pid_part[i]== 2112){
                xt_neutron.push_back(x_part[i]);
                yt_neutron.push_back(y_part[i]);
            }
        }
        x_pt.push_back( x_part[i] );
        y_pt.push_back( y_part[i] );
        
        if(pid_part[i]== 2112){
            x_neutron.push_back(x_part[i]);
            y_neutron.push_back(y_part[i]);
        }
        
        x_mean[2] += x_part[i];
        y_mean[2] += y_part[i];
        nparticipant[2]++;
        
        }
    
    if( nparticipant[2] == 0 ) return false;

    //cal center of mass <x>, <y>
    x_mean[0] /= nparticipant[0];
    y_mean[0] /= nparticipant[0];
    
    x_mean[1] /= nparticipant[1];
    y_mean[1] /= nparticipant[1];
    
    x_mean[2] /= nparticipant[2];
    y_mean[2] /= nparticipant[2];
    //recenter
    for(unsigned j=0; j<x_p.size(); j++){
        x_p[j] -= x_mean[0];
        y_p[j] -= y_mean[0];
    }
    
    for(unsigned j=0; j<x_t.size(); j++){
        x_t[j] -= x_mean[1];
        y_t[j] -= y_mean[1];
    }
    
    for(unsigned j=0; j<x_pt.size(); j++){
        x_pt[j] -= x_mean[2];
        y_pt[j] -= y_mean[2];
    }
    
    float rn[3][NHAR] = {0};
    float rnc_c[3][NHAR] = {0};
    float rnc_s[3][NHAR] = {0};
    
    memset(ecc, 0, sizeof ecc);
    memset(eccP, 0, sizeof eccP);
    memset(eccT, 0, sizeof eccT);
    memset(eccang, 0, sizeof eccang);
    memset(eccangP, 0, sizeof eccangP);
    memset(eccangT, 0, sizeof eccangT);


    npartP = x_p.size();
    npartT = x_t.size();
    eff_npart = x_pt.size();
    
    Np_neutron = xp_neutron.size();
    Nt_neutron = xt_neutron.size();
    N_neutron = x_neutron.size();
    

    for(unsigned j=0; j<x_p.size(); j++){
        
        for(int ih=0; ih<NHAR-1; ih++){
            int mi = ih+2;
            
            float r = sqrt(x_p[j]*x_p[j] + y_p[j]*y_p[j]);
            float phi = atan2(y_p[j], x_p[j]);
            rnc_c[0][ih] += cos((ih+2)*phi)*pow(r,mi);
            rnc_s[0][ih] += sin((ih+2)*phi)*pow(r,mi);
            rn[0][ih] += pow(r, mi);
            
        }
    }
    
    
    for(unsigned j=0; j<x_t.size(); j++){
        
        for(int ih=0; ih<NHAR-1; ih++){
            int mi = ih+2;
            
            float r = sqrt(x_t[j]*x_t[j] + y_t[j]*y_t[j]);
            float phi = atan2(y_t[j], x_t[j]);
            rnc_c[1][ih] += cos((ih+2)*phi)*pow(r,mi);
            rnc_s[1][ih] += sin((ih+2)*phi)*pow(r,mi);
            rn[1][ih] += pow(r, mi);
            
        }
    }
    
    
    for(unsigned j=0; j<x_pt.size(); j++){
        
        for(int ih=0; ih<NHAR; ih++){
            int mi = ih+2;
            
            float r = sqrt(x_pt[j]*x_pt[j] + y_pt[j]*y_pt[j]);
            float phi = atan2(y_pt[j], x_pt[j]);
            rnc_c[2][ih] += cos((ih+2)*phi)*pow(r,mi);
            rnc_s[2][ih] += sin((ih+2)*phi)*pow(r,mi);
            rn[2][ih] += pow(r, mi);
            
        }
    }
    
    
    
    for(int ih=0; ih<NHAR-1; ih++){
        
        eccP[ih] = sqrt(pow(rnc_c[0][ih],2) + pow(rnc_s[0][ih],2));
        eccP[ih] /= rn[0][ih];
        
        eccT[ih] = sqrt(pow(rnc_c[1][ih],2) + pow(rnc_s[1][ih],2));
        eccT[ih] /= rn[1][ih];
        
        ecc[ih] = sqrt(pow(rnc_c[2][ih],2) + pow(rnc_s[2][ih],2));
        ecc[ih] /= rn[2][ih];
        
        eccangP[ih] = atan2(rnc_s[0][ih], rnc_c[0][ih])+PI;
        eccangP[ih] /=(ih+2);
        eccangT[ih] = atan2(rnc_s[1][ih], rnc_c[1][ih])+PI;
        eccangT[ih] /=(ih+2);
        eccang[ih] = atan2(rnc_s[2][ih], rnc_c[2][ih])+PI;
        eccang[ih] /=(ih+2);
        
        if(ecc[ih] >0.7 ) return false;


        h_ecc[ih] ->Fill( ecc[ih] );
    }

    return true;
}

void UUana::eventInfo(){

    h_thetaPT ->Fill( thetaP, thetaT );
    h_phiPT ->Fill( phiP, phiT );
    h_bimp ->Fill( bimp );

}
/*
int UUana::getphibin(float mphi){

    double itv = 2*PI/PBIN;

    int ibin = (mphi + PI)/itv;

    if(ibin<=0) ibin=0;
    if(ibin>=PBIN) ibin = PBIN-1;

    return ibin;
} 

int UUana::getetabin(float meta){

    double itv = 2*6.0/EBIN;
    //if(fabs(meta) >6 ) cout<<"Bug "<<meta<<endl;

    int ibin = (meta + 6)/itv;

    if(ibin<=0) ibin=0;
    else if(ibin>=EBIN) ibin = EBIN-1;

    return ibin;
} 
*/
void UUana::looptrack2(){


    badtrk = 0;

//    memset(phibin,0, sizeof phibin);
//    memset(etabin,0, sizeof etabin);

    memset(qnx, 0, sizeof qnx);
    memset(qny, 0, sizeof qny);
    memset(qn0, 0, sizeof qn0);
    memset(qn, 0, sizeof qn);

    eff_trk = 0;

    memset(mrawvnc, 0, sizeof mrawvnc );
    memset(mrawvns, 0, sizeof mrawvns );
    memset(psiang, 0, sizeof psiang );


    for(int itrk=0; itrk<ntrk; itrk++){

        if(fabs(eta[itrk]) >= 6) continue;
        if(pt[itrk]<0.1) continue;

//        int ibin  = getphibin(phi[itrk]);
//        int ebin  = getetabin(eta[itrk]);

//        phibin[ibin] ++;
//        etabin[ebin] ++;

        //divide the detector into 3 eta ranges, [-6, 0], [0,6], [-6,6]
        //starts from harmonics n=2
        //calculate vn corresponding to Participate plane
        for(int ih=0; ih<NHAR-1; ih++){
            double angg =(ih+2)*phi[itrk] - (ih+2)*eccang[ih];
            mrawvnc[ih] += cos(angg);
            mrawvns[ih] += sin(angg);
        }

        if( eta[itrk] <0 ){
            for(int ih=0; ih<NHAR-1; ih++){
                qnx[0][ih] += cos( (ih+2)*phi[itrk] );
                qny[0][ih] += sin( (ih+2)*phi[itrk] );
            }
            qn0[0]++;

        }else{
            for(int ih=0; ih<NHAR-1; ih++){
                qnx[1][ih] += cos( (ih+2)*phi[itrk] );
                qny[1][ih] += sin( (ih+2)*phi[itrk] );
            }
            qn0[1]++;
        }

        eff_trk++;
    }//end of trk

    for(int ih=0; ih<NHAR-1; ih++){
        qnx[2][ih] = qnx[0][ih] + qnx[1][ih];
        qny[2][ih] = qny[0][ih] + qny[1][ih];
        psiang[ih] = atan2( qny[2][ih], qnx[2][ih]);
    }

    qn0[2] = qn0[0] + qn0[1];

    for(int fb=0; fb<3; fb++){
        if(qn0[fb] <0 ){
            cout<<"Wow you win $5 Million"<<fb<<" "<<qn0[fb]<< endl;
            flag = true;  
        } 
    }

}

/*void UUana::makeCan(){

    //make some monitoring plots

    fout->cd();

    sprintf(name,"h_qn_fullDet");
    can = new TCanvas(name,"",600*3, 600*2);
    can ->Divide(3,2);
    for(int ih=0; ih<NHAR-1; ih++){
        can->cd(ih+1);
        h_qn[2][ih]->Draw();
        gPad->SetLogy();
    }

    can->Write();

    sprintf(name,"h_Psi_fullDet");
    can = new TCanvas(name,"",600*3, 600*2);
    can ->Divide(3,2);
    for(int ih=0; ih<NHAR-1; ih++){
        can->cd(ih+1);
        h_psin[2][ih]->Draw();
    }
    can->Write();

    sprintf(name,"h_mulcan");
    can = new TCanvas(name,"",600*2, 600*2);
    can ->Divide(2,2);
    can->cd(1);
    h_mulraw->Draw();
    gPad->SetLogy();
    can->cd(2);
    h_mul->Draw();
    gPad->SetLogy();
    can->cd(3);
    h_npart->Draw();
    gPad->SetLogy();

    can->Write();
}
*/


void UUana::Loop()
{

    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntries();
    cout<<"We have "<<nentries<<" Events"<<endl;

    if(mevts>0){
        cout<<mevts<<" preset"<<endl;
        nentries = mevts;
    }

    cout<<"We will run "<<nentries<<" Events"<<endl;

    Long64_t goodevt = 0;

    for (Long64_t jentry=0; jentry<nentries;jentry++) {

        flag = false; //boolean of bad event
        eff_npart = 0;
        npartP=0;
        npartT=0;
        Nproton = 0;
        Nneutron = 0;


        if(jentry%2000==0) cout<<jentry<<endl;

        fChain->GetEntry(jentry);
        
//select tip-tip:
/*        float deltheta=30;
        deltheta/=180.0;
        deltheta*=PI;
        bool passevent=true;
        
//        static int cc = 0;
//        if(cc==0){
//            cout<<deltheta<<endl; cc = 1;
//        }
        
        
        if ((thetaP>deltheta) && (thetaP<PI-deltheta)) passevent=false;
        if ((thetaT>deltheta) && (thetaT<PI-deltheta)) passevent=false;
        
        if(!passevent) continue;
 */
 
 
//select body-body:
/*        bool passcut = true;

        float cutTa=18;
        cutTa/=180.0;
        cutTa*=PI;

        float dcutPhi=20;
        dcutPhi/=180.0;
        dcutPhi*=PI;
        
        
        if( fabs(thetaP-PI/2) > cutTa ||  fabs(thetaT-PI/2) > cutTa ) passcut = false;

        double delphiPT=phiP-phiT;
        delphiPT = sin( delphiPT );
        dcutPhi = fabs(sin(dcutPhi));
        if(  fabs(delphiPT) > dcutPhi ) passcut = false;
        
        if( !passcut) continue;
*/
        

        h_multotal->Fill( ntrk );

        if(bimp<0.001 && ntrk<7000 ){
            flag = true;
        }

        bool runecc = calcEcc();
        if(!runecc) flag = true;

        if(flag) continue;

        evtno = jentry;
        looptrack2();

        if(flag) continue;


        for(int ih=0; ih<NHAR-1; ih++){
            
            mrawvnc[ih]/=eff_trk;
            mrawvns[ih]/=eff_trk;   
	}

// v2 cut
//        bool vncut=true;
//        if(mrawvnc[0]<0.08) vncut=false;
//        if(!vncut) continue;

//ecc2 cut
//        bool ecc2cut=true;
//        if(ecc[1]<0.3) ecc2cut=false;
//        if(!ecc2cut) continue;

        
        
        
        h_mulraw->Fill( ntrk );
        
        h_spectator->Fill(NAUM-npartP,NAUM-npartT);

        for(int fb=0; fb<3; fb++){
            for(int ih=0; ih<NHAR-1; ih++){

                qnx[fb][ih] /=  qn0[fb];
                qny[fb][ih] /=  qn0[fb];

            }
        }

        for(int fb=0; fb<3; fb++){
            for(int ih=0; ih<NHAR-1; ih++){
                
                qn[fb][ih] = pow(qnx[fb][ih],2) + pow(qny[fb][ih],2);
                qn[fb][ih] = sqrt(qn[fb][ih]);
            }
        }
        
        //************************************************************************
        //begin to fill the histograms
        //************************************************************************

        for(int itrk=0; itrk<ntrk; itrk++){

            if(fabs(eta[itrk]) > 6) continue;
            if(pt[itrk]<0.1) continue;

            h_eta ->Fill( eta[itrk] );
            h_pt  ->Fill( pt[itrk] );
            h_phi ->Fill( phi[itrk] );

            h_etaphi->Fill(eta[itrk], phi[itrk]);
        }

        for(int fb=0; fb<3; fb++){
            for(int ih=0; ih<NHAR-1; ih++){

                double psin = atan2(qny[fb][ih], qnx[fb][ih]);

                h_qn[fb][ih] ->Fill( qn[fb][ih] );
                h_psin[fb][ih] ->Fill( psin );
                h_qnxy[fb][ih] ->Fill( qnx[fb][ih], qny[fb][ih] );

            }
        }

        for(int ih=0; ih<NHAR-1; ih++){
           
            h_qndxy[ih]->Fill(qnx[0][ih]-qnx[1][ih], qny[0][ih]-qny[1][ih]);
            h_vnc[ih] ->Fill( mrawvnc[ih] );
            h_vns[ih] ->Fill( mrawvns[ih] );


            hecc_vn[ih] ->Fill( ecc[ih],  mrawvnc[ih]);
        }

        for(int ih=0; ih<NHAR-1; ih++){
            hecc_psi[ih] ->Fill(eccang[ih], psiang[ih] );
        }

        h_mul ->Fill( eff_trk );
        h_npart->Fill(eff_npart);


        eventInfo();

        //************************************************************************
        //end of fill histograms
        //************************************************************************

        //if(fill_tree)
        outtree->Fill();

        goodevt++;
    }

    cout<<"We have "<<goodevt<<endl;

    fout->Write();
//    makeCan();
}

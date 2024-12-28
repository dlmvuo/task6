#include <iostream>
#include <cmath>
#include <TRandom.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLegend.h>

int N = 100000;

double E = 1020;
double E_min = 40;
double M_pi = 139.57;
double M_K = 497.61;
double tau = 2.69;

double r = 30;
double l = 50;

void pion()
{
    TRandom3 *rnd = new TRandom3();

    TH1D *Ks_Theta = new TH1D("", "", 100, 0, TMath::Pi());
    Ks_Theta->SetTitle("#theta of Ks");
    
    TH1D *Ks_Phi = new TH1D("", "", 100, 0, TMath::Pi());
    Ks_Phi->SetTitle("#phi of Ks");

    TH1D *Pi_Theta = new TH1D("", "", 100, 0, TMath::Pi());
    Pi_Theta->SetTitle("#theta of #pi_{#pm}");
    
    TH1D *Pi_Phi = new TH1D("", "", 100, 0, TMath::Pi());
    Pi_Phi->SetTitle("#phi of #pi}");
    
    TH1D *Ks_length = new TH1D("", "", 100, 0, 5);
    Ks_length->SetTitle("length Ks");

    TF1 *KsThetaDist = new TF1("KsTheta", "sin(x)*sin(x)*sin(x)", 0, TMath::Pi());

    for (int i = 0; i < N; i++) {
        
        double Ks_theta = KsThetaDist->GetRandom();
        double Ks_phi   = 2 * TMath::Pi() * rnd->Rndm() - TMath::Pi();
        double Ks_cosTheta = TMath::Cos(Ks_theta);
        double Ks_sinTheta = TMath::Sin(Ks_theta);
        double Ks_cosPhi = TMath::Cos(Ks_phi);
        double Ks_sinPhi = TMath::Sin(Ks_phi);
            
        double E_Ks = E/2;
        double p_Ks = TMath::Sqrt(pow(E_Ks,2) - pow(M_K, 2));

        TLorentzVector Ks(0, 0, 0, M_K);

        double KsBx = p_Ks*Ks_sinTheta*Ks_cosPhi / E_Ks;
        double KsBy = p_Ks*Ks_sinTheta*Ks_sinPhi / E_Ks;
        double KsBz = p_Ks*Ks_cosTheta/ E_Ks;

        double length = - TMath::Log(rnd->Rndm()) * tau * p_Ks/E_Ks; 
        TVector3 Length(length*Ks_sinTheta*Ks_cosPhi, length*Ks_sinTheta*Ks_sinPhi,length*Ks_cosTheta);

        if (TMath::Abs(Length.Z()) < l && TMath::Abs(Length.X()) < r) {
            double pi_cosTheta = 2 * rnd->Rndm() - 1;
            double pi_sinTheta = TMath::Sqrt(1 - pow(pi_cosTheta,2));

            double pi_phi   = 2 * TMath::Pi() * rnd->Rndm() - TMath::Pi();
            
            double pi_cosPhi = TMath::Cos(pi_phi);
            double pi_sinPhi = TMath::Sin(pi_phi);
            
            double E_pi = M_K/2;
            double p_pi = TMath::Sqrt(pow(E_pi,2) - pow(M_pi, 2));

            // double Pi_Px = p_pi*pi_sinTheta*pi_cosPhi;
            // double Pi_Py = p_pi*pi_sinTheta*pi_sinPhi;
            // double Pi_Pz = p_pi*pi_cosTheta;

            TLorentzVector pi_plus(p_pi*pi_sinTheta*pi_cosPhi, p_pi*pi_sinTheta*pi_sinPhi, p_pi*pi_cosTheta, E_pi);
            TLorentzVector pi_minus(-p_pi*pi_sinTheta*pi_cosPhi, -p_pi*pi_sinTheta*pi_sinPhi, -p_pi*pi_cosTheta, E_pi);

            TLorentzRotation T;
            T.Boost(KsBx, KsBy, KsBz);

            Ks = T.VectorMultiplication(Ks);
            pi_plus = T.VectorMultiplication(pi_plus);
            pi_minus = T.VectorMultiplication(pi_minus);

            if (pi_plus.P() < E_min || pi_minus.P() < E_min) continue;
           
            Ks_Theta->Fill(Ks_theta);
            Ks_Phi->Fill(Ks_phi);

            Pi_Theta->Fill(pi_plus.Theta());
            Pi_Theta->Fill(pi_minus.Theta());

            Pi_Phi->Fill(pi_plus.Phi());
            Pi_Phi->Fill(pi_minus.Phi());

            Ks_length->Fill(length);

        }
    }

    TCanvas *c = new TCanvas("canvas", "title.png", 1200, 900);
    c->Divide(3, 2);

    c->cd(1);
    Ks_Theta->SetXTitle("#theta");
    Ks_Theta->SetYTitle("N");
    Ks_Theta->Draw();
    
    c->cd(2);
    Ks_Phi->SetXTitle("#phi");
    Ks_Phi->SetYTitle("N");
    Ks_Phi->Draw();

    c->cd(3);
    Pi_Theta->SetXTitle("#theta");
    Pi_Theta->SetYTitle("N");
    Pi_Theta->Draw();

    c->cd(4);
    Pi_Phi->SetXTitle("#phi");
    Pi_Phi->SetYTitle("N");
    Pi_Phi->Draw();    
    
    c->cd(5);
    Ks_length->SetXTitle("cm");
    Ks_length->SetYTitle("N");
    Ks_length->Draw();

}

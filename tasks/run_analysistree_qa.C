/* Copyright (C) 2021 Physikalisches Institut, Eberhard Karls Universitaet Tuebingen, Tuebingen
   SPDX-License-Identifier: GPL-3.0-only
   Authors: Viktor Klochkov [committer], Viktor Klochkov [committer] */

/** @brief run_analysistree_qa
 ** @param filelist    Filefist (text file) of input AnalysisTree files
 ** @param is_single_file  if true, instead of filelist a single ROOT file will be used as input
 **
 ** Macro to run AnalysisTreeQA package (https://github.com/HeavyIonAnalysis/AnalysisTreeQA)
 ** Produces an output ROOT file with specified histograms / TProfiles.
 ** Examples how to add plots could be found here:
 ** https://github.com/HeavyIonAnalysis/AnalysisTreeQA/blob/master/examples/example.cpp
 ** To add event cuts:
 ** task->SetEventCuts(EventCuts);
 ** where EventCuts is AnalysisTree::Cuts* object, for example
 ** GetCbmEventCuts("RecEventHeader") from macro/analysis/common/cuts/cbm_cuts.h
 ** To apply a cut on some branch, for example select only primiry MC-tracks:
 ** task->AddBranchCut(GetCbmMcTracksCuts("SimParticles"));
 ** or apply quality cuts on STS tracks:
 ** task->AddBranchCut(GetCbmTrackCuts("RecTracks"));
 **/

using namespace AnalysisTree;
void EtaSearchQA(QA::Task& task, Cuts* cuts=nullptr);
void BasicTpcTracksQA(QA::Task& task, Cuts* cuts=nullptr);
void BasicMcTracksQA(QA::Task& task, Cuts* cuts=nullptr);
/*void ProtonTpcTracksQA(QA::Task& task, Cuts* cuts=nullptr);
void PionTpcTracksQA(QA::Task& task, Cuts* cuts=nullptr);
void KaonTpcTracksQA(QA::Task& task, Cuts* cuts=nullptr);*/
void EventQA(QA::Task& task);
void FHCalQA(QA::Task& task, Cuts* cuts=nullptr);
void EnergyModulesQA(QA::Task& task, Cuts* cuts=nullptr);

/*void VertexTracksQA(QA::Task& task, std::string branch=std::string("GlobalTracks"), Cuts* cuts=nullptr);
void TofHitsQA(QA::Task& task);
void SimParticlesQA(QA::Task& task, Cuts* cuts=nullptr);
void SimEventHeaderQA(QA::Task& task);
void RecEventHeaderQA(QA::Task& task);
void EfficiencyMaps(QA::Task& task);
void TrackingQA(QA::Task& task);*/

const int kPdgLambda = 10000000;
const int kPdgCharge = 10000;
const int kPdgMass   = 10;

Int_t GetIonCharge(int pdgCode) { return (pdgCode % kPdgLambda) / kPdgCharge; }

void run_analysistree_qa(std::string filelist,std::string namelist, bool is_single_file = false);

const std::string mc_tracks     = "McTracks";
const std::string fhcal_modules = "FHCalModules";
const std::string mc_event      = "McEvent";
const std::string reco_event    = "RecoEvent";
const std::string tpc_tracks    = "TpcTracks";
/*AnalysisTree::ModulePositions fhcal_modules_positions_;
fhcal_modules_positions_ = data_header_->GetModulePositions(0);*/
void run_analysistree_qa(std::string filelist,std::string namelist,bool is_single_file)
{
    if (is_single_file) {
        std::ofstream fl("fl_temp.txt");
        fl << filelist << "\n";
        fl.close();
        filelist = "fl_temp.txt";
    }

    TaskManager* man = TaskManager::GetInstance();
    auto* task = new QA::Task;
/*std::function<bool(std::vector<double>&)> lambda= [&](double n){if(n==0)
return false;};*/
    /* Cuts* event_cuts = new Cuts("event_cuts",{SimpleCut({"RecoEvent.vtx_z"}, [](std::vector<double>& var){
if(var.at(0)==0){ 
return false;
}
else
return true;
})});
    Cuts* vtx_cuts =new Cuts("vtx_cuts",{EqualsCut("McEvent.vtx_z",0)});
*/
    task->SetOutputFileName(namelist+"mpd_BiBi_9.2gev_qa.root");
    
 //   task->SetEventCuts(event_cuts);

 /*   VertexTracksQA(*task, rec_tracks);
    VertexTracksQA(*task, rec_tracks, new Cuts("Primary", {EqualsCut({sim_particles + ".mother_id"}, -1)}));
    VertexTracksQA(*task, sts_tracks);
    VertexTracksQA(*task, sts_tracks, new Cuts("Primary", {EqualsCut({sim_particles + ".mother_id"}, -1)}));
    SimParticlesQA(*task);
    SimParticlesQA(*task, new Cuts("Primary", {EqualsCut({sim_particles + ".mother_id"}, -1)}));
    SimEventHeaderQA(*task);
    RecEventHeaderQA(*task);
    TofHitsQA(*task);
    TrackingQA(*task);*/
//for(int j=0;j<90;j++)
//{
// EnergyModulesQA(*task, new Cuts(Form("Energy_in_Modules %i",j), {EqualsCut("FHCalModules.id",j)}));

//}
   // EventQA(*task);
   // FHCalQA(*task);
   // EtaSearchQA(*task, new Cuts("FHCalProtons", {EqualsCut("McTracks.pid", 2212)}));
   // EtaSearchQA(*task, new Cuts("FHCalNeutrons", {EqualsCut("McTracks.pid", 2112)}));
   // EtaSearchQA(*task, new Cuts("FHCalNuclei", {RangeCut("McTracks.pid", 1e9,1e12)}));
  //  BasicTpcTracksQA(*task, new Cuts("Secondary", {RangeCut("TpcTracks.mc_mother_id",0,1e9)}));
 //   BasicTpcTracksQA(*task);
      BasicTpcTracksQA(*task, new Cuts("Proton",{EqualsCut({"TpcTracks.mc_mother_id"},-1),EqualsCut("TpcTracks.mc_pdg",2212)}));
      BasicTpcTracksQA(*task, new Cuts("Pion",{EqualsCut({"TpcTracks.mc_mother_id"},-1),EqualsCut("TpcTracks.mc_pdg",211)}));
      BasicMcTracksQA(*task, new Cuts("McProton",{EqualsCut({"McTracks.mother_id"},-1),EqualsCut("McTracks.pid",2212)}));
      BasicMcTracksQA(*task, new Cuts("McPion",{EqualsCut({"McTracks.mother_id"},-1),EqualsCut("McTracks.pid",211)}));
/*BasicTpcTracksQA(*task, new Cuts("PrimaryProtonSearch", {EqualsCut({"TpcTracks.mc_mother_id"},-1),EqualsCut("TpcTracks.mc_pdg",2212),RangeCut("TpcTracks.eta", -1,1),RangeCut("TpcTracks.pT", 0.5,3),RangeCut("TpcTracks.tof_mass2", -0.5,0.5)}));
BasicTpcTracksQA(*task, new Cuts("PrimaryKaonSearch", {EqualsCut({"TpcTracks.mc_mother_id"},-1),EqualsCut("TpcTracks.mc_pdg",321),RangeCut("TpcTracks.eta", -1,1),RangeCut("TpcTracks.pT", 0.5,3),RangeCut("TpcTracks.tof_mass2", -0.5,0.1)}));
    BasicTpcTracksQA(*task, new Cuts("Primary", {EqualsCut({"TpcTracks.mc_mother_id"},-1),RangeCut("TpcTracks.eta", -1,1),RangeCut("TpcTracks.pT", 0.5,3)}));
    BasicTpcTracksQA(*task, new Cuts("PrimaryProton", {EqualsCut({"TpcTracks.mc_mother_id"},-1),EqualsCut("TpcTracks.mc_pdg",2212),RangeCut("TpcTracks.eta", -1,1),RangeCut("TpcTracks.pT", 0.5,3)}));
    BasicTpcTracksQA(*task, new Cuts("PrimaryKaon", {EqualsCut({"TpcTracks.mc_mother_id"},-1),EqualsCut("TpcTracks.mc_pdg",321),RangeCut("TpcTracks.eta", -1,1),RangeCut("TpcTracks.pT", 0.5,3)}));
    BasicTpcTracksQA(*task, new Cuts("PrimaryPion", {EqualsCut({"TpcTracks.mc_mother_id"},-1),EqualsCut("TpcTracks.mc_pdg",211),RangeCut("TpcTracks.eta", -1,1),RangeCut("TpcTracks.pT", 0.5,3)}));
    BasicTpcTracksQA(*task, new Cuts("TrackSelectionAllchi2", {RangeCut("TpcTracks.chi2", 0,100),RangeCut("TpcTracks.nhits", 16,100),RangeCut("TpcTracks.dca_x", -0.51,0.51),RangeCut("TpcTracks.dca_y", -0.51,0.51),RangeCut("TpcTracks.dca_z", -0.51,0.51),EqualsCut("TpcTracks.charge",1),EqualsCut("TpcTracks.tof_flag",6),RangeCut("TpcTracks.eta", -1,1),RangeCut("TpcTracks.pT", 0.5,3)}));
    BasicTpcTracksQA(*task, new Cuts("TrackSelectionAll", {RangeCut("TpcTracks.nhits", 16,100),RangeCut("TpcTracks.dca_x", -0.51,0.51),RangeCut("TpcTracks.dca_y", -0.51,0.51),RangeCut("TpcTracks.dca_z", -0.51,0.51),EqualsCut("TpcTracks.charge",1),EqualsCut("TpcTracks.tof_flag",6),RangeCut("TpcTracks.eta", -1,1),RangeCut("TpcTracks.pT", 0.5,3)}));
    BasicTpcTracksQA(*task, new Cuts("TrackSelectionPion", {RangeCut("TpcTracks.chi2", 0,100),RangeCut("TpcTracks.nhits", 16,100),RangeCut("TpcTracks.dca_x", -0.51,0.51),RangeCut("TpcTracks.dca_y", -0.51,0.51),RangeCut("TpcTracks.dca_z", -0.51,0.51),EqualsCut("TpcTracks.mc_pdg",211),EqualsCut("TpcTracks.tof_flag",6),RangeCut("TpcTracks.eta", -1,1),RangeCut("TpcTracks.pT", 0.5,3)}));
    BasicTpcTracksQA(*task, new Cuts("TrackSelectionKaon", {RangeCut("TpcTracks.chi2", 0,100),RangeCut("TpcTracks.nhits", 16,100),RangeCut("TpcTracks.dca_x", -0.51,0.51),RangeCut("TpcTracks.dca_y", -0.51,0.51),RangeCut("TpcTracks.dca_z", -0.51,0.51),EqualsCut("TpcTracks.mc_pdg",321),EqualsCut("TpcTracks.tof_flag",6),RangeCut("TpcTracks.eta", -1,1),RangeCut("TpcTracks.pT", 0.5,3)}));
    BasicTpcTracksQA(*task, new Cuts("TrackSelectionProton", {RangeCut("TpcTracks.chi2", 0,100),RangeCut("TpcTracks.nhits", 16,100),RangeCut("TpcTracks.dca_x", -0.51,0.51),RangeCut("TpcTracks.dca_y", -0.51,0.51),RangeCut("TpcTracks.dca_z", -0.51,0.51),EqualsCut("TpcTracks.mc_pdg",2212),EqualsCut("TpcTracks.tof_flag",6),RangeCut("TpcTracks.eta", -1,1),RangeCut("TpcTracks.pT", 0.5,3)}));


BasicTpcTracksQA(*task, new Cuts("TrackSelectionKaonSearch", {RangeCut("TpcTracks.chi2", 0,230),RangeCut("TpcTracks.nhits", 16,100),RangeCut("TpcTracks.dca_x", -0.51,0.51),RangeCut("TpcTracks.dca_y", -0.51,0.51),RangeCut("TpcTracks.dca_z", -0.51,0.51),EqualsCut("TpcTracks.mc_pdg",321),EqualsCut("TpcTracks.tof_flag",6),RangeCut("TpcTracks.eta", -1,1),RangeCut("TpcTracks.pT", 0.5,3),RangeCut("TpcTracks.tof_mass2", -0.5,0.1)}));
    BasicTpcTracksQA(*task, new Cuts("TrackSelectionProtonSearch", {RangeCut("TpcTracks.chi2", 0,230),RangeCut("TpcTracks.nhits", 16,100),RangeCut("TpcTracks.dca_x", -0.51,0.51),RangeCut("TpcTracks.dca_y", -0.51,0.51),RangeCut("TpcTracks.dca_z", -0.51,0.51),EqualsCut("TpcTracks.mc_pdg",2212),EqualsCut("TpcTracks.tof_flag",6),RangeCut("TpcTracks.eta", -1,1),RangeCut("TpcTracks.pT", 0.5,3),RangeCut("TpcTracks.tof_mass2", -0.5,0.5)}));*/
   /* ProtonTpcTracksQA(*task, new Cuts("ProtonTpcCuts",{RangeCut("TpcTracks.chi2", 0,230),RangeCut("TpcTracks.nhits", 16,100),EqualsCut("TpcTracks.mc_pdg",2212),RangeCut("TpcTracks.dca_x", -0.51,0.51),RangeCut("TpcTracks.dca_y", -0.51,0.51),RangeCut("TpcTracks.dca_z", -0.51,0.51)}));
    ProtonTpcTracksQA(*task, new Cuts("ProtonPrimary",{EqualsCut({"TpcTracks.mc_mother_id"},-1),EqualsCut("TpcTracks.mc_pdg",2212)}));
    PionTpcTracksQA(*task, new Cuts("PionTpcCuts",{RangeCut("TpcTracks.chi2", 0,230),RangeCut("TpcTracks.nhits", 16,100),EqualsCut("TpcTracks.mc_pdg",211),RangeCut("TpcTracks.dca_x", -0.51,0.51),RangeCut("TpcTracks.dca_y", -0.51,0.51),RangeCut("TpcTracks.dca_z", -0.51,0.51)}));
    PionTpcTracksQA(*task, new Cuts("PionPrimary",{EqualsCut({"TpcTracks.mc_mother_id"},-1),EqualsCut("TpcTracks.mc_pdg",211)}));
    KaonTpcTracksQA(*task, new Cuts("KaonTpcCuts",{RangeCut("TpcTracks.chi2", 0,230),RangeCut("TpcTracks.nhits", 16,100),EqualsCut("TpcTracks.mc_pdg",321),RangeCut("TpcTracks.dca_x", -0.51,0.51),RangeCut("TpcTracks.dca_y", -0.51,0.51),RangeCut("TpcTracks.dca_z", -0.51,0.51)}));
    KaonTpcTracksQA(*task, new Cuts("KaonPrimary",{EqualsCut({"TpcTracks.mc_mother_id"},-1),EqualsCut("TpcTracks.mc_pdg",321)}));
    ProtonTpcTracksQA(*task, new Cuts("ProtonTpcCutsmass2=0",{RangeCut("TpcTracks.chi2", 0,230),RangeCut("TpcTracks.nhits", 16,100),EqualsCut("TpcTracks.mc_pdg",2212),RangeCut("TpcTracks.dca_x", -0.51,0.51),RangeCut("TpcTracks.dca_y", -0.51,0.51),RangeCut("TpcTracks.dca_z", -0.51,0.51),EqualsCut("TpcTracks.tof_mass2",0)}));
    PionTpcTracksQA(*task, new Cuts("PionTpcCutsmass2=0",{RangeCut("TpcTracks.chi2", 0,230),RangeCut("TpcTracks.nhits", 16,100),EqualsCut("TpcTracks.mc_pdg",211),RangeCut("TpcTracks.dca_x", -0.51,0.51),RangeCut("TpcTracks.dca_y", -0.51,0.51),RangeCut("TpcTracks.dca_z", -0.51,0.51),EqualsCut("TpcTracks.tof_mass2",0)}));
    KaonTpcTracksQA(*task, new Cuts("KaonTpcCutsmass2=0",{RangeCut("TpcTracks.chi2", 0,230),RangeCut("TpcTracks.nhits", 16,100),EqualsCut("TpcTracks.mc_pdg",321),RangeCut("TpcTracks.dca_x", -0.51,0.51),RangeCut("TpcTracks.dca_y", -0.51,0.51),RangeCut("TpcTracks.dca_z", -0.51,0.51),EqualsCut("TpcTracks.tof_mass2",0)}));*/
    man->AddTask(task);
    
    man->Init({filelist}, {"aTree"});
    man->Run(-1);
    man->Finish();

    if (is_single_file) {
        // -----   Finish   -------------------------------------------------------
        std::cout << " Test passed" << std::endl;
        std::cout << " All ok " << std::endl;
    }
}
/*void FHCalQA(QA::Task& task)
   auto module_pos=fhcal_modules_positions_.GetChannel("FHCalModules.id");
    task.AddH2("x, cm", {fhcal_modules, "
}*/

void EventQA(QA::Task& task)
{
    task.AddH1({"Impact parameter, fm",{mc_event, "B"}, {100, 0, 20}});
    task.AddH1({"Reaction plane angle, rad",{mc_event, "PhiRp"},{100,-4,4}});
    task.AddH1({"Reco_vtx_x, cm", {reco_event, "vtx_x"}, {QA::gNbins, -1, 1}});
    task.AddH1({"Reco_vtx_z, cm", {reco_event, "vtx_z"}, {QA::gNbins, -100, 100}});
    task.AddH1({"Reco_vtx_y, cm", {reco_event, "vtx_y"}, {QA::gNbins, -1, 1}});
    task.AddH1({"Mc_vtx_x, cm", {mc_event, "vtx_x"}, {QA::gNbins, -1, 1}});
    task.AddH1({"Mc_vtx_z, cm", {mc_event, "vtx_z"}, {QA::gNbins, -100, 100}});
    task.AddH1({"Mc_vtx_y, cm", {mc_event, "vtx_y"}, {QA::gNbins, -1, 1}});
}
void EtaSearchQA(QA::Task& task, Cuts* cuts)
{
    task.AddH1({"pseudorapidity", {mc_tracks, "eta"}, {1000,-10,10}},cuts);
}
void EnergyModulesQA(QA::Task& task, Cuts* cuts)
{
    task.AddH1({"signal, GeV",{fhcal_modules, "signal"}, {1000, 0, 2}}, cuts);
}    
void FHCalQA(QA::Task& task, Cuts* cuts)
{

    task.AddH2({"module id, ", {fhcal_modules, "id"}, {100, 0, 100}},
               {"Phi, ", {fhcal_modules, "phi"}, {100, -8, 8}}, cuts);
    task.AddH2({"module id, ", {fhcal_modules, "id"}, {100, 0, 100}},
               {"Signal, ", {fhcal_modules, "signal_eta_signed"}, {200, -2, 2}}, cuts);
}

void BasicTpcTracksQA(QA::Task& task, Cuts* cuts)
{
 /*   task.AddH1({"DCA_{x}, cm", {tpc_tracks, "dca_x"}, {QA::gNbins, -5, 5}}, cuts);
    task.AddH1({"DCA_{y}, cm", {tpc_tracks, "dca_y"}, {QA::gNbins, -5, 5}}, cuts);
    task.AddH1({"DCA_{z}, cm", {tpc_tracks, "dca_z"}, {QA::gNbins, -5, 5}}, cuts);
    task.AddH1({"pT,   GeV/c", {tpc_tracks, "pT"},  {200, 0, 2}}, cuts);
    task.AddH1({"chi2", {tpc_tracks, "chi2"}, {500, 0, 500}}, cuts);
    task.AddH1({"pseudorapidity", {tpc_tracks, "eta"},  {200, -2, 2}}, cuts);
    task.AddH1({"number of hits", {tpc_tracks, "nhits"},  {60, 0, 60}}, cuts);
        task.AddH1({"mother_id, ", {tpc_tracks, "mc_mother_id"}, {2001, -1, 2000}}, cuts);
    task.AddH1({"probability of protons, ", {tpc_tracks, "pid_prob_proton"}, {100, 0, 1}}, cuts);
    task.AddH1({"probability of pions, ", {tpc_tracks, "pid_prob_pion"}, {100, 0, 1}}, cuts);
    task.AddH1({"probability of kaons, ", {tpc_tracks, "pid_prob_kaon"}, {100, 0, 1}}, cuts);*/
 /*   task.AddH1({"Mc particles, ", {tpc_tracks, "mc_pdg"}, {8000, -4000, 4000}}, cuts); */
/*    task.AddH2({"mass2, GeV/c^2", {tpc_tracks, "tof_mass2"}, {1000, -0.3, 2}},
               {"dE/dx, ", {tpc_tracks, "dedx"}, {1000, 0, 3000}}, cuts);
    task.AddH2({"p, GeV/c", {tpc_tracks, "p"}, {1000, 0, 5}},
               {"dE/dx, ", {tpc_tracks, "dedx"}, {1000, 0, 3000}}, cuts);
    task.AddH2({"p, GeV/c", {tpc_tracks, "p"}, {1000, 0, 5}},
               {"mass2, GeV/c^2 ", {tpc_tracks, "tof_mass2"}, {1000, -0.3, 2}}, cuts);
*/
    task.AddH2({"pseudorapidity, ", {tpc_tracks, "eta"}, {200, -2, 2}},
               {"pT, GeV/c", {tpc_tracks, "pT"}, {300, 0,3}}, cuts);

  //  task.AddH1({"tof_flag, ", {tpc_tracks, "tof_flag"}, {10, 0, 10}}, cuts);
   
}
void BasicMcTracksQA(QA::Task& task, Cuts* cuts)
{    
 //   task.AddH1({"pT,   GeV/c", {mc_tracks, "pT"},  {200, 0, 2}}, cuts);
  //  task.AddH1({"pseudorapidity", {mc_tracks, "eta"},  {200, -2, 2}}, cuts);
    task.AddH2({"pseudorapidity, ", {mc_tracks, "eta"}, {200, -2, 2}},
               {"pT, GeV/c", {mc_tracks, "pT"}, {300, 0,3}}, cuts);
}
/*void VertexTracksQA(QA::Task& task, std::string branch, Cuts* cuts)
{
    AddTrackQA(&task, branch, cuts);
    if (!sim_particles.empty()) { AddTracksMatchQA(&task, branch, sim_particles, cuts); }

    task.AddH1({"pT,   GeV/c", {tpc_tracks, "pT"},  {QA::gNbins, 0, 5}}, cuts);
    task.AddH1({"pseudorapidity", {tpc_tracks, "eta"},  {QA::gNbins, -1, 1}}, cuts);
    task.AddH1({"number of hits", {tpc_tracks, "nhits"},  {100, 0, 100}}, cuts);
    task.AddH1({"mother_id, ", {tpc_tracks, "mc_mother_id"}, {2001, -1, 2000}}, cuts);
}*/
/*void ProtonTpcTracksQA(QA::Task& task, Cuts* cuts)
{   
    task.AddH2({"mass2, GeV/c^2", {tpc_tracks, "tof_mass2"}, {1000, -0.3, 2}},
               {"dE/dx, ", {tpc_tracks, "dedx"}, {1000, 0, 3000}}, cuts);
    task.AddH2({"p, GeV/c", {tpc_tracks, "p"}, {1000, 0, 5}},
               {"dE/dx, ", {tpc_tracks, "dedx"}, {1000, 0, 3000}}, cuts);
    task.AddH2({"p, GeV/c", {tpc_tracks, "p"}, {1000, 0, 5}},
               {"mass2, GeV/c^2 ", {tpc_tracks, "tof_mass2"}, {1000, -0.3, 2}}, cuts);
      
    task.AddH2({"pseudorapidity, ", {tpc_tracks, "eta"}, {200, -1, 1}},
               {"pT, GeV/c", {tpc_tracks, "pT"}, {200, 0,2}}, cuts);   
}

void PionTpcTracksQA(QA::Task& task, Cuts* cuts)
{
    task.AddH2({"mass2, GeV/c^2", {tpc_tracks, "tof_mass2"}, {1000, -0.3, 2}},
               {"dE/dx, ", {tpc_tracks, "dedx"}, {1000, 0, 3000}}, cuts);
    task.AddH2({"p, GeV/c", {tpc_tracks, "p"}, {1000, 0, 5}},
               {"dE/dx, ", {tpc_tracks, "dedx"}, {1000, 0, 3000}}, cuts);
    task.AddH2({"p, GeV/c", {tpc_tracks, "p"}, {1000, 0, 5}},
               {"mass2, GeV/c^2 ", {tpc_tracks, "tof_mass2"}, {1000, -0.3, 2}}, cuts);

    task.AddH2({"pseudorapidity, ", {tpc_tracks, "eta"}, {200, -1, 1}},
               {"pT, GeV/c", {tpc_tracks, "pT"}, {200, 0,2}}, cuts);

}
void KaonTpcTracksQA(QA::Task& task, Cuts* cuts)
{
    task.AddH2({"mass2, GeV/c^2", {tpc_tracks, "tof_mass2"}, {1000, -0.3, 2}},
               {"dE/dx, ", {tpc_tracks, "dedx"}, {1000, 0, 3000}}, cuts);
    task.AddH2({"p, GeV/c", {tpc_tracks, "p"}, {1000, 0, 5}},
               {"dE/dx, ", {tpc_tracks, "dedx"}, {1000, 0, 3000}}, cuts);
    task.AddH2({"p, GeV/c", {tpc_tracks, "p"}, {1000, 0, 5}},
               {"mass2, GeV/c^2 ", {tpc_tracks, "tof_mass2"}, {1000, -0.3, 2}}, cuts);

    task.AddH2({"pseudorapidity, ", {tpc_tracks, "eta"}, {200, -1, 1}},
               {"pT, GeV/c", {tpc_tracks, "pT"}, {200, 0,2}}, cuts);

}*/

/*void VertexTracksQA(QA::Task& task, std::string branch, Cuts* cuts)
{
    AddTrackQA(&task, branch, cuts);
    if (!sim_particles.empty()) { AddTracksMatchQA(&task, branch, sim_particles, cuts); }

    Variable chi2_over_ndf("chi2_ndf", {{branch, "chi2"}, {branch, "ndf"}},
                           [](std::vector<double>& var) { return var.at(0) / var.at(1); });

    task.AddH1({"DCA_{x}, cm", {branch, "dcax"}, {QA::gNbins, -5, 5}}, cuts);
    task.AddH1({"DCA_{y}, cm", {branch, "dcay"}, {QA::gNbins, -5, 5}}, cuts);
    task.AddH1({"DCA_{z}, cm", {branch, "dcaz"}, {QA::gNbins, -5, 5}}, cuts);
    task.AddH1({"z_{extr}, cm", {branch, "z_first"}, {QA::gNbins, -10, 190}}, cuts);

    task.AddH1({"NDF", {branch, "ndf"}, {30, 0, 30}}, cuts);
    task.AddH1({"N_{hits}", {branch, "n_hits"}, {30, 0, 30}}, cuts);
    task.AddH1({"#chi^{2}_{vertex}", {branch, "vtx_chi2"}, {500, 0, 100}}, cuts);
    task.AddH1({"#chi^{2}/NDF", chi2_over_ndf, {QA::gNbins, 0, 100}}, cuts);
    task.AddH2({"DCA_{x}, cm", {branch, "dcax"}, {QA::gNbins, -10, 10}},
               {"DCA_{y}, cm", {branch, "dcay"}, {QA::gNbins, -10, 10}}, cuts);

    task.AddH2({"q/p first [GeV/c]", {branch, "qp_first"}, {QA::gNbins, -5, 5}},
               {"q/p last [GeV/c]", {branch, "qp_last"}, {QA::gNbins, -5, 5}}, cuts);

    task.AddH2({"t_{x} first", {branch, "tx_first"}, {QA::gNbins, -1, 1}},
               {"t_{x} last", {branch, "tx_last"}, {QA::gNbins, -1, 1}}, cuts);

    task.AddH2({"t_{y} first", {branch, "ty_first"}, {QA::gNbins, -1, 1}},
               {"t_{y} last", {branch, "ty_last"}, {QA::gNbins, -1, 1}}, cuts);
}

void TrackingQA(QA::Task& task){
    task.AddH2({"q/p STS [GeV/c]", {sts_tracks, "qp_first"}, {QA::gNbins, -5, 5}},
               {"q/p Global [GeV/c]", {rec_tracks, "qp_first"}, {QA::gNbins, -5, 5}});

    task.AddH2({"q/p STS [GeV/c]", {sts_tracks, "qp_last"}, {QA::gNbins, -5, 5}},
               {"q/p Global [GeV/c]", {rec_tracks, "qp_last"}, {QA::gNbins, -5, 5}});

    task.AddH2({"t_{x} STS", {sts_tracks, "tx_first"}, {QA::gNbins, -1, 1}},
               {"t_{x} Global", {rec_tracks, "tx_first"}, {QA::gNbins, -1, 1}});

    task.AddH2({"t_{x} STS", {sts_tracks, "tx_last"}, {QA::gNbins, -1, 1}},
               {"t_{x} Global", {rec_tracks, "tx_last"}, {QA::gNbins, -1, 1}});

    task.AddH2({"t_{y} STS", {sts_tracks, "ty_first"}, {QA::gNbins, -1, 1}},
               {"t_{y} Global", {rec_tracks, "ty_first"}, {QA::gNbins, -1, 1}});

    task.AddH2({"t_{y} STS", {sts_tracks, "ty_last"}, {QA::gNbins, -1, 1}},
               {"t_{y} Global", {rec_tracks, "ty_last"}, {QA::gNbins, -1, 1}});

    Variable sim_particles_qp("sim_particles_qp", {{sim_particles, "p"}, {sim_particles, "pid"}},
                              [](std::vector<double>& var) {
                                  auto particle = TDatabasePDG::Instance()->GetParticle(var.at(1));
                                  double charge=1.0;
                                  if( particle )
                                      charge = particle->Charge() / 3.0;
                                  else
                                      charge = GetIonCharge(var.at(1));
                                  return charge / var.at(0);
                              });

    Variable sim_particles_tx("sim_particles_tx", {{sim_particles, "px"}, {sim_particles, "pz"}, {sim_particles, "pid"}},
                              [](std::vector<double>& var) {
                                  auto particle = TDatabasePDG::Instance()->GetParticle(var.at(2));
                                  double charge=1.0;
                                  if( particle )
                                      charge = particle->Charge() > 0 ? +1.0 : -1.0;
                                  return var.at(0) / var.at(1);
                              });

    Variable sim_particles_ty("sim_particles_ty", {{sim_particles, "py"}, {sim_particles, "pz"}, {sim_particles, "pid"}},
                              [](std::vector<double>& var) {
                                  auto particle = TDatabasePDG::Instance()->GetParticle(var.at(2));
                                  double charge=1.0;
                                  if( particle )
                                      charge = particle->Charge() > 0 ? +1.0 : -1.0;
                                  return var.at(0) / var.at(1);
                              });

    task.AddH2({"q/p STS [GeV/c]", {sts_tracks, "qp_first"}, {QA::gNbins, -5, 5}},
               {"q/p GEN [GeV/c]", sim_particles_qp, {QA::gNbins, -5, 5}});

    task.AddH2({"q/p STS [GeV/c]", {sts_tracks, "qp_last"}, {QA::gNbins, -5, 5}},
               {"q/p GEN [GeV/c]", sim_particles_qp, {QA::gNbins, -5, 5}});

    task.AddH2({"t_{x} STS", {sts_tracks, "tx_first"}, {QA::gNbins, -1, 1}},
               {"t_{x} GEN", sim_particles_tx, {QA::gNbins, -1, 1}});

    task.AddH2({"t_{x} STS", {sts_tracks, "tx_last"}, {QA::gNbins, -1, 1}},
               {"t_{x} GEN", sim_particles_tx, {QA::gNbins, -1, 1}});

    task.AddH2({"t_{y} STS", {sts_tracks, "ty_first"}, {QA::gNbins, -1, 1}},
               {"t_{y} GEN", sim_particles_ty, {QA::gNbins, -1, 1}});

    task.AddH2({"t_{y} STS", {sts_tracks, "ty_last"}, {QA::gNbins, -1, 1}},
               {"t_{y} GEN", sim_particles_ty, {QA::gNbins, -1, 1}});

    task.AddH2({"q/p Global [GeV/c]", {rec_tracks, "qp_first"}, {QA::gNbins, -5, 5}},
               {"q/p GEN [GeV/c]", sim_particles_qp, {QA::gNbins, -5, 5}});

    task.AddH2({"q/p Global [GeV/c]", {rec_tracks, "qp_last"}, {QA::gNbins, -5, 5}},
               {"q/p GEN [GeV/c]", sim_particles_qp, {QA::gNbins, -5, 5}});

    task.AddH2({"t_{x} Global", {rec_tracks, "tx_first"}, {QA::gNbins, -1, 1}},
               {"t_{x} GEN", sim_particles_tx, {QA::gNbins, -1, 1}});

    task.AddH2({"t_{x} Global", {rec_tracks, "tx_last"}, {QA::gNbins, -1, 1}},
               {"t_{x} GEN", sim_particles_tx, {QA::gNbins, -1, 1}});

    task.AddH2({"t_{y} Global", {rec_tracks, "ty_first"}, {QA::gNbins, -1, 1}},
               {"t_{y} GEN", sim_particles_ty, {QA::gNbins, -1, 1}});

    task.AddH2({"t_{y} Global", {rec_tracks, "ty_last"}, {QA::gNbins, -1, 1}},
               {"t_{y} GEN", sim_particles_ty, {QA::gNbins, -1, 1}});
}

void TofHitsQA(QA::Task& task)
{
    task.AddH1({"TOF hit matching radius (cm)", {tof400_hits, "matching_radius"}, {QA::gNbins, 0, 30}});
    task.AddH1({"TOF hit x-position (cm)", {tof400_hits, "x"}, {QA::gNbins, -200, 200}});
    task.AddH1({"TOF hit y-position (cm)", {tof400_hits, "y"}, {QA::gNbins, -100, 100}});
    task.AddH1({"TOF hit z-position (cm)", {tof400_hits, "z"}, {QA::gNbins, 440, 490}});
    task.AddH1({"TOF hit time (nc)", {tof400_hits, "time"}, {QA::gNbins, 15.0, 30.0}});
    task.AddH1({"TOF hit beta (1/c)", {tof400_hits, "beta"}, {QA::gNbins, -0.2, 1.8}});
    task.AddH1({"TOF hit mass2 (GeV^{2}/c^{4})", {tof400_hits, "mass2"}, {QA::gNbins, -0.5, 4.5}});

    task.AddH2({"TOF hit x-position (cm)", {tof400_hits, "x"}, {QA::gNbins, -200, 200}},
               {"TOF hit y-position (cm)", {tof400_hits, "y"}, {QA::gNbins, -100, 100}});

    task.AddH1({"TOF hit matching radius (cm)", {tof700_hits, "matching_radius"}, {QA::gNbins, 0, 30}});
    task.AddH1({"TOF hit x-position (cm)", {tof700_hits, "x"}, {QA::gNbins, -200, 200}});
    task.AddH1({"TOF hit y-position (cm)", {tof700_hits, "y"}, {QA::gNbins, -100, 100}});
    task.AddH1({"TOF hit z-position (cm)", {tof700_hits, "z"}, {QA::gNbins, 580, 680}});
    task.AddH1({"TOF hit time (nc)", {tof700_hits, "time"}, {QA::gNbins, 0.0, 30.0}});
    task.AddH1({"TOF hit beta (1/c)", {tof700_hits, "beta"}, {QA::gNbins, -0.2, 1.8}});
    task.AddH1({"TOF hit mass2 (GeV^{2}/c^{4})", {tof700_hits, "mass2"}, {QA::gNbins, -0.5, 4.5}});

    task.AddH2({"TOF hit x-position (cm)", {tof700_hits, "x"}, {QA::gNbins, -200, 200}},
               {"TOF hit y-position (cm)", {tof700_hits, "y"}, {QA::gNbins, -100, 100}});

    Variable qp_global("qp_reco", {{rec_tracks, "charge"}, {rec_tracks, "p"}},
                       [](std::vector<double>& qp) { return qp.at(0) * qp.at(1); });

    task.AddH2({"q#timesp, GeV/c", qp_global, {QA::gNbins, -5, 5}},
               {"m^{2}, GeV^{2}/c^{4}", {tof400_hits, "mass2"}, {QA::gNbins, -1.0, 5}});
    task.AddH2({"q#timesp, GeV/c", qp_global, {QA::gNbins, -5, 5}},
               {"m^{2}, GeV^{2}/c^{4}", {tof700_hits, "mass2"}, {QA::gNbins, -1.0, 5}});

    task.AddH2({"q#timesp, GeV/c", qp_global, {QA::gNbins, -5, 5}},
               {"m^{2}, GeV^{2}/c^{4}", {tof400_hits, "beta"}, {QA::gNbins, -0.2, 1.2}});
    task.AddH2({"q#timesp, GeV/c", qp_global, {QA::gNbins, -5, 5}},
               {"m^{2}, GeV^{2}/c^{4}", {tof700_hits, "beta"}, {QA::gNbins, -0.2, 1.2}});

    task.AddH2({"Z_{last} Global", {rec_tracks, "z_last"}, {QA::gNbins, 300, 700}},
               {"Z_{TOF-400} Global", {tof400_hits, "z"}, {QA::gNbins, 300, 700}});
    task.AddH2({"Z_{last} Global", {rec_tracks, "z_last"}, {QA::gNbins, 300, 700}},
               {"Z_{TOF-700} Global", {tof700_hits, "z"}, {QA::gNbins, 300, 700}});
        return 0.5 * log((e + var[1]) / (e - var[1]));
    });

    Variable pion_y("y-y_{beam}", {{rec_tracks, "p"}, {rec_tracks, "pz"}}, [y_beam, pi_mass](std::vector<double>& var) {
        const float e = sqrt(pi_mass * pi_mass + var[0] * var[0]);
        return 0.5 * log((e + var[1]) / (e - var[1]));
    });

    Cuts* mc_protons   = new Cuts("McProtons", {EqualsCut({sim_particles + ".pid"}, 2212)});
    Cuts* mc_pions_neg = new Cuts("McPionsNeg", {EqualsCut({sim_particles + ".pid"}, -211)});
    Cuts* mc_pions_pos = new Cuts("McPionsPos", {EqualsCut({sim_particles + ".pid"}, -211)});

    task.AddH2({"#it{y}_{Lab}", {sim_particles, "rapidity"}, {QA::gNbins, -1, 5}},
               {"p_{T}, GeV/c", {sim_particles, "pT"}, {QA::gNbins, 0, 2}}, mc_protons);
    task.AddH2({"#it{y}_{Lab}", proton_y, {QA::gNbins, -1, 5}}, {"p_{T}, GeV/c", {rec_tracks, "pT"}, {QA::gNbins, 0, 2}},
               mc_protons);

    task.AddH2({"#it{y}_{Lab}", {sim_particles, "rapidity"}, {QA::gNbins, -1, 5}},
               {"p_{T}, GeV/c", {sim_particles, "pT"}, {QA::gNbins, 0, 2}}, mc_pions_neg);
    task.AddH2({"#it{y}_{Lab}", pion_y, {QA::gNbins, -1, 5}}, {"p_{T}, GeV/c", {rec_tracks, "pT"}, {QA::gNbins, 0, 2}},
               mc_pions_neg);

    task.AddH2({"#it{y}_{Lab}", {sim_particles, "rapidity"}, {QA::gNbins, -1, 5}},
               {"p_{T}, GeV/c", {sim_particles, "pT"}, {QA::gNbins, 0, 2}}, mc_pions_pos);
    task.AddH2({"#it{y}_{Lab}", pion_y, {QA::gNbins, -1, 5}}, {"p_{T}, GeV/c", {rec_tracks, "pT"}, {QA::gNbins, 0, 2}},
               mc_pions_pos);
}*/

int main(int argc, char** argv)
{
    if (argc <= 1) {
        std::cout << "Not enough arguments! Please use:" << std::endl;
        std::cout << "   ./cbm_qa filelist" << std::endl;
        return -1;
    }

    const std::string filelist = argv[1];
    const std::string namelist = argv[2];
    run_analysistree_qa(filelist,namelist);
    return 0;
}

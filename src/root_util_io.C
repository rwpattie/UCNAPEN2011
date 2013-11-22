// #ifndef __ROOT_UTILS_IO_CC__
// #define __ROOT_UTILS_IO_CC__
#include <TFile.h>
#include <TTree.h>
#include <stdio.h>
#include <iostream>

extern "C" {
  extern struct{
    float decs[89];
    float costheta[12];
  } hbooku_;
}

extern "C" {
  extern struct{
    double trgeast[10];
    double trgwest[10];
  } trgt_;
}

extern "C" {
 extern struct{
   double tracktime;
   double abor;
   double ae;
   double ax;
   double ay;
   double az;
   double au;
   double av;
   double aw;
 } abort_;
}

extern "C" {
 extern struct{
  double epr;
  double xpr;
  double ypr;
  double zpr;
  double upr;
  double vpr;
  double wpr;
  double ptof;
 }proton_;
}

extern "C" {
 extern struct{
   double phte;
   double phtw;
   double phten;
   double phtwn;
   double wi;
   double wim;
   double time;
 } npebl_;
}

extern "C" {
 extern struct{
   int npar0;
   int kpar0;
   double ei;
   double x0;
   double y0;
   double z0;
   double u0;
   double v0;
   double w0;
   double ptype0;
 } init_;
}

extern "C" {
  extern struct{
    double epe,epw,eper,eped,epwd,epwr,egea;
    double egwa,egedf,egedb,egwdf,egwdb,efle;
    double eflw,etube,eal,eme1,eme2,emw1,emw2;
    double phw,phen,phe,phwn,eholder,eean,eec1;
    double eec2,ewan,ewc1,ewc2,ecol;
  }edep_;
}

extern "C" {
  extern struct{
    double emx,emy,wmx,wmy;
    double eposx,eposy,wposx,wposy;
  }multi_;
}

extern "C"{
  extern struct{
    double w[12]; 
  }cos_;
}


TFile *fout;
TTree *tree;

extern "C" int openrootfile_(char *hbookfile,int ll)
{
    hbookfile[ll--] = '\0';  // ADD A NULL CHARACTER TO THE END
    printf("ROOT File to open is : %s\n",hbookfile);
    fout = new TFile(hbookfile,"RECREATE");
    tree = new TTree("tree","Output ntuple from PENELOPE");
    tree->Branch("init",&init_,"Npar/I:Kpar:E/D:X:Y:Z:U:V:W:Ptype");
    tree->Branch("edep",&edep_,"Epe/D:Epw:Eper:Eped:Epwd:Epwr:Egea:Egwa:"
			       "Egedf:Egedb:Egwdf:Egwdb:Efle:Eflw:"
			       "Etube:Eal:Eme1:Eme2:Emw1:Emw2:Phw:Phen:"
			       "Phe:Phwn:Eholder:Eean:Eec1:Eec2:Ewan:"
			       "Ewc1:Ewc2:Ecol");
    tree->Branch("trgt"  ,&trgt_,"trgeast[10]/D:trgwest[10]");
    tree->Branch("proton",&proton_,"Ep/D:Xp:Yp:Zp:Up:Vp:Wp:PTOF");
    tree->Branch("cos"   ,&cos_,"W1[12]/D");
    tree->Branch("multi" ,&multi_,"Emx/D:Emy:Wmx:Wmy:Eposx:Eposy:Wposx:Wposy");
    tree->Branch("abort" ,&abort_,"Time/D:Abor:Ae:Ax:Ay:Az:Au:Av:Aw");
    
    return(1);
}

extern "C" void filltree_()
{
  tree->Fill(); 
}

extern "C" void closerootfile_()
{
  fout->Write();
  fout->Close();
}

// #endif
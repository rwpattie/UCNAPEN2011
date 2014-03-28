// #ifndef __ROOT_UTILS_IO_CC__
// #define __ROOT_UTILS_IO_CC__
// #include <Math/Random.h>
#include <TFile.h>
#include <TTree.h>
#include <stdio.h>
#include <iostream>

extern "C" {
  extern struct{
    float decs[96];
    float costheta[16];
  } hbooku_;
}

extern "C" {
  extern struct{
    double trge[10];
    double trgw[10];
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
    double eec2,ewan,ewc1,ewc2,ecol,eapd,eapdh;
  }edep_;
}

extern "C" {
  extern struct{
    double emx,emy,wmx,wmy;
    double eposx,eposy,wposx,wposy;
    double exsci,eysci,wxsci,wysci;
    double tapd;
  }multi_;
}

extern "C"{
  extern struct{
    double w1[16]; 
  }cos_;
}


TFile *fout;
TTree *tree;

char stuff[80];
char junk[80];
char edood[80];
char wdood[80];
char edude[80];
char wdude[80];

extern "C" int openrootfile_(char *hbookfile,int ll)
{
    hbookfile[ll--] = '\0';  // ADD A NULL CHARACTER TO THE END
    printf("ROOT File to open is : %s\n",hbookfile);
    fout = new TFile(hbookfile,"RECREATE");
    tree = new TTree("tree","Output ntuple from PENELOPE");
    tree->Branch("Npar",&init_.npar0,"Npar/I");
    tree->Branch("Kpar",&init_.kpar0,"Kpar/I");
    tree->Branch("E",&init_.ei,"E/D");
    tree->Branch("X",&init_.x0,"X/D");
    tree->Branch("Y",&init_.y0,"Y/D");
    tree->Branch("Z",&init_.z0,"Z/D");
    tree->Branch("U",&init_.u0,"U/D");
    tree->Branch("V",&init_.v0,"V/D");
    tree->Branch("W",&init_.w0,"W/D");
    tree->Branch("Ptype",&init_.ptype0,"Ptype/D");
    tree->Branch("Epe",&edep_.epe,"Epe/D");
    tree->Branch("Epw",&edep_.epw,"Epw/D");
    tree->Branch("Eper",&edep_.eper,"Eper/D");
    tree->Branch("Eped",&edep_.eped,"Eped/D");
    tree->Branch("Epwr",&edep_.eper,"Epwr/D");
    tree->Branch("Epwd",&edep_.eped,"Epwd/D");
    tree->Branch("Egea",&edep_.egea,"Egea/D");
    tree->Branch("Egwa",&edep_.egwa,"Egwa/D");
    tree->Branch("Egedf",&edep_.egedf,"Egedf/D");
    tree->Branch("Egedb",&edep_.egedb,"Egedb/D");
    tree->Branch("Egwdf",&edep_.egwdf,"Egwdf/D");
    tree->Branch("Egwdb",&edep_.egwdb,"Egwdb/D");
    tree->Branch("Efle",&edep_.efle,"Efle/D");
    tree->Branch("Eflw",&edep_.eflw,"Eflw/D");
    tree->Branch("Etube",&edep_.etube,"Etube/D");
    tree->Branch("Eal",&edep_.eal,"Eal/D");
    tree->Branch("Eme1",&edep_.eme1,"Eme1/D");
    tree->Branch("Eme2",&edep_.eme2,"Eme2/D");
    tree->Branch("Emw1",&edep_.emw1,"Emw1/D");
    tree->Branch("Emw2",&edep_.emw2,"Emw2/D");
    tree->Branch("Phe",&edep_.phe,"Phe/D");
    tree->Branch("Phen",&edep_.phen,"Phen/D");
    tree->Branch("Phw",&edep_.phw,"Phw/D");
    tree->Branch("Phwn",&edep_.phwn,"Phwn/D");
    tree->Branch("Eean",&edep_.eean,"Eean/D");
    tree->Branch("Ewan",&edep_.ewan,"Ewan/D");
    tree->Branch("Eec1",&edep_.eec1,"Eec1/D");
    tree->Branch("Eec2",&edep_.eec2,"Eec2/D");
    tree->Branch("Ewc1",&edep_.ewc1,"Ewc1/D");
    tree->Branch("Ewc2",&edep_.ewc2,"Ewc2/D");
    tree->Branch("Eholder",&edep_.eholder,"Eholder/D");
    tree->Branch("Ecol",&edep_.ecol,"Ecol/D");
    tree->Branch("Eapd",&edep_.eapd,"Eapd/D");
    tree->Branch("Eapdh",&edep_.eapdh,"Eapdh/D");
    for (int i=1;i<11;i++){
        sprintf(edood,"trge%d",i);
        sprintf(edude,"trge%d/D",i);
        tree->Branch(edood, &trgt_.trge[i-1], edude);
        sprintf(wdood,"trgw%d",i);
        sprintf(wdude,"trgw%d/D",i);
        tree->Branch(wdood, &trgt_.trgw[i-1], wdude);
    }
    tree->Branch("Ep",&proton_.epr,"Ep/D");
    tree->Branch("Xp",&proton_.xpr,"Xp/D");
    tree->Branch("Yp",&proton_.ypr,"Yp/D");
    tree->Branch("Zp",&proton_.zpr,"Zp/D");
    tree->Branch("Up",&proton_.upr,"Up/D");
    tree->Branch("Vp",&proton_.vpr,"Vp/D");
    tree->Branch("Wp",&proton_.wpr,"Wp/D");
    tree->Branch("PTOF",&proton_.ptof,"PTOF/D");
    for (int i=1;i<17;i++){
        sprintf(stuff,"W%d",i);
        sprintf(junk,"W%d/D",i);
        tree->Branch(stuff, &cos_.w1[i-1], junk);
    }
    tree->Branch("emx",&multi_.emx,"emx/D");
    tree->Branch("emy",&multi_.emy,"emy/D");
    tree->Branch("wmx",&multi_.wmx,"wmx/D");
    tree->Branch("wmy",&multi_.wmy,"wmy/D");
    tree->Branch("eposx",&multi_.eposx,"eposx/D");
    tree->Branch("eposy",&multi_.eposy,"eposy/D");
    tree->Branch("wposx",&multi_.wposx,"wposx/D");
    tree->Branch("wposy",&multi_.wposy,"wposy/D");
    tree->Branch("exsci",&multi_.exsci,"exsci/D");
    tree->Branch("eysci",&multi_.eysci,"eysci/D");
    tree->Branch("wxsci",&multi_.wxsci,"wxsci/D");
    tree->Branch("wysci",&multi_.wysci,"wysci/D");
    tree->Branch("Tapd",&multi_.tapd,"Tapd/D");
    tree->Branch("Time",&abort_.tracktime,"Time/D");
    tree->Branch("Abor",&abort_.abor,"Abor/D");     
    tree->Branch("Ae",&abort_.ae,"Ae/D");    
    tree->Branch("Ax",&abort_.ax,"Ax/D");     
    tree->Branch("Ay",&abort_.ay,"Ay/D");
    tree->Branch("Az",&abort_.az,"Az/D");
    tree->Branch("Au",&abort_.au,"Au/D");
    tree->Branch("Av",&abort_.av,"Av/D");
    tree->Branch("Aw",&abort_.aw,"Aw/D");

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

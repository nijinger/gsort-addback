TCutG *load_masslcut(){
//========= Macro generated from object: masslcut/Graph
//========= by ROOT version5.34/36
   
   TCutG *cutg = new TCutG("masslcut",18);
   cutg->SetVarX("fthetaL[0]");
   cutg->SetVarY("dT[0]");
   cutg->SetTitle("Graph");
   cutg->SetFillColor(1);
   cutg->SetPoint(0,0.378832,-120.649);
   cutg->SetPoint(1,0.463686,-78.7816);
   cutg->SetPoint(2,0.586861,-41.4003);
   cutg->SetPoint(3,0.74927,-8.50475);
   cutg->SetPoint(4,0.950912,7.94304);
   cutg->SetPoint(5,1.2365,48.3149);
   cutg->SetPoint(6,1.29854,54.2959);
   cutg->SetPoint(7,1.29945,24.3908);
   cutg->SetPoint(8,1.12701,-12.9905);
   cutg->SetPoint(9,0.943613,-48.8766);
   cutg->SetPoint(10,0.785766,-108.687);
   cutg->SetPoint(11,0.670803,-190.926);
   cutg->SetPoint(12,0.573175,-297.089);
   cutg->SetPoint(13,0.414416,-397.271);
   cutg->SetPoint(14,0.362409,-350.918);
   cutg->SetPoint(15,0.341423,-222.326);
   cutg->SetPoint(16,0.341423,-222.326);
   cutg->SetPoint(17,0.378832,-120.649);
   return cutg;
}

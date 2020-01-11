TCutG *load_massrcut(){
//========= Macro generated from object: massrcut/Graph
//========= by ROOT version5.34/36
   
   TCutG *cutg = new TCutG("massrcut",23);
   cutg->SetVarX("fthetaR[0]");
   cutg->SetVarY("dT[0]");
   cutg->SetTitle("Graph");
   cutg->SetFillColor(1);
   cutg->SetPoint(0,0.352235,259.146);
   cutg->SetPoint(1,0.368072,175.411);
   cutg->SetPoint(2,0.43054,117.097);
   cutg->SetPoint(3,0.519402,72.2389);
   cutg->SetPoint(4,0.619702,46.8196);
   cutg->SetPoint(5,0.734078,28.8766);
   cutg->SetPoint(6,0.863412,12.4288);
   cutg->SetPoint(7,0.983948,1.96203);
   cutg->SetPoint(8,1.12032,-20.4668);
   cutg->SetPoint(9,1.23998,-44.3908);
   cutg->SetPoint(10,1.29365,-66.8196);
   cutg->SetPoint(11,1.30244,-42.8956);
   cutg->SetPoint(12,1.23294,-8.50475);
   cutg->SetPoint(13,1.13088,18.4098);
   cutg->SetPoint(14,1.01034,31.8671);
   cutg->SetPoint(15,0.897725,61.7722);
   cutg->SetPoint(16,0.781589,99.1535);
   cutg->SetPoint(17,0.688328,173.916);
   cutg->SetPoint(18,0.604745,257.65);
   cutg->SetPoint(19,0.522041,380.261);
   cutg->SetPoint(20,0.469252,471.472);
   cutg->SetPoint(21,0.382149,456.519);
   cutg->SetPoint(22,0.352235,259.146);
   return cutg;
}
